! Subrutine for calculating the local values of the viscosity, the density and the dimentionles parameter EPS
    ! Different sets of equations are solved depending on the value of the lub_param parameter
    ! See more about the different options in the lubrication_def subroutine where the parameters are initialized
	SUBROUTINE Newtonian(SS, NYs, T40, Tcvot, NN, t, k, M_conv)
        implicit none 
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_Contact_mat.h'
        
        include     'inc_EDA_contact.h'
        include     'inc_Grid.h'
        include     'inc_G0DT.h'
        include     'inc_Holmes.h'
        
	    include     'inc_Method.h'
        include     'inc_Outp.h'
        include     'inc_Ref.h'
        include     'inc_Rho.h'
        include     'inc_RLarsson.h'
        
        include     'inc_Shear_lim.h'
	    include     'inc_Temp_param.h'
        include     'inc_Temp_reduction.h'
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'
	    include     'inc_Y_liu.h'
        include     'inc_Cores_used.h'
        ! Input
        integer     NN, NYs, SS, t, k, M_conv
        real        T40, Tcvot
        ! Calculations
        integer     I,J, JJ, temp_iter, err, printer
        real        EDA_1, EDA_2, EDA_3
        real        EDA22, EDA33, EDA1
        real        pg, Tg, YF 
        real        temp_0, u_slip, h_real, tau
        real        tau_x_f, tau_y_f, tau_x_c, tau_y_c, tau_lim_c
        real        dPdx, dPdy, dpds
        integer     J0, J1, J2, I0, I1, I2
        real        eda01,eda02,eda03,eda04,delta_s
        real        tau_gamma_r, dxss_Phb
        real        epst
        logical     temp_1_1_initialized
        ! Output
        SAVE      /Current/                    ! Current timestep               ! ??? If save is not called will the changes to the common block variables not be registered? 
        SAVE      /CurrentT/
        
    ! Calculate the viscocity, the density and the dimentionles parameter EPS -----------------------------------------------------------------------------------------
    IF( Lub_param .EQ. 1) THEN                                         !D_IN: /Ref/ -> Lub_param
        ! Newtonian acc X.Tans paper ref 30 31, Cylinder
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,EDAx,alpha,Pref,Z,P,PH,EDAy,RO,RA1,RA2,xi,contact,EDA_cont,EPSx,H,ENDA,EPSy,EDA0,EpsT,t) &
        !$OMP&            PRIVATE(J,I)
        DO J=1,NN,SS
            DO I=1,NX,SS                                               !D_IN: /Grid/ ->Nx

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(alpha*Pref/Z*(-1+(1+P(I,J)*PH/Pref)**Z))      !D_IN: /Visc/ -> alpha, PH, Pref, Z; /CurrentP/ -> P(I,J) !D_OUT:  EDAx(I,J) -> /Current/
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to some reasonable? values. Intended to increase convergense. 
                EDAy(I,J) = EDAx(I,J)                                     !D_OUT: EDAy(I,J) -> /Current/
                
            
                ! D-H Formulation acc X.Tan
                 !RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J))          !D_IN: /Rho/ -> RA1, RA2; !D_OUT: RO -> /Current/
                 xi(I,J)=0.0                                        !D_OUT: xi(I,J) -> /Current/     ! Should be able to remove xi becouse we're never changing it
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    !$OMP CRITICAL (C_3)
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J),  'EDAO = ',EDA0                  !D_IN: /Visc/ -> EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)             !CALL: current time step t
                    !$OMP END CRITICAL (C_3)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    !$OMP CRITICAL (C_4)
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT       !??? EpsT is not initialized anywhere
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                    !$OMP END CRITICAL (C_4)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity            !D_IN: /Contact_mat/ -> Contact(I,J); /EDA_contact/ -> EDA_cont
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))                                !D_IN: /Visc/ -> ENDA; /CurrentH/ -> H(I,J); !D_OUT: Epsx -> /Current/ 
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))                      !D_OUT: Epsy -> /Current/
                
                
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
    ELSE IF( Lub_param .EQ. 11) THEN
        ! Barus equation with extra Z
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,EDAx,alpha,Pref,Z,P,PH,EDAy,RO,RA1,RA2,xi,contact,EDA_cont,EPSx,H,ENDA,EPSy,EDA0,EpsT,t) &
        !$OMP&            PRIVATE(J,I) 
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(alpha*(P(I,J)*PH)**Z)
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
            
                ! D-H Formulation acc X.Tan
                 !RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
            ENDDO
        ENDDO 
        !$OMP END PARALLEL DO
    ELSEIF( Lub_param .EQ. 2) THEN
        ! Newtonian acc Gohars book Elastohydrodynamics
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,PH,Z,H,P,EDAx,EDAy,RA1,RA2,RO,xi,EDA0,contact,EDA_cont,ENDA,t,EPSx,EPSy,EpsT) &
        !$OMP&            PRIVATE(I,J,EDA_1,EDA_2,EDA_3) 
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc Gohar
                EDA_1=5.1*PH*P(I,J)/10**(9)
                EDA_2=(1.0+EDA_1)**Z
                EDA_3=(log(EDA0)+9.67)
                EDAx(I,J)=EXP(EDA_3*(EDA_2-1.0))
                
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
            
                ! D-H Formulation acc X.Tan
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    !$OMP CRITICAL (C_5)
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                    !$OMP END CRITICAL (C_5)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    !$OMP CRITICAL (C_6)
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                    !$OMP END CRITICAL (C_6)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                ENDDO 
        ENDDO
        !$OMP END PARALLEL DO
        
    ELSE IF( Lub_param .EQ. 3)THEN
        ! Newtonian input parameters ac Holmes et al Transient EHL point contact analysis 2003, Ball 
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,PH,Z,H,P,EDAx,EDAy,RA1,RA2,RO,xi,EDA0,contact,EDA_cont,ENDA,t,EPSx,EPSy,EpsT,kH,xH) &
        !$OMP&            PRIVATE(I,J) 
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc Holmes
                 EDAx(I,J)=EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1)) 
                 EDAy(I,J)=EDAx(I,J) !EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1)) 
                
                ! D-H Formulation acc Homes et al
                 !RO(I,J)=(1 + gH*PH*P(I,J))/(1+lH*P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
            ENDDO
    
        ENDDO
        !$OMP END PARALLEL DO
        
    ELSE IF( Lub_param .EQ. 4 .and. lub_temp .GE. 0)THEN      ! No temperature adjustments
        ! R.Larssons 2000 formulation
        ! This part is copied to EDA_calc
        if ( temp_param == 2 ) then                           !D_IN: /temp_reduction/ -> temp_param
            ub= 2*um-ua                                       !D_IN: /Outp/ -> Um, Ua; !D_OUT: Ub -> /Outp/   (Not saved)
            u_slip=ua-ub
        endif
        printer=0
        dxss_Phb=1.0/(2*DX*SS)*Ph/b                           !D_IN: /Grid/ -> Dx; /Outp/ -> b
        
        temp_1_1_initialized = .false.
        
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,Temp_param,temp,temp_fac,RL_G0,S0,Dz,Cz,P,PH,U_slip,EDAx,p_red,temp_red,H,b,Rx,NYs,term1,term2,term3,term4,dxss_Phb,tau_gamma,tau_lim,lim_factor,t,contact,RO,EpsT0,RL_c,RA1,RA2,RL_T0,EDA_cont,EPSx,EPSy,ENDA,temp_1_1_initialized) &
        !$OMP&            PRIVATE(I,J,err,RL_Ta,EDA0,EDA1,h_real,J0,J2,J1,JJ,I0,I2,I1,dPdx,tau,dpdy,tau_gamma_r,tau_lim_c,dpds,EpsT,Z)                                                               &
        !$OMP&            REDUCTION(+:printer)
        DO J=1,NN,SS
            DO I=1,NX,SS
899             err=0                                         ! !!! Branch target 899
                If ( Temp_param == 3) then
                    IF((I .NE. 1) .AND. (J .NE. 1) .AND. (.NOT. temp_1_1_initialized)) then   !Assuming the 1,1 iteration is the first iteration assigned to one of the threads. This is correct when using the STATIC schedule (which we do here since its the default)
                        DO WHILE (.NOT. temp_1_1_initialized)
                        !$OMP FLUSH (temp_1_1_initialized) 
                        ENDDO
                    ENDIF
                    RL_Ta=temp(1,1) +( temp(I,J) - temp(1,1))*temp_fac      !D_IN: /CurrentT/ -> temp; /Temp_reduction/ -> temp_fac; !D_OUT: RL_ta -> /RLarsson/  (Not saved)
                    IF( RL_Ta .LT. temp(1,1)) then
                        RL_ta=temp(1,1)
                    endif
                else
                    RL_Ta=temp(I,J)
                endif
                
                ! Roelands equation acc R.Larsson
                EDA0 = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)             ! Eq (3)         !D_IN: /RLarsson/ -> RL_G0, S0   !D_OUT: EDA0 -> /Visc/  (Not saved)
                Z=Dz+Cz*log10(1.0+RL_Ta/135)                                               !D_IN: /RLarsson/ -> Dz, Cz      !D_OUT: Z -> /Visc/     (Not saved)
                
                EDA1=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                
                IF( U_slip == 0) then
                    EDAx(i,j) = EDA1
                    
                elseif ( Temp_param == 1) then
                    ! Reducing eta
                    EDAx(i,j) = EDA1 *  1/ (1 + temp_fac * max( 0.0, (P(i,j)-p_red)/p_red) * ( max(0.0 , (RL_Ta-temp_red)/temp_red ) ) )     !D_IN: /Temp_reduction/ -> p_red, temp_red
                    
                elseif( temp_param == 2 .or. temp_param==6) then
                    ! applying a shear limmit  based on the pressure according to data from Höglund 1989 and 1999 with the temperature reduction of the shear limit
                    h_real=H(I,J)*b**2/Rx                            !D_IN: /outp/ -> Rx
                    
                    ! Define nodenumbers for past and next nodes
                    J0=J-SS
                    IF( J0 .LE. 0)  J0=J+SS
                    J2=J-2*SS
                    IF( J2 .LE. 0)  J2=J0
                    J1=J+1*SS
                    IF( J1 .GE. NYs) J1=NYs-SS
                    JJ=NYs+1-J
                    IF( JJ .LE. 1)  JJ=1+SS
                
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I2=I-2*SS
                    IF( I2 .LE. 0)  I2=I0
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                    
                    dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))* dxss_Phb ! /(2*DX*SS)*Ph/b      !D_IN: /Method/ -> term1, term2, term3, term4
                    tau = abs(H_real/2*dpdx) + abs(EDA1*U_slip/H_real) 
                    
                    dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))* dxss_Phb 
                    tau= sqrt( tau**2 + (H_real/2*dpdy)**2)
                    
                    if (temp_red .GT. 0) then                                       ! Take calculate temperature reductio of the shear limit
                        tau_gamma_r = tau_gamma - temp_fac * RL_TA     ! tau_gamma and temp_fac anready includes PH        !D_IN: /Temp_reduction/ -> tau_gamma
                        if( tau_gamma_r .LT. temp_red) then                             ! Just to make sure that not negative
                            tau_gamma_r = temp_red
                        endif
                        
                    else
                        tau_gamma_r = tau_gamma
                    endif
                    
                    if (temp_param==2)then
                        tau_lim_c = tau_lim + tau_gamma_r * P(i,j)                  !D_IN: /Temp_reduction/ -> tau_lim
                    elseif ( temp_param == 6) then
                        tau_lim_c = tau_lim + tau_gamma_r * P(i,j) + tau*p_red
                    endif
                    
                    if (tau .gt. tau_lim_c) then
                        dpds= sqrt( dPdx**2 + dPdy**2)
                        
                        if( tau_lim_c - abs(H_real/2*dpds) .LT. lim_factor*tau_lim ) then       !D_IN: /Temp_reduction/ -> lim_factor
                            
                            ! To increase numerical stability the maximum pressure of the neighbouring nodes is used
                            if (temp_param==2)then
                                tau_lim_c = tau_lim + tau_gamma_r *max( P(i0,j),P(i,j),P(i1,j), P(i,j0), P(i,j1))
                            elseif ( temp_param == 6) then
                                tau_lim_c = tau_lim + tau_gamma_r *max( P(i0,j),P(i,j),P(i1,j), P(i,j0), P(i,j1)) + tau*p_red
                            endif
                            
                            if( tau_lim_c - abs(H_real/2*dpds) .LT. lim_factor*tau_lim ) then
                                EDA1 = abs( (lim_factor*tau_lim)*H_real/U_slip)    ! Adding a factor of 2 here in hope of increased numerical stability
                                printer=printer+1
                            else
                                EDA1 = abs( (tau_lim_c - abs(H_real/2*dpds))*H_real/U_slip)
                            endif
                            
                        else
                            EDA1 = abs( (tau_lim_c - abs(H_real/2*dpds))*H_real/U_slip)
                        endif
                    endif
                    
                    EDAx(i,j) = EDA1
                elseif( temp_param== 4 .or. temp_param== 5) then ! affects the temp_calc_ij subroutine
                    EDAx(i,j) = EDA1
                else
                    EDAx(i,j) = EDA1
                endif
                
                ! D-H Formulation acc R.Larsson
                EpsT       = EpsT0*exp(-RL_c*P(I,J)*PH)               !D_IN: /RLarsson/ -> EpsT0, RL_c
                RO(I,J)    = (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0))      !D_IN: /RLarsson/ -> RL_T0
                !xi(I,J)    = 0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) .or. (EDAx(i,j) .gt. 1e20) ) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_6)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_6)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 899                           ! !!! Branch to 899
                    else
                        !$OMP CRITICAL (C_7)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_7)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_8)
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_8)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 899
                    else
                        !$OMP CRITICAL (C_9)
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_9)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
                
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                IF ((I .EQ. 1) .AND. (J .EQ. 1)) then   
                    temp_1_1_initialized = .true.
                    !$OMP FLUSH (temp_1_1_initialized)   
                ENDIF
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO

     if (printer .gt. 200 .and. k == 1) then
            WRITE(4,*)'Bad equation for shear limmit ', printer ,' number of times'
            Call Stop_to_large_out(t)
            if ( printer .gt. 1800 .and. SS == 1 .and. M_conv .gt. 10) then  ! This is a clear sign that the code is diverging. Howevere, at the first itterations of a timestep it might not be true. 
                WRITE(4,*) 'Breaking due to many bad shear limit equations'
                WRITE(4,*) 'The time is t=',t
                OPEN(20,FILE='STOP_Too_many_bad_Shear.DAT',STATUS='UNKNOWN')
                stop
            endif
            
            
    endif
    
    ELSE IF(lub_param .EQ. 4 .and. lub_temp == -1) then ! simple Themperature rise     
        ub= 2*um-ua
        u_slip=ua-ub
        temp_1_1_initialized = .false.

        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,NYs,DX,Rx,temp_1_1_initialized,temp,RL_G0,S0,Dz,Cz,EDAx,EDAy,PH,b,P,EpsT0,RL_c,RO,RA1,RA2,RL_T0,xi,H,contact,EDA_cont,term1,term2,term3,term4,u_slip,shear_max,shear_min,temp_max,t,EPSx,EPSy,um,ENDA) &
        !$OMP&            PRIVATE(J,I,temp_iter,temp_0,RL_Ta,EDA0,Z,EpsT,J0,J1,J2,JJ,I0,I1,I2,dPdx,dPdy,h_real,tau_x_f,tau_y_f,tau_x_c,tau_y_c)
        DO J=1,NN,SS
            DO I=1,NX,SS
                
                IF((I .NE. 1) .AND. (J .NE. 1) .AND. (.NOT. temp_1_1_initialized)) then   !Assuming the 1,1 iteration is the first iteration assigned to one of the threads. This is correct when using the STATIC schedule (which we do here since its the default)
                    DO WHILE (.NOT. temp_1_1_initialized)
                        !$OMP FLUSH (temp_1_1_initialized) 
                    ENDDO
                ENDIF
                
                temp_iter=0             ! Resetting counter
                temp_0=temp(1,1)        ! The global temperature
 
777             RL_Ta=temp(I,J)         ! Extract the temperature                  ! !!! Branch target 777
                
                ! Roelands equation acc R.Larsson
                EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135.0)**S0)             ! Eq (3)
                Z=Dz+Cz*log10(1.0+RL_Ta/135.0)
                
                EDAx(i,j)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                RO(I,J)=(1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    !$OMP CRITICAL (C_10)
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    !$OMP END CRITICAL (C_10)
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    !$OMP CRITICAL (C_11)
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    !$OMP END CRITICAL (C_11)
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) then
                    EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                    
                ELSE IF( temp_iter .LE. 10 .AND. SS .LE. 2) THEN     ! Might be good to introduce this after the coursest solution is obtained. 
                    
                    ! Define nodenumbers for past and next nodes
                    J0=J-SS
                    IF( J0 .LE. 0)  J0=J+SS
                    J2=J-2*SS
                    IF( J2 .LE. 0)  J2=J0
                    J1=J+1*SS
                    IF( J1 .GE. NYs) J1=NYs-SS
                    JJ=NYs+1-J
                    IF( JJ .LE. 1)  JJ=1+SS
                
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I2=I-2*SS
                    IF( I2 .LE. 0)  I2=I0
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                
                    dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))/(2*DX*SS)*Ph/b
                    dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))/(2*DX*SS)*Ph/b
                
                    h_real=H(I,J)*b**2/Rx
                    tau_x_f = abs(-h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)
                    tau_y_f = abs(-h_real/2*dPdy)
                    tau_x_c = abs( h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)
                    tau_y_c = abs( h_real/2*dPdy)
                    
                    
                    if( (tau_x_f .GT. shear_max .or. tau_y_f .gt. shear_max .or. tau_x_c .gt. shear_max .OR. tau_y_c .GT. shear_max) .AND. temp(i,j) .LT. temp_max) then                !D_IN: /shear_lim/ -> shear_max, temp_max
                        
                        temp(I,J)=temp(I,J)+5
                        RL_Ta=temp(I,J)
                        ! Roelands equation acc R.Larsson
                        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135.0)**S0)             ! Eq (3)

                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                        ! D-H Formulation acc R.Larsson
                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                        RO(I,J)=(1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                
                
                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                            !$OMP CRITICAL (C_12)
                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                            !$OMP END CRITICAL (C_12)
                            EDAx(I,J)=0.1
                            Call Stop_to_large_out(t)
                        ENDIF

                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                            !$OMP CRITICAL (C_13)
                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                            !$OMP END CRITICAL (C_13)
                            RO(I,J)=1.0
                            Call Stop_to_large_out(t)
                        ENDIF
        
                        IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                        
                        
                        
                        temp_iter=temp_iter+1
                        GO TO 777                      ! !!! Branch point 777
                    elseif( (tau_x_f .LT. shear_min .AND. tau_y_f .LT. shear_min .AND. tau_x_c .LT. shear_min .AND. tau_y_c .LT. shear_min) .AND. temp(i,j) .GT. temp_0) then      !D_IN: /shear_lim/ -> shear_min 
                        
                        temp(I,J)=temp(I,J)-2
                        RL_Ta=temp(I,J)
                        ! Roelands equation acc R.Larsson
                        EDA0 = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)             ! Eq (3)
                        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx

                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                        ! D-H Formulation acc R.Larsson
                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                        RO(I,J)=(1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                
                
                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                            !$OMP CRITICAL (C_14)
                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                            !$OMP END CRITICAL (C_14)
                            EDAx(I,J)=0.1
                            Call Stop_to_large_out(t)
                        ENDIF

                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                            !$OMP CRITICAL (C_15)
                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                            !$OMP END CRITICAL (C_15)
                            RO(I,J)=1.0
                            Call Stop_to_large_out(t)
                        ENDIF
        
                        IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))

                        
                        temp_iter=temp_iter+1
                        GO TO 777                      ! !!! Branch point 777
                    endif
                ENDIF
    
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))

                IF ((I .EQ. 1) .AND. (J .EQ. 1)) then   
                    temp_1_1_initialized = .true.
                    !$OMP FLUSH (temp_1_1_initialized)   
                ENDIF

            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
        !Mirroring the temperature
        DO J=1,NN,SS
            JJ=NYs-J+1
            DO I=1,NX,SS   
                temp(I,JJ)=temp(I,J)
            ENDDO
        ENDDO   
        
        
    ELSE IF( Lub_param .EQ. 5) THEN
        ! Yasutomi lubrication
        EDA22=Yedag*EXP(Yalfag*1)      !Done outside to avoide making every thread calculate this every iteration
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,PH,RA1,RA2,Tcvot,Temp,T40,P,Yedag,Yalfag,t,EDAx,EDA0,EDA22,RO,contact,EDA_cont,H,ENDA,EPSx,EPSy) &
        !$OMP&            PRIVATE(J,I,Tg,pg,YF,EDA1,EDA33)
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Yasutomi equations
                Tg=Tg0+YA1*LOG(1.0+YA2*PH*P(I,J))           !D_IN: /Yasutomi/ -> Tg0, YA1, YA2; 
                pg=1.0/YA2*(EXP((1.0/YA1)*((T40*(1.0-Tcvot)+Tcvot*Temp(I,J))-Tg0))-1.0)     
                YF=1.0-YB1*LOG(1.0+YB2*PH*P(I,J))           !D_IN: /Yasutomi/ -> Yb1, Yb2
                
                ! Since for some specific combinations of p and YF EDA33 gest biger than EDA22 eaven do p < pg this formulation should be better. 
                IF (P(I,J)*PH .GT. pg) THEN
                    EDA1=Yedag*EXP(Yalfag*(P(I,J)*PH-pg))   !D_IN: /Yasutomi/ -> Yedag
                ELSE if(Temp(i,j) .LT. TG) then
                    EDA1=Yedag*EXP(Yalfag*(P(I,J)*PH-pg))   !D_IN: /Yasutomi/ -> Yalfag
                    
                    !$OMP CRITICAL (C_16)
                    WRITE(4,*)'Selfcomponed case for temp < Tg. Temp and Tg are', Temp(i,j), Tg
                    WRITE(4,*)'at i,j, =', i,j
                    !$OMP END CRITICAL (C_16)
                    Call Stop_to_large_out(t)
                    
                ELSE
                    EDA33=Yedag*10**(-(YC1*((T40*(1-Tcvot)+Tcvot*Temp(I,J))-Tg)*YF) / (YC2+((T40*(1-Tcvot)+Tcvot*Temp(I,J))-Tg)*YF)) !D_IN: /Yasutomi/ -> YC1, YC2                    
                    
                    if (EDA22 .lt. EDA33) then
                    EDA1 = EDA22
                    else
                    EDA1 = EDA33
                    endif
                ENDIF
                
                EDAx(I,J)=EDA1                                              ! No Non-newtonian reduction. 
                
                 IF (EDAx(I,J) .LE. 0.001*EDA0 ) THEN ! With themp increase outside of pressure zone EDAx can decrease
                     IF( J .GT. 5) then
                        !$OMP CRITICAL (C_17)
                        WRITE(4,*)'BAD EDAx. EDAx lt 0.1*EDA0. EDAx=', EDAx(I,J), 'For I J = ', I ,J ,'EDA1 = ', EDA1
                        WRITE(4,*)'EDA33 = ', EDA33, ' EDA22 = ', EDA22 ,'Tcvot = ', Tcvot, ' Temp(i,j) = ', temp(i,j)
                        !$OMP END CRITICAL (C_17)
                     ENDIF
                     
                    EDAx(I,J)=0.001*EDA0
                    Call Stop_to_large_out(t)
                ENDIF
                
                IF( isnan(EDAX(i,j))) then
                    EDAX(i,j)=EDA1
                endif
                IF( (EDAX(i,j)) .gt. 1e15) then
                    !temp(1,1)=temp(i,j)
                    !P(1,1)=P(i,j)
                    EDAX(i,j)=1e15
                endif
            
                ! D-H Formulation acc P.Ehret D.Dowsin and C.M. Taylor
                RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                 
            ENDDO
        ENDDO
        !$OMP END  PARALLEL DO 
    ELSE IF( Lub_param .EQ. 6) THEN
        ! Newtonian acc N Deolalikers contact paper from 2008
        ! Does not exist in Lubrication_def
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NX,NN,SS,EDAy,EDAx,RO,xi,P,H,contact,PH,Pref,Z,EDA0,RA1,RA2,t,EPSx,EPSy,EpsT,EDA_cont,ENDA) &
        !$OMP&            PRIVATE(J,I)
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(log(EDA0+9.67)*((1+P(I,J)*PH/Pref)**Z-1)) ! Should include EDA0 due to fomulation in main
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
                ! D-H Formulation acc X.Tan and N Deolaliker
                ! RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    !$OMP CRITICAL (C_18)
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    !$OMP END CRITICAL (C_18)
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    !$OMP CRITICAL (C_19)
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    !$OMP END CRITICAL (C_19)
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
    ELSE IF( Lub_param .EQ. 7 .or. Lub_param .EQ. 8)THEN      
        ! Wang 2015 formulation  with temperature calculations or Liu 2005 formulation  with temperature calculations
         eda01=log(EDA0)+9.67
         !$OMP PARALLEL DO IF(use_multiple_cores) &
         !$OMP&            SHARED(NN,SS,NX,NYs,temp,RL_T0,S0,PH,Z,P,H,eda01,EDAx,RO,EpsT0,RA1,RA2,EDA0,t,EpsT,contact,EDA_cont,ENDA,EPSx,EPSy) &
         !$OMP&            PRIVATE(I,J,J0,J1,J2,JJ,I0,err,RL_Ta,eda02,eda03)
         DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J2=J-2*SS
            IF( J2 .LE. 0)  J2=J0
            J1=J+1*SS
            IF( J1 .GE. NYs) J1=NYs-SS
            JJ=NYs+1-J
            IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
                err=0
888             RL_Ta=temp(I,J)                      ! Branch target 888 - might be replaced with a Fortran equivalent to Continue

                
                eda02=(1.0+5.1E-9*P(I,J)*PH)**Z
                eda03=((RL_Ta-138.0)/(RL_T0-138.0))**(-S0)
                
                EDAx(I,J)=EDA0*exp(eda01 * (-1.0 + eda02 * eda03))
                
                !EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                RO(I,J)=1.0 + RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)) - EpsT0*(RL_Ta-RL_T0) 
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_20)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_20)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 888             ! !!! Branch to 888
                    else
                        !$OMP CRITICAL (C_21)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_21)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_22)
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_22)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 888            ! !!! Branch to 888
                    else
                        !$OMP CRITICAL (C_23)
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_23)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)
            enddo
         enddo
         !$OMP END PARALLEL DO
         
    ELSE IF( Lub_param .EQ. 9)THEN
        eda01=log(EDA0)+9.67
        ! Bruyerer 2012 formulation
        ! Also used for lub_param=10, Hartinger 2008
        !$OMP PARALLEL DO IF(use_multiple_cores) &
        !$OMP&            SHARED(NN,SS,NX,NYs,contact,EDA_cont,EPSx,EPSy,temp,P,H,RO,eda01,RL_T0,Z,PH,EDAx,EDA0,RA1,RA2,EpsT0,t,EpsT,ENDA,S0) &
        !$OMP&            PRIVATE(J,J0,J1,J2,JJ,I,I0,err,RL_Ta,eda02,eda03,eda04,delta_s)
         DO J=1,NN,SS
                    J0=J-SS
                    IF( J0 .LE. 0)  J0=J+SS
                    J2=J-2*SS
                    IF( J2 .LE. 0)  J2=J0
                    J1=J+1*SS
                    IF( J1 .GE. NYs) J1=NYs-SS
                    JJ=NYs+1-J
                    IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
855             err=0
                RL_Ta=temp(I,J)                 !Branch target 855
                
                eda02=(1.0+P(I,J)*PH/(1.98*1E8))**Z
                eda03=((RL_Ta-138.0)/(RL_T0-138.0))
                
                delta_s=eda01*eda02*S0/(RL_T0-138.0)
                
                eda04=delta_s*(RL_Ta-RL_T0)
                
                EDAx(I,J)=EDA0*exp(eda01 * ((-1.0 + eda02 )* eda03)-eda04)
                
                !EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                RO(I,J)=(5.9*1E8 + RA1*PH*P(I,J))/(5.9*1E8 + RA2*PH*P(I,J)) - EpsT0*(RL_Ta-RL_T0) 
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_24)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_24)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 855                   ! !!! branch to 855
                    else
                        !$OMP CRITICAL (C_25)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_25)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_26)
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_26)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 855                   ! !!! Branch to 855
                    else
                        !$OMP CRITICAL (C_27)
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_27)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)
            enddo
         enddo
        !$OMP END PARALLEL DO
         
    ELSE IF( Lub_param .EQ. 12 .or. Lub_param .EQ. 13)THEN      
        ! Wang 2015 formulation  with temperature calculations or Liu 2005 formulation  with temperature calculations
         eda01 = log(EDA0)+9.67
         !$OMP PARALLEL DO IF(use_multiple_cores) &
         !$OMP&            SHARED(NN,SS,NYs,NX,eda01,RA1,RA2,PH,RL_T0,EDA0,P,H,EDAx,temp,RO,EPSx,EPSy,contact,EDA_cont,ENDA,Z,S0,t,EpsT) &
         !$OMP&            PRIVATE(J,J0,J1,J2,JJ,I,I0,err,RL_Ta,eda02,eda03)
         DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J2=J-2*SS
            IF( J2 .LE. 0)  J2=J0
            J1=J+1*SS
            IF( J1 .GE. NYs) J1=NYs-SS
            JJ=NYs+1-J
            IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
889             err=0
                RL_Ta=temp(I,J)                             ! !!! Branch target 889

                eda02 = (1.0+5.1E-9*P(I,J)*PH)**Z
                eda03 = S0*(RL_Ta-RL_T0)       
                EDAx(I,J) = EDA0*exp(eda01 * (-1.0 + eda02) - eda03)

                ! D-H Formulation acc R.Larsson
                RO(I,J)=(1.0 + RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*( 1- EpsT0*(RL_Ta-RL_T0)) 
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_28)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_28)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 889                           ! !!! Branching point
                    else
                        !$OMP CRITICAL (C_29)
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_29)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_30)
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_30)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 889
                    else
                        !$OMP CRITICAL (C_31)
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_31)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)
            enddo
         enddo
         !$OMP END PARALLEL DO
    else
        WRITE(4,*)'Bad lubrication number in subroutine Newtonian. lub_param ='
        WRITE(4,*) lub_param
        stop 'Bad Lubrication'
    ENDIF
    
    RETURN
    END
