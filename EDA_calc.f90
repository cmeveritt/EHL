! Subroutine for updating the EDA = epsilon = density * film thickness^3 / (viscosity * lambda) term in the Reynolds equation
    ! Where lambda = 12 * average spped * radius^2 / contact halfwidth^3 * Hertzian pressure)
    ! based on the average temperature. Created to stabilize the thermal simulations
    ! Idealy EDA would only be calculated in one place to make sure that the same formulation is used 
	SUBROUTINE EDA_calc(SS, NYs, T40, Tcvot, NN, t, k, dPdx, dPdy, P_ave)
        implicit none 
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_Contact_mat.h'
        include     'inc_EDA_contact.h'
        include     'inc_EDA_average.h'
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
        integer     NN, NYs, SS, t, k
        real        T40, Tcvot
        ! Calculations
        integer     I, J, I0, I1, J0, J1, JJ, I2, J2, im, jm
        integer     temp_iter, err, printer
        real        EDA_1, EDA_2, EDA_3
        real        EDA22, EDA33, EDA1
        real        pg, Tg, YF 
        real        temp_0, u_slip, h_real, tau
        real        tau_x_f, tau_y_f, tau_x_c, tau_y_c, tau_lim_c
        real        P_ave(NX,NY), dPdx(NX, NY), dPdy(NX, NY), dpds					    !D_IN: /Grid/ -> Nx
        real        eda01,eda02,eda03,eda04,delta_s
        real        tau_gamma_r, dxss_Phb
        real        epst
        logical::   temp_1_1_initialized = .true.
        ! Output
        SAVE      /Current_EDAx_a/
        
        
    ! Calculate the viscocity, the density and the dimentionles parameter EPS 
    ! Only implemented for Lub_param = 4, R.Larssons 2000 formulation
    IF( Lub_param .EQ. 4 .and. lub_temp .GE. 0)THEN                                     !D_IN: /Ref/ -> lub_param; /temp_param/ -> lub_temp
        if ( temp_param == 2 ) then                                                     !D_IN: /temp_reduction/ -> temp_param
            ub= 2*um-ua                                                                 !D_IN: /Outp/ -> Um, Ua;       !D_OUT: Ub -> /Outp/
            u_slip=ua-ub
        endif
        
        printer=0
        
        ! Calcualte the distance between nodes i m. DX is dimensionless
        dxss_Phb=1.0/(2*DX*SS)*Ph/b                                                     !D_IN: /Grid/ -> Dx; /Visc/ -> Ph; /Outp/ -> B
        
        ! Calculate for the whole half surface
        !$OMP PARALLEL DO IF(use_multiple_cores)&
        !$OMP&            SHARED(NN,SS,NX,temp_param,temp,temp_fac,RL_G0,S0,Dz,Cz,PH,b,Rx,P,EDAx_a,p_red,temp_red,H,dpdx,dpdy,U_slip,tau_gamma,tau_lim,lim_factor,t,contact,EDA_cont) &
        !$OMP&            PRIVATE(J,J0,J1,J2,JJ,I,I0,I1,I2,err,RL_Ta,EDA0,Z,EDA1,h_real,tau,tau_gamma_r,dpds,tau_lim_c) &
        !$OMP&            REDUCTION(+:printer)
        DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J2=J-2*SS
            IF( J2 .LE. 0)  J2=J0
            J1=J+1*SS
            IF( J1 .GT. NN) J1=J-SS
            JJ=NYs+1-J
            IF( JJ .LE. 1)  JJ=1+SS
            
            DO I=1,NX,SS
899                 err=0                                                             
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I2=I-2*SS
                    IF( I2 .LE. 0)  I2=I0
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                    
                IF((I .NE. 1) .AND. (J .NE. 1) .AND. (.NOT. temp_1_1_initialized)) then   !Assuming the 1,1 iteration is the first iteration assigned to one of the threads. This is correct when using the STATIC schedule (which we do here since its the default)
                        DO WHILE (.NOT. temp_1_1_initialized)
                        !$OMP FLUSH (temp_1_1_initialized) 
                        ENDDO
                ENDIF
                
                ! Scaling the temperatures with temp_fac
                If ( Temp_param == 3) then
                    RL_Ta=temp(1,1) +( temp(I,J) - temp(1,1))*temp_fac                 !D_OUT: RL_Ta -> /RLarsson/(Not saved however)                   !D_IN: /CurrentT/ -> temp; /Temp_reduction/ -> temp_fac
                    IF( RL_Ta .LT. temp(1,1)) then
                        RL_ta=temp(1,1)
                    endif
                else
                    ! Calculate an somewhat averaged temperature based on the temperature in the neibouring nodes. 
                    RL_Ta=( temp(I0,J) + temp(I,J0) + 4*temp(I,J) + temp(I1,J) + temp(I1,J)) / 8.0
                endif
                
                ! Roelands equation for the viscosity EDA1 acc R.Larsson 2000
                EDA0 = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)             
                Z=Dz+Cz*log10(1.0+RL_Ta/135)
                EDA1=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(i,j)*PH)**Z))														                !D_IN: /Visc/ -> EDA0, Ph, Z; /CurrentP/ -> P(I,J)
               
                ! If slip a shear limit or simliar could be used to decrease the effective viscosity
                IF( U_slip == 0) then
                    EDAx_a(i,j) = EDA1   
                                                                                                                                                        !D_OUT: EDAx_a(I,J) -> /Current_EDAx_a/                                            
                ! Reducing eta    
                elseif ( Temp_param == 1) then
                    
                    EDAx_a(i,j) = EDA1 *  1/ (1 + temp_fac * max( 0.0, (P(i,j)-p_red)/p_red) * ( max(0.0 , (RL_Ta-temp_red)/temp_red ) ) )              !D_IN: /Temp_reduction/ -> temp_red, p_red; 
                
                ! applying a shear limmit  based on the pressure according to data from Hoglund 1989 and 1999 with the temperature reduction of the shear limit                                                                                                                                                              
                elseif( temp_param == 2 .or. temp_param==6) then
                    
                    h_real=H(I,J)*b**2/Rx                                                                                                               !D_IN: /outp/ -> b, Rx; /CurrentH/ -> H(I,J)
                    tau = abs(H_real/2*dpdx(i,j)) + abs(EDA1*U_slip/H_real) 
                    tau= sqrt( tau**2 + (H_real/2*dpdy(i,j))**2)
                    
                    if (temp_red .GT. 0) then                           ! Take calculate temperature reductio of the shear limit
                        tau_gamma_r = tau_gamma - temp_fac * RL_TA      ! tau_gamma and temp_fac anready includes PH                                    !D_IN: /Temp_reduction/ -> tau_gamma
                        if( tau_gamma_r .LT. temp_red) then             ! Just to make sure that not negative
                            tau_gamma_r = temp_red
                        endif
                    else
                        tau_gamma_r = tau_gamma
                    endif
                    
                    if (temp_param==2)then
                        tau_lim_c = tau_lim + tau_gamma_r * P(i,j)                                                                                      !D_IN: /Temp_reduction/ -> tau_lim
                    elseif ( temp_param == 6) then
                        tau_lim_c = tau_lim + tau_gamma_r * P(i,j) + tau*p_red
                    endif
                    
                    ! If the shear stresses are too high. Then they have to be limited and therefor has the viscosity also to be reduced. Since the viscosity is what is used in the Reynolds equation. 
                    if (tau .gt. tau_lim_c) then
                        dpds= sqrt( dPdx(i,j)**2 + dPdy(i,j)**2)
                        if( tau_lim_c - abs(H_real/2*dpds) .LT. lim_factor*tau_lim ) then                                                               !D_IN: /Temp_reduction/ -> lim_factor
                            
                            ! To increase numerical stability the maximum pressure of the neighbouring nodes is used
                            if (temp_param==2)then
                                tau_lim_c = tau_lim + tau_gamma_r *max( P(i0,j),P(i,j),P(i1,j), P(i,j0), P(i,j1))
                            elseif ( temp_param == 6) then
                                tau_lim_c = tau_lim + tau_gamma_r *max( P(i0,j),P(i,j),P(i1,j), P(i,j0), P(i,j1)) + tau*p_red
                            endif
                            
                            if( tau_lim_c - abs(H_real/2*dpds) .LT. lim_factor*tau_lim ) then
                                EDA1 = abs( (lim_factor*tau_lim)*H_real/U_slip)    
                                printer=printer+1
                            else
                                EDA1 = abs( (tau_lim_c - abs(H_real/2*dpds))*H_real/U_slip)
                            endif
                            
                        else
                            EDA1 = abs( (tau_lim_c - abs(H_real/2*dpds))*H_real/U_slip)
                        endif
                    endif
                    
                    EDAx_a(i,j) = EDA1
                elseif( temp_param== 4 .or. temp_param== 5) then ! affects the tempt_calc_ij subroutine
                    EDAx_a(i,j) = EDA1
                else
                    EDAx_a(i,j) = EDA1
                endif

                IF (EDAx_a(I,J) .LE. 0.0 .OR. isnan(EDAx_a(I,J)) .or. (EDAx_a(I,J) .GT. 1e20)  ) THEN
                    err=err+1
                    if (err==1) then
                        !$OMP CRITICAL (C_32)
                        WRITE(4,*)'BAD EDAx_a =', EDAx_a(I,J), '. For I J = ', I ,J , 'P = ', P(i,j), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_32)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 899
                    else
                        !$OMP CRITICAL (C_33)
                        WRITE(4,*)'BAD EDAx_a =', EDAx_a(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(i,j), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        !$OMP END CRITICAL (C_33)
                        err=1
                        EDAx_a(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)                     !CALL: No comment 
                ENDIF
                
                ! If metal contact, defined by contact = 2, then ensure that the viscosity is at least equal to EDA_cont
                IF(contact(I,J) .EQ. 2 .and. EDAx_a(I,J) .LT. EDA_cont) EDAx_a(I,J)=EDA_cont                                                            !D_IN: /Contact_mat/ -> contact(i,j); /EDA_contact/ -> EDA_cont
                
                IF ((I .EQ. 1) .AND. (J .EQ. 1)) then   
                    temp_1_1_initialized = .true.
                    !$OMP FLUSH (temp_1_1_initialized)   
                ENDIF
                
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
        if (printer .gt. 20 .and. k == 1) then
            WRITE(4,*)'Bad equation for shear limmit in EDA_clac ', printer ,' number of times'
            Call Stop_to_large_out(t)                             !CALL: No comment
            if ( printer .gt. 200 .and. SS == 1) then
                WRITE(4,*) 'Breaking due to many bad shear limit equations in EDA_clac'
                WRITE(4,*) 'The time is t=',t
                OPEN(20,FILE='STOP_Too_many_bad_Shear.DAT',STATUS='UNKNOWN')
                stop
            endif          
        endif
   
    else
        WRITE(4,*)'Bad lubrication number in subroutine EDA_clac. lub_param ='
        WRITE(4,*) lub_param
        !stop 'Bad Lubrication'
    ENDIF
    
    RETURN
    END
