! Subroutine for updating the material parameters of the lubricant based on the Liu shear thinning model. 
    !This subroutine does not work properly
	SUBROUTINE Temp_Liu_shear_thin(SS, NYs, T40, Tcvot, NN, t)
        implicit none
        COMMON /Visc/       ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters                    
        COMMON /Grid/       NX,NY,X0,XE,DX                                      ! Grid parameters
        COMMON /Ref/        PAIAK, meth, tmeth, Geom, Lub_param, asp_shape, contact_alg, shear_thin      ! Coise of equations based on referense
        COMMON /RLarsson/   EpsT0,RL_Ta,S0, RL_G0, Dz, Cz,RL_T0, RL_c           ! Parameters for ref 10, R. Larssons formulation
        COMMON /outp/       W0,EE,RX,Um, Ua, B, U2_over_Us, BY, Ry 
        COMMON  /Method/    term1,term2,term3,term4                             ! Controling the numerical method
        COMMON  /Current/   RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                   ! Current timestep
        COMMON  /CurrentP/  P
        COMMON  /CurrentH/  H
        COMMON  /Rho/           RA1,RA2  
        COMMON  /COMT/          Temp(1:601,1:601)
        COMMON  /Y_Liu/     L_n, L_G, L_h_limit, L_iter, L_stab                 ! Lubrication and convergence parameters acc Y.Liu 2007 
        COMMON  /Contact_mat/   contact
        COMMON  /shear_nbr/ shear_incre, shear_iter, temp_iter_t
        ! Input
        real        P(1:601,1:601), H(1:601,1:601)
        real        um
        integer     NX, SS, t
        integer     contact(1:601,1:601)
        integer     term1,term2,term3,term4 
        real        dpdx, dpdy, RA1,RA2 
        real        L_n, L_G, L_h_limit, L_iter, L_stab
        real        EpsT0,RL_Ta,S0, RL_G0, Dz, Cz,RL_T0, RL_c, EpsT
        ! Calculations
        real*8      gauss(10,3), xii, wi
        integer     I,J, k,JJ, iter, temp_iter
        integer     J1,J0, J2, I0,I1, I2, IA, II
        real        EDA_cont, temp_0
        real        Px, Py, ta, tb, fp
        real        I0p, I1p, I2p, J0p, J1p, J2p, I0t, I1t, I2t, J0t, J1t, J2t 
        real        dtau(2), tt
        real        dus,dvs, tt_old
        real        ty, tx, txp, te
        real        aa(2), bb(2), cc(2)
        real        Ub, Us, Ul 
        real        uss, vs, h_real, limit    
        real        tau_x_f, tau_y_f, tau_x_c, tau_y_c
        ! Other
        real        W(1:601,1:601),Wside(1:601,1:601)                      ! Current timestep
        Integer     NY, NN, ink,  tmeth, NYs, Geom, Lub_param, meth, contact_alg
        Integer     asp_shape , shear_thin, temp_iter_t
        Real        X0, XE, PAIAK
        Real        RX, W0, Ua, BY, Ry 
        Real        Z,  EDA0, Pref, alpha
        Real        EE, DX, U2_over_Us
        real        ENDA, B, PH, HM0r
        real*8      temp
        real        A1, A2, A3, Tau, EDA1, tau1, tau2
        real        T40, Tcvot
        ! Output
        real        RO(1:601,1:601),EPSx(1:601,1:601),EPSy(1:601,1:601),EDAx(1:601,1:601),EDAy(1:601,1:601),xi(1:601,1:601)
        integer     shear_iter
        real        shear_incre
        SAVE        /Current/                    ! Current timestep
        SAVE        /shear_nbr/
        SAVE        /COMT/      
        
        WRITE(4,*) 'Warning! Temp_Liu_shear_thin does not work properly. This part of the code is under construction'
        Ub          =2*um-Ua
        EDA_cont=EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*3.0*PH)**Z))  
        
        uss=ua-ub
        vs=0

        
        !% Gausian positions
        !% https://pomax.github.io/bezierinfo/legendre-gauss.html
        !%       number  weight          xi coordinate
        gauss(1,1:3) = (/   1.00,	0.2955242247147529,	-0.1488743389816312 /)
        gauss(2,1:3) = (/   2.00,	0.2955242247147529,	0.1488743389816312  /)
        gauss(3,1:3) = (/   3.00,	0.2692667193099963,	-0.4333953941292472 /)
        gauss(4,1:3) = (/   4.00,	0.2692667193099963,	0.4333953941292472  /)
        gauss(5,1:3) = (/   5.00,	0.2190863625159820,	-0.6794095682990244 /)
        gauss(6,1:3) = (/   6.00,	0.2190863625159820,	0.6794095682990244  /)
        gauss(7,1:3) = (/   7.00,	0.1494513491505806,	-0.8650633666889845 /)
        gauss(8,1:3) = (/   8.00,	0.1494513491505806,	0.8650633666889845  /)
        gauss(9,1:3) = (/   9.00,	0.0666713443086881,	-0.9739065285171717 /)
        gauss(10,1:3)= (/   10.00,  0.0666713443086881,	0.9739065285171717  /)
        
        
       !NonNewtonian Liu imput for Lub_param = 60
        ! Add so that if contact has occured the viscosity gets really high. Othervise the code will have a hard time converging for contact. 
        DO J=1,NN,SS
            DO I=1,NX,SS
            temp_iter=0
            temp_0=temp(1,1)                        ! The global temperature
                
777         IF(contact(I,J) .EQ. 2 ) THEN
                !If contact no shear thinning
                !IF( EDAy(I,J) .LT. EDA_cont) ty=1    
                !IF( EDAx(I,J) .LT. EDA_cont) tx=1    !If contact no shear thinning
            ELSE
            
            ! Ensure that the lubrication can not cool off downstream. 
            IF( L_stab .EQ. 2 .AND. t .LE. -2 .and. I .GT. 1 ) Then
                temp(I,j)=max(temp(i,j),temp(I-SS,J))       ! Do not update if higher temp 
                temp_0=temp(I-SS,J)                         ! The minimum allowed temperature
            ENDIF

                 
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
                
                                
                h_real=H(I,J)*b**2/Rx
                limit=3e8/L_G


                
                EDA1=EDAx(i,j)
                
                dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))/(2*DX*SS)
                dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))/(2*DX*SS)
                
               ! %Dimensionless components
                dUs=uss*EDA1/(h_real*L_G);
                dVs=vs*EDA1/(h_real*L_G);
                Px=h_real/L_G*dPdx*Ph;
                Py=h_real/L_G*dPdy*Ph;
                
                !%start guess
                ta=dUs
                tb=dVs
                
                tt=1
                
                ! IF pure rolling
                if (dUs .eq. 0.0  .and. dVs .eq. 0.0)then
                    ! IF no preassure gradient
                    if( abs(px) .lt. 1e-3 .and. abs(py) .lt. 1e-3)then
                        tx=1
                        ty=1
                    ! If pressure gradient
                    else 
                        txp=0;

                        DO k=1,10
                            xii  = gauss(k,3)
                            CALL f_pure_roll(0.5*xii, PX, Py, L_n, fp)
                            wi  = gauss(k,2)
                            txp = txp + wi* 0.5 *fp
                        enddo
                        tx=txp*12
                        ty=tx
                    endif
                    
                ! If sliding
                else
                    iter=1;
                    DO while(tt .gt. 1e-6 .and. iter .lt. L_iter)
                        I0p=0
                        I1p=0
                        I2p=0
        
                        J0p=0
                        J1p=0
                        J2p=0
        
                        do k=1,10
                            xii = gauss(k,3);
                            wi = gauss(k,2);
                            
                            CALL ICALC(te,ta,tb,0.5*xii,L_n,Px,Py, I0t, I1t, I2t) 
                            CALL JCALC(te,ta,tb,0.5*xii,L_n,Px,Py, J0t, J1t, J2t) 
                            
                            I0p = I0p+  wi* 0.5* I0t
                            I1p = I1p+  wi* 0.5* I1t
                            I2p = I2p+  wi* 0.5* I2t
            
                            J0p = J0p+  wi* 0.5* J0t
                            J1p = J1p+  wi* 0.5* J1t
                            J2p = J2p+  wi* 0.5* J2t
                        enddo
                        if(( abs(Px) .LT. 1e-4 .and. abs(Py) .LT. 1e-4 .and. abs(I1p) .LT. 1e-4) .OR. abs( I1p) .LT. 1e-5 ) then ! te does almost not depend on xii and therefore
                            I1p=0.0
                        endif
                        
                        
                        aa(1)= I0p+ ta*J0p*ta+ 2*ta*Px*J1p+         Px**2*J2p;
                        aa(2)=      ta*J0p*tb+  (ta*Py+tb*Px)*J1p+  Px*Py*J2p;
                        bb(1)= aa(2);
                        bb(2)= I0p+ tb*J0p*tb+  2*tb*Py*J1p+        Py**2*J2p;
        
                        cc(1)= dUs- (ta*I0p+ Px*I1p);
                        cc(2)= dVs- (tb*I0p+ Py*I1p);
        
                        dtau(1) = ( cc(1)/aa(1) - bb(1)*cc(2) / (bb(2)*aa(1)) ) / ( 1-bb(1)*aa(2) / (bb(2)*aa(1)) )
                        dtau(2) = ( cc(2)/bb(2) - aa(2)*cc(1) / (bb(2)*aa(1)) ) / ( 1-bb(1)*aa(2) / (bb(2)*aa(1)) )
                        
                        if ( abs(dtau(1)) .gt. abs(dtau(2))) then
                            dtau(2) = (cc(2)-aa(2)*dtau(1))/bb(2)
                        else
                            dtau(1) = (cc(1)-bb(1)*dtau(2))/aa(1)
                        endif
                        
        
                        tt=sqrt(dtau(1)**2+dtau(2)**2)
                        
                        ta=ta+dtau(1)
                        tb=tb+dtau(2)
        
                        iter=iter+1;
                        tt_old=tt;
        
                    enddo
    
                    tau_x_f=abs(ta*L_G-0.5*h_real*px*L_G)
                    tau_y_f=abs(tb*L_G-0.5*h_real*py*L_G)
                    tau_x_c=abs(ta*L_G+0.5*h_real*px*L_G)
                    tau_y_c=abs(tb*L_G+0.5*h_real*py*L_G)
                    
                    IF( temp_iter .LE. 10 .AND. SS .LE. 2) THEN     ! Might be good to introduce this after the coursest solution is obtained. 
                    if( ( tau_x_f .GT. 1e9 .or. tau_y_f .gt. 1e9 .or. tau_x_c .gt. 1e9 .OR. tau_y_c .GT. 1e9) .AND. temp(i,j) .LT. 200) then
                        
                        temp(I,J)=temp(I,J)+5
                        RL_Ta=temp(I,J)
                        ! Roelands equation acc R.Larsson
                        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135)**S0)             ! Eq (3)
                        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx

                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                        ! D-H Formulation acc R.Larsson
                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                        RO(I,J)=(1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)))*(1-EpsT*(RL_Ta-RL_T0)) 
                
                
                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                            EDAx(I,J)=0.1
                            Call Stop_to_large_out(t)
                        ENDIF

                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                            RO(I,J)=1.0
                            Call Stop_to_large_out(t)
                        ENDIF
        
                        IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                        
                        
                        iter=0
                        GO TO 777
                        temp_iter=temp_iter+1
                    elseif( (tau_x_f .LT. 0.1e9 .AND. tau_y_f .gt. 0.1e9 .AND. tau_x_c .gt. 0.1e9 .AND. tau_y_c .GT. 0.1e9) .AND. temp(i,j) .GT. temp_0) then
                        
                        temp(I,J)=temp(I,J)-2
                        RL_Ta=temp(I,J)
                        ! Roelands equation acc R.Larsson
                        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135)**S0)             ! Eq (3)
                        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx

                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                        ! D-H Formulation acc R.Larsson
                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                        RO(I,J)=(1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)))*(1-EpsT*(RL_Ta-RL_T0)) 
                
                
                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                            EDAx(I,J)=0.1
                            Call Stop_to_large_out(t)
                        ENDIF

                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                            RO(I,J)=1.0
                            Call Stop_to_large_out(t)
                        ENDIF
        
                        IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                        
                        
                        iter=0
                        GO TO 777
                        temp_iter=temp_iter+1
                        
                    endif
                    ENDIF
                    
                    IF(temp_iter .GE. 10) temp_iter_t=temp_iter_t+1     ! The counter variable to be printed out
                        
                    !if (iter .lt. 100) then
                        if( abs(px) .LT. 1e-6) px=sign(1e-6,px)
                        tx=abs(12*(ta*I1p/Px+I2p));
                        if( abs(py) .LT. 1e-6) py=sign(1e-6,py)
                        ty=abs(12*(tb*I1p/Py+I2p));
                    !else 
                    !    tx=1;
                    !    ty=1;
                    !endif
    
                endif
                
                if( iter .GE. L_iter) then
                    Shear_iter=shear_iter+1             !Counter for how many times the shear thinning does not converge
                    shear_incre=max(shear_incre,tt)     !The biggest nonconverged shear increment this itteration
                endif
                
                
                IF (tx .LT. 1.0 .OR. isnan(tx) ) THEN
                    IF (tx .LT. 0.9 .OR. isnan(tx) ) then
                        WRITE(4,*)'BAD tx, tx = ', tx, 'and EPSx =', EPSx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                        Call Stop_to_large_out(t)
                    endif
                    tx=1
                    
                ENDIF
                IF (ty .LT. 1.0 .OR. isnan(ty) ) THEN
                    IF (ty .LT. 0.9 .OR. isnan(ty) ) WRITE(4,*)'BAD ty, ty = ', ty, 'and EPSy = ', EPSy(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    ty=1
                    Call Stop_to_large_out(t)
                ENDIF
                
                EPSx(I,J)=tx*EPSx(I,J)! RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=ty*EPSx(I,J)! RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                                
            ENDIF
                 
                
        ENDDO
        ENDDO
        
        !Mirroring the temperature
        DO J=1,NN,SS
            JJ=NYs-J+1
            DO I=1,NX,SS   
                temp(I,JJ)=temp(I,J)
            ENDDO
        ENDDO
        
    
    RETURN
    END