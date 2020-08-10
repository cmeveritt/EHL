! Subrutine for updating the film  the viscocity and the density based on the Liu shear thinning equations. 
    ! This subroutine is not compleated and does not work correctly
	SUBROUTINE Liu_shear_thin(SS, NYs, T40, Tcvot, NN)
        implicit none
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
        include     'inc_Contact_mat.h'
        include     'inc_Grid.h'
        include     'inc_Method.h'
        include     'inc_Outp.h'
        include     'inc_Ref.h'
        include     'inc_RLarsson.h'
        include     'inc_Shear_nbr.h'
        include     'inc_Visc.h'
        include     'inc_Y_liu.h'                
        ! Input
        integer     SS
        ! Calculations
        real        dpdx, dpdy
        real*8      gauss(10,3), xii, wi
        integer     I,J, k,JJ, iter
        integer     J1,J0, J2, I0,I1, I2, IA, II
        real        EDA_cont
        real        Px, Py, ta, tb, fp
        real        I0p, I1p, I2p, J0p, J1p, J2p, I0t, I1t, I2t, J0t, J1t, J2t 
        real        dtau(2), tt
        real        dus,dvs, tt_old
        real        ty, tx, txp, te
        real        aa(2), bb(2), cc(2)
        real        Us, Ul 
        real        uss, vs, h_real, limit    
        ! Other
        Integer     NN, ink, NYs
        real        Tau, EDA1, tau1, tau2
        real        T40, Tcvot
        ! Output
        SAVE        /Current/                   
        SAVE        /shear_nbr/
        
        WRITE(4,*) 'Warning, The Liu shearthinning subroutine does not work properly' 

        EDA_cont=EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*3.0*PH)**Z)) !D_IN: /Visc/ -> EDA0 
        
        ! Slip velocity
        uss=ua-ub                     !D_IN: /Outp/ -> Ua, Ub
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
        
        
        ! NonNewtonian Liu imput for Lub_param = 60
        ! Add so that if contact has occured the viscosity gets really high. Othervise the code will have a hard time converging for contact. 
        DO J=1,NN,SS
            DO I=1,NX,SS                                        !D_IN: /Grid/ -> Nx
                
                 IF(contact(I,J) .EQ. 2 ) THEN                  !D_IN: /Contact_mat/ -> contact
                    !If contact no shear thinning
                 ELSE
                     
                ! Define nodenumbers for past and next nodes
                J0=J-SS
                IF( J0 .LE. 0)  J0=J+SS
                J2=J-2*SS
                IF( J2 .LE. 0)  J2=J0
                J1=J+1*SS
                IF( J1 .GE. NY) J1=NY-SS
                JJ=NY+1-J
                IF( JJ .LE. 1)  JJ=1+SS
                I0=I-1*SS
                IF( I0 .LE. 0)  I0=1
                I2=I-2*SS
                IF( I2 .LE. 0)  I2=I0
                I1=I+1*SS
                IF( I1 .GE. NX) I1=NX
                
                ! Real film thickness in m
                h_real=H(I,J)*b**2/Rx    !D_IN: /CurrentH/ -> H(I,J); /Outp/ -> b,Rx
                
                limit=3e8/L_G            !D_IN: /Y_liu/ -> L_G                 
                EDA1=EDAx(i,j)           !D_IN: /Current/ -> EDAx(i,j)
                
                ! Pressure derivatives in RD and TD
                dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))/(2*DX*SS)             !D_IN: /Grid/ -> Dx; /Method/ -> term1, term2, term3, term4; /CurrentP/ -> P
                dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))/(2*DX*SS)
                
                ! Dimensionless components
                dUs=uss*EDA1/(h_real*L_G);
                dVs=vs*EDA1/(h_real*L_G);
                Px=h_real/L_G*dPdx*Ph;
                Py=h_real/L_G*dPdy*Ph;
                
                ! Start guess
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
                            CALL f_pure_roll(0.5*xii, PX, Py, L_n, fp)                 !D_IN: /Y_Liu/ -> L_n
                            wi  = gauss(k,2)
                            txp = txp + wi* 0.5 *fp
                        enddo
                        tx=txp*12
                        ty=tx
                    endif
                    
                ! If sliding
                else
                    iter=1;
                    DO while(tt .gt. 1e-6 .and. iter .lt. L_iter)                      !D_IN: /Y_liu/ -> L_iter
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

                        aa(1)= I0p+ ta*J0p*ta+ 2*ta*Px*J1p+        Px**2*J2p;
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
                    Shear_iter=shear_iter+1             !Counter for how many times the shear thinning does not converge     !D_OUT: Shear_iter -> /shear_nbr/                     
                    shear_incre=max(shear_incre,tt)     !The biggest nonconverged shear increment this itteration            !D_OUT: Shear_incre -> /shear_nbr/ 
                endif
                
                ! Add function for EDAx(x) and EDAY(tb) ?
               
                EPSx(I,J)=tx*EPSx(I,J)! RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))                                                   !D_OUT: EPSx(I,J) -> /Current/
                EPSy(I,J)=ty*EPSx(I,J)! RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))                                                   !D_OUT: EPSy(I,J) -> /Current/
                
                IF (tx .LT. 1.0 .OR. isnan(tx) ) THEN
                    IF (tx .LT. 0.9 .OR. isnan(tx) ) then
                        WRITE(4,*)'BAD tx, tx = ', tx, 'and EPSx =', EPSx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                        Call Stop_to_large_out(0) ! THe zero is used since the time is not defined here
                    endif
                    tx=1
                    
                ENDIF
                IF (ty .LT. 1.0 .OR. isnan(ty) ) THEN
                    IF (ty .LT. 0.9 .OR. isnan(ty) ) WRITE(4,*)'BAD ty, ty = ', ty, 'and EPSy = ', EPSy(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    ty=1
                    Call Stop_to_large_out(0) ! THe zero is used since the time is not defined here
                ENDIF
                                
            ENDIF
                 
                
        ENDDO
    ENDDO
    
    RETURN
    END

    ! Subfunctions
    SUBROUTINE ICALC(te,ta,tb,xi,n,Px,Py, I0, I1, I2)
    implicit none 
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f
    real*8  xi
    
            CALL I0p(te,ta,tb,xi,n,Px,Py, I0) 
            CALL I1p(te,ta,tb,xi,n,Px,Py, I1) 
            CALL I2p(te,ta,tb,xi,n,Px,Py, I2) 
    end
    
    SUBROUTINE JCALC(te,ta,tb,xi,n,Px,Py, J0, J1, J2)
    implicit none 
    real    te,ta,tb,n,Px,Py, J0, J1, J2,f
    real*8  xi
    
            CALL J0p(te,ta,tb,xi,n,Px,Py, J0) 
            CALL J1p(te,ta,tb,xi,n,Px,Py, J1) 
            CALL J2p(te,ta,tb,xi,n,Px,Py, J2) 
    end
    
    ! To calculate the effective shear stress
    SUBROUTINE Y_te(ta,tb,xi,Px,Py, te) 
    implicit none 
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f
    real*8  xi
    
        te= sqrt((ta+xi*Px)**2+(tb+xi*Py)**2)
    end
    
    ! To calculate f(te)
    SUBROUTINE Y_f(te,ta,tb,xi,n,Px,Py, f) 
    implicit none
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f
    real*8  xi
    
        CALL Y_te(ta,tb,xi,Px,Py, te) 
        
        f=(1+te**2)**((1-n)/(2*n))
    end
    
    !To calculate f(te)/te
    SUBROUTINE Y_f2(te,ta,tb,xi,n,Px,Py, f2)
    implicit none
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f2
    real*8  xi
    
        CALL Y_te(ta,tb,xi,Px,Py, te) 
        
        f2=  (1-n)/(n)*(1+te**2)**((1-3*n)/(2*n))
    end
    
    ! To Calculate the I functions
    SUBROUTINE I0p(te,ta,tb,xi,n,Px,Py, I0)
    implicit none
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f
    real*8  xi
    
               CALL Y_f(te,ta,tb,xi,n,Px,Py, f) 
               I0 = 1*f
    end
    
    SUBROUTINE I1p(te,ta,tb,xi,n,Px,Py, I1) 
    implicit none
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f
    real*8  xi
    
                CALL Y_f(te,ta,tb,xi,n,Px,Py, f) 
                I1 = xi* f
    end
    
    SUBROUTINE I2p(te,ta,tb,xi,n,Px,Py, I2)
    implicit none
    real    te,ta,tb,n,Px,Py, I0, I1, I2,f
    real*8  xi
    
                CALL Y_f(te,ta,tb,xi,n,Px,Py, f) 
                I2 = xi**2*f
    end

    ! To calculate the J functions
    SUBROUTINE J0p(te,ta,tb,xi,n,Px,Py, J0) 
    implicit none
    real    te,ta,tb,n,Px,Py, J0, J1, J2,f2
    real*8  xi
    
                CALL Y_f2(te,ta,tb,xi,n,Px,Py, f2) 
                J0 = 1*f2
    end
    
    SUBROUTINE J1p(te,ta,tb,xi,n,Px,Py, J1) 
    implicit none
    real    te,ta,tb,n,Px,Py, J0, J1, J2,f2
    real*8  xi
    
                CALL Y_f2(te,ta,tb,xi,n,Px,Py, f2) 
                J1=xi* f2
    end
    
    SUBROUTINE J2p(te,ta,tb,xi,n,Px,Py, J2) 
    implicit none
    real    te,ta,tb,n,Px,Py, J0, J1, J2,f2
    real*8  xi
    
                CALL Y_f2(te,ta,tb,xi,n,Px,Py, f2) 
                J2=xi**2* f2
    end
                
    SUBROUTINE f_pure_roll(xi, PX, Py, n, fp)
    implicit none
    real    te,ta,tb,n,Px,Py, I0, I1, I2, fp
    real*8  xi
    
            fp= xi**2*(1+xi**2*(Px**2+Py**2))**((1-n)/(2*n));
    end
    

