! Subrutine for updating the film thicknes, the viscocity and the density based on the given pressure
	SUBROUTINE HREE(H00,t)
        implicit none 
        COMMON      /HREEp/Emat, pg, Tg,YF                                       ! Parameters for HREE subrutine
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Setup/X,Y, P0                                       ! Parameters whihc do not change over time
        COMMON      /G0DT/G0,DT                                              ! G0 and DT
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real        Emat(1:300,1:300), pg(1:300,1:300), Tg(1:300,1:300),YF(1:300,1:300) 
        real        X(1:300),Y(1:300), P0(1:300,1:300)
        real        AK
        Integer     NX, NY, NN, ink, ref, t, tmeth
        Integer     I,J, JJ
        integer     Ntime, MK_stat, MK_time, KK
        integer     PK, deformation, loadupd
        real        load, loadOLD, Sload, SloadOLD, loadDiff, DW
        Real        X0, XE, PAI, PAIAK, Plim
        Real        RX, W0, Ua, Ub, RA1, RA2
        Real        Z,  EDA0, Pref, alpha
        Real        Elast1, Elast2, EE
        Real        asph_real, aspw_real, ER_stat, ER_time
        Real        DX, DT, width
        Real        asph, aspw,US, G0
        real        ENDA
        real        B, PH, H00, U, SRR
        real        tauc_real, tauc, tauS, taua, taua2,xilim
        real        Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
        real        kH, xH, gH, lH, L, M
        real        AB, WB, RAD, temp, position, Hmin, HM0r, H0
        real        DWold, bump, asp, A1, A2, A3, Tau, EDA1
        real        COSH, SINH
        COMMON /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
        COMMON /NonNew/ tauc, taua, taua2, xilim                               ! Non Newtonian viscosity parameters
        COMMON /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag    ! Lubrication parameters acc Yasutomi
        COMMON /Rho/RA1,RA2                                             ! Density parameters
        COMMON /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
        COMMON /Itt/MK_stat, MK_time, ER_stat, ER_time, Ntime, KK   ! Itteration paramters
        COMMON /asp/asph,aspw                                           ! Apserity parameters
        COMMON /Holmes/kH, xH, gH, lH                                   ! Viscosity and denisty param acc Holems et al
        COMMON /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referense
        COMMON /COMT/Temp(1:300,1:300)
        COMMON /COMAK2D/AK(0:300,0:300)
        COMMON /COMH/RAD(1:300,1:300)            
        SAVE      /Current/                    ! Current timestep
        DATA PAI/3.14159265/
        
        NN=(NY+1)/2
        PK=0
15      PK=PK+1
        ink=0
        
        Plim=15                     ! Limit on maximum allowed pressure. To help the solution converge
        
        ! Check that the pressure is not diverging --------------------------------------------
        DO J=1,NY
            DO I=1,NX
                IF(P(I,J) .GT. Plim)THEN         ! If to high pressure
                    !WRITE(4,*)'BAD P(I,J)', P(I,J), 'For I J = ', I ,J , 'H(I,J) = ', H(I,J)
                    P(I,J)=Plim                 ! Extra check added by Carl 
                    ink=1                       ! Set ink to 1 to later increase H00 to release the contact pressure. 
                ENDIF
            ENDDO
        ENDDO
        
        ! Calculate the elastic deformation --------------------------------------------------
        CALL VI
        ! Calculate the elastic deformation from the sidepressures if a cylinder
        IF( ref .EQ. 1 .OR. ref .EQ. 2 .OR. ref .EQ. 5 .OR. ref .EQ. 6) THEN
            CALL VIside
        ELSE
            Wside=0.0*W
        ENDIF

    ! If using equation for asperity. Then this gives the position. 
    position=X0+(XE-X0)*t/Ntime ! = Ua/US*t/Ntime*Tend
        
    ! Choice of asperity chape ------------------------------------------------------------------------------------------------ 
    IF( ref .EQ. 1) THEN
        ! Newtonian acc X.Tans paper ref 30 31, Cylinder
        DO  J=1,NN
            DO  I=1,NX 
                asp=asph*(10**(-10*(X(I)-position)**2))!*(10**(-5*(Y(J))**2))
                bump=asp*cos(2*PAI*(X(I)-position))!*cos(2*PAI*Y(J))
                
                ! + or - depending on wake or asperity
                H0=RAD(I,J)+W(I,J)+bump+Wside(I,J)!-Emat(I+NX+20-t,J)
                IF(H0.LT.HMIN)HMIN=H0
                H(I,J)=H0    
            ENDDO
        ENDDO
        
    ELSE IF( ref .EQ. 3)THEN
        ! Input parameters ac Holmes et al Transient EHL point contact analysis 2003, Ball 
        DO  J=1,NN
            DO  I=1,NX 
                Ab=asph
                Wb=aspw
                bump=Ab*10**(-10*((X(I)-position)/Wb)**2)*cos(2*PAI*(X(I)-position)/Wb)
                
                ! + or - depending on wake or asperity
                H0=RAD(I,J)+W(I,J)-bump+Wside(I,J)!-Emat(I+NX+20-t,J)
                IF(H0.LT.HMIN)HMIN=H0
                H(I,J)=H0    
            ENDDO
        ENDDO 
        
    ELSE IF(ref .EQ. 6) THEN   ! Dave asperity as line for ref 6
        DO  J=1,NN
            DO  I=1,NX 
                
                bump=0
                IF ( (X(I)-position)**2/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(X(I)-position)**2+asph/2;
                ENDIF
                
                ! + or - depending on wake or asperity
                H0=RAD(I,J)+W(I,J)-bump+Wside(I,J)!-Emat(I+NX+20-t,J)
                IF(H0.LT.HMIN)HMIN=H0
                H(I,J)=H0    
            ENDDO
        ENDDO   
        
    ELSE    ! Dave asperity for ref 2 and 5. And for ref 4 which does not have an own asperity
        DO  J=1,NN
            DO  I=1,NX 
                
                bump=0
                IF ( ((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*PAI/(aspw))+asph/2;
                ENDIF
                
                ! + or - depending on wake or asperity
                H0=RAD(I,J)+W(I,J)-bump+Wside(I,J)!-Emat(I+NX+20-t,J)
                IF(H0.LT.HMIN)HMIN=H0
                H(I,J)=H0    
            ENDDO
        ENDDO   
    !Else If using asperity not shaped by equations
        !Emat(1:2*NX+20,1:NY)=0.0
        ! Back part of E
        !Emat(NX,NY/2-15:NY/2+15)=asph/2
        !Emat(NX+1:NX+3,NY/2-15:NY/2+15)=asph
        !Emat(NX+4,NY/2-15:NY/2+15)=asph/2
        ! Middel bar of E
        !Emat(NX+13,NY/2-3:NY/2+3)=asph/2
        !Emat(NX+4:NX+12,NY/2-3)=asph/2
        !Emat(NX+4:NX+12,NY/2-2:NY/2+2)=asph
        !Emat(NX+4:NX+12,NY/2+3)=asph/2
        ! Edge bar
        !Emat(NX+16,NY/2-16:NY/2-9)=asph/2
        !Emat(NX:NX+15,NY/2-16)=asph/2
        !Emat(NX+4:NX+15,NY/2-9)=asph/2
        !Emat(NX+4:NX+15,NY/2-15:NY/2-10)=asph
        ! Edge bar
        !Emat(NX+16,NY/2+9:NY/2+16)=asph/2
        !Emat(NX:NX+15,NY/2+16)=asph/2
        !Emat(NX+4:NX+15,NY/2+9)=asph/2
        !Emat(NX+4:NX+15,NY/2+10:NY/2+15)=asph
        !---------------------------------------------------------------
    ENDIF
    
    deformation=0
    DO J=1,NN
        DO I=1,NX
                H(I,J)=H00+H(I,J)
                ! If Contat
                IF (H(I,J) .LT. 5E-6) THEN
                    IF (I .GE. 2) THEN
                        P(I,J)=P(I,J)-1.0*(H(I,J)-5E-6)/(AK(0,0)*PAIAK*DX)
                        P(I,J+1)=P(I,J+1)-0.0*(H(I,J)-5E-6)/(AK(1,0)*PAIAK*DX)
                    ENDIF
                    H(I,J)=5E-6
                    deformation=1
                ENDIF
                IF(P(I,J).LT.0.0)P(I,J)=0.0
                IF(P(I,J).GT.20.0)P(I,J)=20.0
        ENDDO
    ENDDO
    
    IF( deformation .EQ. 1) GOTO 67
    
    ! Calculate the viscocity, the density and the dimentionles parameter EPS -----------------------------------------------------------------------------------------
    IF( ref .LT. 3 .OR. ref .EQ. 6) THEN
        ! Newtonian acc X.Tans paper ref 30 31, Cylinder
        DO J=1,NN
            DO I=1,NX

                ! Roelands equation acc X.Tan
                EDAx(I,J)=EXP(alpha*Pref/Z*(-1+(1+P(I,J)*PH/Pref)**Z)) 
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
            
                ! D-H Formulation acc X.Tan
                 !RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                IF (EDAx(I,J) .LE. 0.0) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'tau = ', tau, 'EDA1 = ', EDA1
                    EDAx(I,J)=0.1
                ENDIF
                IF (EDAy(I,J) .LE. 0.0) THEN
                    WRITE(4,*)'BAD EDAy', EDAy(I,J), 'For I J = ', I ,J, 'tau = ', tau, 'EDA1 = ', EDA1
                    EDAy(I,J)=0.1
                ENDIF
                
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
        ENDDO
    ENDDO

    ELSE IF( ref .EQ. 3)THEN
        ! Newtonian input parameters ac Holmes et al Transient EHL point contact analysis 2003, Ball 
        DO J=1,NN
            DO I=1,NX

                ! Roelands equation acc Holmes
                 EDAx(I,J)=EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1)) 
                 EDAy(I,J)=EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1)) 
                
                ! D-H Formulation acc Homes et al
                 !RO(I,J)=(1 + gH*PH*P(I,J))/(1+lH*P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                 IF (EDAx(I,J) .LE. 0.0) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'tau = ', tau, 'EDA1 = ', EDA1
                    EDAx(I,J)=0.1
                ENDIF
                IF (EDAy(I,J) .LE. 0.0) THEN
                    WRITE(4,*)'BAD EDAy', EDAy(I,J), 'For I J = ', I ,J, 'tau = ', tau, 'EDA1 = ', EDA1
                    EDAy(I,J)=0.1
                ENDIF
                
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
        ENDDO
    ENDDO
    
    ELSE !Nonnewtoinan imput for case 4 and 5
        DO J=1,NN
            DO I=1,NX
           
                ! Yasutomi equations
                Tg(I,J)=Tg0+YA1*LOG(1+YA2*PH*P(I,J))
                pg(I,J)=1/YA2*(EXP((1/YA1)*(Temp(I,J)-Tg0))-1)
                YF(I,J)=1-YB1*LOG(1+YB2*PH*P(I,J))
                
                IF (P(I,J)*PH .GT. pg(I,J)) THEN
                    EDA1=Yedag*EXP(Yalfag*(P(I,J)*PH-pg(I,J)))
                ELSE
                    EDA1=Yedag*10**(-(YC1*(Temp(I,J)-Tg(I,J))*YF(I,J))/(YC2+(Temp(I,J)-Tg(I,J))*YF(I,J)))
                ENDIF

                tau=LOG(EDA1/H(I,J)/tauc+SQRT((EDA1/H(I,J)/tauc)**2+1))                 ! Absolute value since does not matter which direction the sliding is. tauc includes abs(a**2/(R*(Ua-Ub)))

                IF (tau .EQ. 0.0) Then
                    EDAx(I,J)=EDA1
                    EDAy(I,J)=EDA1
                ELSE
                    EDAx(I,J)=EDA1/(COSH(tau))                                          ! Rewritten to get better numbers. !EDA1/(f_tau+tau*df_tau)
                    EDAy(I,J)=EDA1*tau/(SINH(tau))                                      ! f_tau
                ENDIF
                
                IF (EDAx(I,J) .LE. 0.0) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'tau = ', tau, 'EDA1 = ', EDA1
                    EDAx(I,J)=0.1
                ENDIF
                IF (EDAy(I,J) .LE. 0.0) THEN
                    WRITE(4,*)'BAD EDAy', EDAy(I,J), 'For I J = ', I ,J, 'tau = ', tau, 'EDA1 = ', EDA1
                    EDAy(I,J)=0.1
                ENDIF
                
                
                ! D-H Formulation acc P.Ehret D.Dowsin and C.M. Taylor
                RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                xi(I,J)=taua*P(I,J)+taua2*P(I,J)**2 
                IF(xi(I,J) .GT. xilim) xi(I,J)=xilim
                IF(xi(I,J) .LT. 0.0) xi(I,J)=0.0
                
                
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
        ENDDO
    ENDDO
    ENDIF
    
    
    
    ! Mirror all updated parameters around the X-axis------------------------------------------------------------------------------------------------------------------
 67     DO J=NN+1,NY
            JJ=NY-J+1
            DO I=1,NX
                P(I,J)=P(I,JJ)
                H(I,J)=H(I,JJ)
                RO(I,J)=RO(I,JJ)
                xi(I,J)=xi(I,JJ)
                EDAx(I,J)=EDAx(I,JJ)
                EDAy(I,J)=EDAy(I,JJ)
                EPSx(I,J)=EPSx(I,JJ)
    	        EPSy(I,J)=EPSy(I,JJ)
            ENDDO
        ENDDO
        
            
         
        IF(PK .GT. 1 .AND. PK .LT. 10) GOTO 15
        
    
    RETURN
    END
