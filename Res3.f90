 ! Subrutine for finding the residual of the Reynolds equation. Pressure increments are then calculated based on this residual and the evaluated derivatives due to pressure changes 
    ! Res 1-3 are based upon that the Reynolds equation can be adjusted if metal contact occures. These are howevere not fully developed yet. 
	SUBROUTINE Res3(I, D2, D2p, MM, NN, J0, J1, JJ, J, C1, C2, kvot, DT2, SS, Pcenter)
    	implicit none
        COMMON      /CTRA4/ D, A
        COMMON      /ITERp/ ID
        COMMON      /Grid2/DX1,DX2,DX3,DX4
        COMMON      /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0                           ! Lubrication parameters
        COMMON      /NonNew/ tauc, taua, taua2,xilim , xi_param                                 ! Non Newtonian viscosity parameters
        COMMON      /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag                       ! Lubrication parameters acc Yasutomi
        COMMON      /Rho/RA1,RA2                                                                ! Density parameters
        COMMON      /Grid/NX,NY,X0,XE,DX                                                        ! Grid parameters
        COMMON      /Ref/PAIAK, meth, tmeth, Geom, Lub_param, asp_shape, contact_alg            ! Coise of equations based on referense     /HREEp/Emat, pg, Tg,YF
        COMMON      /COMAK2D/AK(0:601,0:601)                                                    ! Parameters for HREE subrutine
        COMMON      /Current/RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                                  ! Current timestep
        COMMON      /CurrentP/P
        COMMON      /CurrentH/H
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast   ! Past
        COMMON      /Past2/Ppast2,Hpast2                                                        ! the timestep two time units ago. 
        COMMON      /Residual/A5dxx,A5dx,A5dt                                                   ! Parameters for the residual
        COMMON      /AKparam/ AK00, AK10, AK20, AK30, BK00, BK10, BK20, BK30
        real         P(1:601,1:601),H(1:601,1:601),RO(1:601,1:601),EPSx(1:601,1:601),EPSy(1:601,1:601),EDAx(1:601,1:601),EDAy(1:601,1:601),xi(1:601,1:601),W(1:601,1:601),Wside(1:601,1:601)                      ! Current timestep
        real         Ppast(1:601,1:601),Hpast(1:601,1:601),ROpast(1:601,1:601),EPSxpast(1:601,1:601),EPSypast(1:601,1:601),EDAxpast(1:601,1:601),EDAypast(1:601,1:601),xipast(1:601,1:601),Wpast(1:601,1:601)  ! Past
        real         A5dxx(1:601,1:601),A5dx(1:601,1:601),A5dt(1:601,1:601)
        real         Ppast2(1:601,1:601),Hpast2(1:601,1:601)
        real         D(1:601), A(1:5,1:601), ID(1:601)
        real        Pcenter(NX,2)
        real        AK, RA1, RA2
        Integer     NX, NY, NN, t, tmeth, SS, Geom, Lub_param, asp_shape, meth, contact_alg
        integer     term1, term2, term3, term4
        integer     MM,K, J,JJ,J1,J0, I, I0,I1, I2, IA, II
        Real        PAI
        Real        Z,  EDA0, Pref, alpha
        Real        DX, X0, XE
        real        ENDA, A1, A2, A3, ALFA, G, AHM, HM0, UTL, HM0r, PH, H00
        real        tauc_real, tauc, tauS, taua, taua2,xilim, xi_param
        real        Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
        real        kH, xH, gH, lH, L, M
        real        kvot, DT2
        real        dxidP, dxidPp, dROdP, dROdPp, dHdx, dHdxp, dPdx, dPdxp, d2ROdP2, d2xidP2 
        real        K00, K10, K20, K30, PAIAK
        real        DX1,DX2,DX3,DX4
        real        D1, D2, D3,D4, D5,D1p, D2p, D3p,D4p, D5p
        real        P1,P2,P3,P4,P5,P1p,P2p,P3p,P4p,P5p, P6, P6p
        real        AK00,AK10,AK20,AK30,BK00,BK10,BK20,BK30
        real        AII5dxx, AII5dxxp, AII5dx, AII5dxp, AII5dx1, AII5dx1p, AII5dx2, AII5dx2p, AII5dx3, AII5dx3p, AII5dt1, AII5dt2
        real        AII4dxx, AII4dx1, AII4dx2, AII4dx3, AII4dx4, AII4dt
        real        AII3dxx, AII3dx1, AII3dx2, AII3dx3, AII3dx4, AII3dt
        real        AII2dxx, AII2dx1, AII2dx2, AII2dx3, AII2dx4, AII2dt
        real                AII1dx1, AII1dx2, AII1dx3, AII1dx4, AII1dt
        real       C1, C2
        SAVE      /Current/                                     ! Current timestep
        SAVE      /Residual/                                    ! Parameters for the residual
        SAVE      /CTRA4/
        SAVE      /ITERp/
        
        ! 1st dx
        term1   =  0
        term2   = -2
        term3   =  2
        term4   =  0
        
                ! Uppdate the presure for NX=2 untill there is zero pressure near the exit. at NX=1 BCs restric this pressure to be 0
               
                    I0=I-1*SS
                    I2=I-2*SS
                    IF (I2 .LE. 0) I2=I0
                    I1=I+1*SS
                    II=I
                                        
                    ! Al nodes are treated the same. 
                    ! IF(J.EQ.NN.AND.ID(I).EQ.1)P(I,J)=P(I,J)-0.5*C2*D(I) ! If at the central node. Update the pressure before calculating
                    
                    
                    ! For current time step -------------------------------------------------------------------------------                   
                    D1=D2
                    D2=0.5*(EPSx(I1,J)+EPSx(I,J))
                    D4=0.5*(EPSy(I,J0)+EPSy(I,J))
                    D5=0.5*(EPSy(I,J1)+EPSy(I,J))
                    P1=P(I0,JJ)
                    P2=P(I1,JJ)
                    P3=P(I,JJ)
                    P4=P(I,JJ+SS)               ! Opposite sign here on SS since where on the JJ side of the mirror
                    P5=P(I,JJ-SS)
                    P6=P(I2,JJ)
                    D3=D1+D2+D4+D5
                    
                    If (JJ .EQ. NN+SS) P5=Pcenter(I,2)
                    If (JJ .EQ. NN) THEN
                        P1=Pcenter(I0,2)
                        P2=Pcenter(I1,2)
                        P3=Pcenter(I,2)
                        P5=Pcenter(I,1)
                        P6=Pcenter(I2,2)
                    ENDIF

                    ! For previus timestep. Here the updating of the current pressures does not matter. The past pressure is fixed. -------------------------------
                    D1p=D2p
                    D2p=0.5*(EPSxpast(I1,J)+EPSxpast(I,J))
                    D4p=0.5*(EPSypast(I,J0)+EPSypast(I,J))
                    D5p=0.5*(EPSypast(I,J1)+EPSypast(I,J))
                    P1p=Ppast(I0,JJ)
                    P2p=Ppast(I1,JJ)
                    P3p=Ppast(I,JJ)
                    P4p=Ppast(I,JJ+SS)
                    P5p=Ppast(I,JJ-SS)
                    P6p=Ppast(I2,JJ)
                    D3p=D1p+D2p+D4p+D5p
                    
                    ! Depending on the parameters in the 2nd order derivative term, choose which method to use for the pressureuppdiating. 
                    IF(D1.GE.DX4)GOTO 30
                    IF(D2.GE.DX4)GOTO 30
                    IF(D4.GE.DX4)GOTO 30
                    IF(D5.GE.DX4)GOTO 30
                    ID(I)=1
                    ! Same for all nodes
                    !IF(J.EQ.NN)P5=P4
                    !IF(J.EQ.NN)P5p=P4p
                    
                    ! Methode 1 -----------------------------------
                    K00=BK00
                    K10=BK10
                    K20=BK20
                    K30=BK30
                    
                    ! d(kvot*(D1*P1+D2*P2+D4*P4+D5*P5-D3*P3))/dP(*,j)
                    AII2dxx =  0.0! D1+0.25*D3
                    AII3dxx =  0.0! -1.25*D3
                    AII4dxx =  0.0! D2+0.25*D3
                    
                    GOTO 40
                    ! Methode 2 -----------------------------------
30                  ID(I)=0
                    ! Use updated values forP4 and P5. 
                    P4=P(I,J0)
                    P4p=Ppast(I,J0)
                    !IF(J.EQ.NN)P5=P4
                    !IF(J.EQ.NN)P5p=P4p
                    
                    K00=AK00
                    K10=AK10
                    K20=AK20
                    K30=AK30
                    
                    ! d(kvot*(D1*P1+D2*P2+D4*P4+D5*P5-D3*P3))/dP(*,j)
                    AII2dxx =   0.0!D1
                    AII3dxx =   0.0!-D3
                    AII4dxx =   0.0!D2
                    
                    ! Derivatives --------------------------------------------------------------------------------------------------------------------
40                  dHdx=(term4*H(I1,JJ)+term3*H(I,JJ)+term2*H(I0,JJ)+term1*H(I2,JJ))/(2*DX*SS)       
                    dPdx=(term4*P2+term3*P3+term2*P1+term1*P6)/(2*DX*SS)
                    IF (xi(I,J) .GE. 0.95*xilim )Then               !.OR. xi(I,J) .LT. 0.01*xilim ej inkluderat för att xi har en derivata även när xi=0
                        dxidP=0.0
                        d2xidP2=0.0
                    ELSE
                        dxidP=taua+2*taua2*P(I,J)
                        d2xidP2=2*taua2
                    ENDIF
                    dROdP=RA1*PH/(RA2*P(I,JJ)*PH+1)**2
                    d2ROdP2=-2*RA1*PH*RA2*PH/(RA2*PH*P(I,JJ)+1)**3
                    
                    dHdxp=(term4*Hpast(I1,JJ)+term3*Hpast(I,JJ)+term2*Hpast(I0,JJ)+term1*Hpast(I2,JJ))/(2*DX*SS)       
                    dPdxp=(term4*P2p+term3*P3p+term2*P1p+term1*P6p)/(2*DX*SS)
                    IF (xipast(I,JJ) .GE. 0.95*xilim)Then
                        dxidPp=0.0
                    ELSE
                        dxidPp=taua+2*taua2*Ppast(I,JJ)
                    ENDIF
                    dROdPp=RA1*PH/(RA2*Ppast(I,JJ)*PH+1)**2

                    ! Start calculations ---------------------------------------------------------------------------------------------------------------------
                    
                    ! Calculate the residual of the node based on Reynolds equation
                    
                    ! d/dx ( eps * d/dx( P))
                    AII5dxx   = 0.0!kvot*      (D1*P1+D2*P2+D4*P4+D5*P5-D3*P3)
                    AII5dxxp  = 0.0!(1.0-kvot)*(D1p*P1p+D2p*P2p+D4p*P4p+D5p*P5p-D3p*P3p)
                    
                    ! d((1-xi)*ROH)/dx = dH/dx*(1-xi)*RO-H*RO*dxidP*dPdx+H*(1-xi)*dRdP*DPdx        at crrent time
                    AII5dx1  =  dHdx*       (1-xi(I,J))*    RO(I,J)
                    AII5dx2  =  H(I,J)*    (-dxidP*dPdx)*   RO(I,J)
                    AII5dx3  =  H(I,J)*     (1-xi(I,J))*    dROdP*dPdx
                    
                    ! d((1-xi)*ROH)/dx at previous time
                    AII5dx1p  =  dHdxp*         (1-xipast(I,J))*    ROpast(I,J)
                    AII5dx2p  =  Hpast(I,J)*    (-dxidPp*dPdxp)*    ROpast(I,J)
                    AII5dx3p  =  Hpast(I,J)*    (1-xipast(I,J))*    dROdPp*dPdxp
                    AII5dx   =  kvot*      (AII5dx1     +AII5dx2     +AII5dx3)    
                    AII5dxp  =  (1.0-kvot)*(AII5dx1p    +AII5dx2p    +AII5dx3p )
                    
                    ! d(ROH)/dt is set to 0
                    IF( tmeth .LE. 2) Then
                        AII5dt1   = 0.0!2*RO(I,J)*(H(I,J)-Hpast(I,J))
                        AII5dt2  =  0.0!2*(dROdP+dROdPp)/2*(P(I,JJ)-Ppast(I,JJ))*H(I,J)
                    elseif( tmeth .eq. 3) then
                        AII5dt1   = 0.0!2*RO(I,J)*  (1.5*H(I,J)  -2.0*Hpast(I,J)  +0.5*Hpast2(I,J))
                        AII5dt2  =  0.0!2*dROdP*    (1.5*P(I,JJ) -2.0*Ppast(I,JJ) +0.5*Ppast2(I,J)) *H(I,J)
                    elseIF( tmeth .EQ. 4) Then
                        AII5dt1   = 0.0!2*RO(I,J)*(H(I,J)-Hpast(I,J))
                        AII5dt2  =  0.0!2*dROdP*(P(I,JJ)-Ppast(I,JJ))*H(I,J)
                    endif
                    
                    
                    
                    ! A(II+5) = Equation (14.11) rewritten to -L(i,j)=gamma in Equation (13.18) + time dep term d(RO*H)/dT
                    A5dxx(I,J)  =   0!-DX3*      ( AII5dxx  +   AII5dxxp)                                             
                    A5dx(I,J)   =             ( AII5dx   +   AII5dxp )                     
                    A5dt(I,J)   =   0! DT2*       ( AII5dt1   +   AII5dt2 )
                    A(5,II)     =   A5dxx(I,J)  +A5dx(I,J)  +A5dt(I,J)
                    
                    ! Since contact we assume
                    dROdP       =   0.0
                    dxidP       =   0.0
                    d2ROdP2     =   0.0
                    ! We want dHdx     to be   0.0
                    
                    
                    !------------------------------------------------------------------------------------------------
                    ! Pressure derivatives of the residual
                    ! The 1/2 is added later in the derivatives of P and H
                    ! -d((dH/dx*(1-xi)*RO)/dP(*,J)                                          
                    AII1dx1  =   -PAIAK*DX*SS*(term4*K30 + term3*K20 + term2*K10 + term1*K00)*(1-xi(I,J))*RO(I,J)/(2*DX*SS)
                    AII2dx1  =   -PAIAK*DX*SS*(term4*K20 + term3*K10 + term2*K00 + term1*K10)*(1-xi(I,J))*RO(I,J)/(2*DX*SS)
                    AII3dx1  =   -PAIAK*DX*SS*(term4*K10 + term3*K00 + term2*K10 + term1*K20)*(1-xi(I,J))*RO(I,J)/(2*DX*SS)       - dHdx*(-dxidP)*RO(I,J)     - (1-xi(I,J))*dROdP*dHdx
                    AII4dx1  =   -PAIAK*DX*SS*(term4*K00 + term3*K10 + term2*K20 + term1*K30)*(1-xi(I,J))*RO(I,J)/(2*DX*SS)   
                    
                    ! d(H(I,J)*RO(I,J)*dxidP*dPdx)/dP(*,j)
                    AII1dx2  =   PAIAK*DX*SS*K20*RO(I,J)*dxidP*dPdx                                                                            + H(I,J)*RO(I,J)*dxidP*(term1)/(2*DX*SS)
                    AII2dx2  =   PAIAK*DX*SS*K10*RO(I,J)*dxidP*dPdx                                                                            + H(I,J)*RO(I,J)*dxidP*(term2)/(2*DX*SS)
                    AII3dx2  =   PAIAK*DX*SS*K00*RO(I,J)*dxidP*dPdx        + H(I,J)*dROdP*dxidP*dPdx       + H(I,J)*RO(I,J)*d2xidP2*dPdx       + H(I,J)*RO(I,J)*dxidP*(term3)/(2*DX*SS)
                    AII4dx2  =   PAIAK*DX*SS*K10*RO(I,J)*dxidP*dPdx                                                                            + H(I,J)*RO(I,J)*dxidP*(term4)/(2*DX*SS)
                    
                    ! -d(H(I,J)*(1-xi(I,J))*dROdP*dPdx)/dP(*,j)
                    AII1dx3  =  - PAIAK*DX*SS*(1-xi(I,J))*K20*dROdP*dPdx                                                                       - H(I,J)*(1-xi(I,J))*dROdP*(term1)/(2*DX*SS)
                    AII2dx3  =  - PAIAK*DX*SS*(1-xi(I,J))*K10*dROdP*dPdx                                                                       - H(I,J)*(1-xi(I,J))*dROdP*(term2)/(2*DX*SS)
                    AII3dx3  =  - PAIAK*DX*SS*(1-xi(I,J))*K00*dROdP*dPdx   + H(I,J)*dxidP*dROdP*dPdx       -H(I,J)*(1-xi(I,J))*d2ROdP2*dPdx    - H(I,J)*(1-xi(I,J))*dROdP*(term3)/(2*DX*SS)
                    AII4dx3  =  - PAIAK*DX*SS*(1-xi(I,J))*K10*dROdP*dPdx                                                                       - H(I,J)*(1-xi(I,J))*dROdP*(term4)/(2*DX*SS)
                    
                    if(tmeth .LE. 2) THEN
                    ! -d(d(2*RO(I,J)*H(I,J))/dt)/dP(*,j) = -d( 2*RO(I,J)*(H(I,J)-Hpast(I,J))*DT2 +  2*(dROdP+dROdPp)/2*(P(I,J)-Ppast(I,J))*DT2*H(I,J) )/dP(*,j)
                    AII1dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K20)  +                                                                                (dROdP+dROdPp)/2*P(I,JJ)*(PAIAK*DX*SS*K20)   )    *DT2 ! DT2 contains a 2
                    AII2dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10)  +                                                                                (dROdP+dROdPp)/2*P(I,JJ)*(PAIAK*DX*SS*K10)   )    *DT2 
                    AII3dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K00)  +  dROdP*H(I,J)  + d2ROdP2*P(I,JJ)*H(I,J)/2 + (dROdP+dROdPp)/2*1.0*H(I,J)    +   (dROdP+dROdPp)/2*P(I,JJ)*(PAIAK*DX*SS*K00)   )    *DT2 
                    AII4dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10)  +                                                                                (dROdP+dROdPp)/2*P(I,JJ)*(PAIAK*DX*SS*K10)   )    *DT2 
                    ELSEIF(tmeth .EQ. 3) THEN
                    ! -d(d(2*RO(I,J)*H(I,J))/dt)/dP(*,j) = -d( 2*RO(I,J)*(1.5*H(I,J)  -2.0*Hpast(I,J)  +0.5*Hpast2(I,J))+   2*dROdP*(1.5*P(I,JJ) -2.0*Ppast(I,JJ) +0.5*Ppast2(I,J)) *H(I,J) )/dP(*,j)
                    AII1dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K20*1.5)  +                                                                                                                                                (dROdP)*(1.5*P(I,JJ) -2.0*Ppast(I,J) +0.5*Ppast2(I,J))*(PAIAK*DX*SS*K20)   )    *DT2 ! DT2 contains a 2
                    AII2dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10*1.5)  +                                                                                                                                                (dROdP)*(1.5*P(I,JJ) -2.0*Ppast(I,J) +0.5*Ppast2(I,J))*(PAIAK*DX*SS*K10)   )    *DT2 
                    AII3dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K00*1.5)  +  dROdP*(1.5*H(I,J) -2.0*Hpast(I,J)  +0.5*Hpast2(I,J))  + d2ROdP2*(1.5*P(I,JJ)-2.0*Ppast(I,J)+0.5*Ppast2(I,J))*H(I,J) + (dROdP)*1.0*H(I,J)    + (dROdP)*(1.5*P(I,JJ) -2.0*Ppast(I,J) +0.5*Ppast2(I,J))*(PAIAK*DX*SS*K00)   )    *DT2 
                    AII4dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10*1.5)  +                                                                                                                                                (dROdP)*(1.5*P(I,JJ) -2.0*Ppast(I,J) +0.5*Ppast2(I,J))*(PAIAK*DX*SS*K10)   )    *DT2 
                    elseif(tmeth .EQ. 4) THEN
                    ! -d(d(2*RO(I,J)*H(I,J))/dt)/dP(*,j) = -d( 2*RO(I,J)*(H(I,J)-Hpast(I,J))+ 2*dROdP*(P(I,JJ)-Ppast(I,JJ))*H(I,J) )/dP(*,j)
                    AII1dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K20)  +                                                                                            dROdP*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K20)   )    *DT2 ! DT2 contains a 2
                    AII2dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10)  +                                                                                            dROdP*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K10)   )    *DT2 
                    AII3dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K00)  +  dROdP*(H(I,J)-Hpast(I,J))  + d2ROdP2*(P(I,JJ)-Ppast(I,J))*H(I,J) + dROdP*1.0*H(I,J)  +    dROdP*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K00)   )    *DT2 
                    AII4dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10)  +                                                                                            dROdP*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K10)   )    *DT2                                                                                                         
                    ENDIF
                    
                    ! A(II+1) = dL(i,j)/dP(i-2,j)
                    A(1,II)=                            kvot*(AII1dx1 + AII1dx2 +AII1dx3)  +   AII1dt                                
                    ! A(II+2)= dL(i,j)/dP(i-1,j) 
                    A(2,II) =   kvot*DX3*AII2dxx     +  kvot*(AII2dx1 + AII2dx2 +AII2dx3)  +   AII2dt                   
                    ! A(II+3) = dL(i,j)/dP(i,j)
                    A(3,II) =   kvot*DX3*AII3dxx    +   kvot*(AII3dx1 + AII3dx2 +AII3dx3)  +   AII3dt 
                    ! A(II+4) = dL(i,j)/dP(i+1,j)
                    A(4,II) =   kvot*DX3*AII4dxx    +   kvot*(AII4dx1 + AII4dx2 +AII4dx3)  +   AII4dt      
                    
50              CONTINUE
        RETURN
    END