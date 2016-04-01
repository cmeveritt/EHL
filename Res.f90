 ! Subrutine for itteration according to Reynolds equation of the timedependent solution
	SUBROUTINE Res(D2, D2p, MM, J0,J1, JJ, J,C1,C2,kvot, DT2)
    	implicit none
        COMMON      /CTRA4/ D, A
        COMMON      /ITERp/ ID
        COMMON      /Grid2/DX1,DX2,DX3,DX4
        COMMON      /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
        COMMON      /NonNew/ tauc, taua, taua2,xilim                               ! Non Newtonian viscosity parameters
        COMMON      /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag    ! Lubrication parameters acc Yasutomi
        COMMON      /Rho/RA1,RA2                                             ! Density parameters
        COMMON      /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
        COMMON      /Method/term1,term2,term3,term4                          ! Controling the numerical method
        COMMON      /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referenseCOMMON      /HREEp/Emat, pg, Tg,YF
        COMMON      /COMAK2D/AK(0:300,0:300)                                       ! Parameters for HREE subrutine
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast  ! Past
        COMMON      /Residual/A5dxx,A5dx,A5dt                                   ! Parameters for the residual
        COMMON      /AKparam/ AK00, AK10, AK20, AK30, BK00, BK10, BK20, BK30
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        real         A5dxx(1:300,1:300),A5dx(1:300,1:300),A5dt(1:300,1:300)
        real         D(1:300), A(1:5,1:300), ID(1:300)
        real        AK, RA1, RA2
        Integer     NX, NY, NN, ref, t, tmeth
        integer     term1, term2, term3, term4
        integer     MM,K, J,JJ,J1,J0, I, I0,I1, I2, IA, II
        Real        PAI
        Real        Z,  EDA0, Pref, alpha
        Real        DX, X0, XE
        real        ENDA, A1, A2, A3, ALFA, G, AHM, HM0, UTL, HM0r, PH, H00
        real        tauc_real, tauc, tauS, taua, taua2,xilim
        real        Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
        real        kH, xH, gH, lH, L, M
        real        kvot, DT2
        real        dxidP, dxidPp, dROdP, dROdPp, dHdx, dHdxp, dPdx, dPdxp, d2ROdP2, d2xidP2 
        real        K00, K10, K20, K30, PAIAK
        real        DX1,DX2,DX3,DX4
        real        D1, D2, D3,D4, D5,D1p, D2p, D3p,D4p, D5p
        real        P1,P2,P3,P4,P5,P1p,P2p,P3p,P4p,P5p
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
        
                ! Uppdate the presure for NX=2 untill there is zero pressure near the exit. at NX=1 BCs restric this pressure to be 0
                DO 50 I=2,MM
                    I0=I-1
                    I2=I-2
                    IF (I2 .EQ. 0) I2=I0
                    I1=I+1
                    II=I
                    
                    ! For current time step                    
                    D1=D2
                    D2=0.5*(EPSx(I1,J)+EPSx(I,J))
                    D4=0.5*(EPSy(I,J0)+EPSy(I,J))
                    D5=0.5*(EPSy(I,J1)+EPSy(I,J))
                    P1=P(I0,JJ)
                    P2=P(I1,JJ)
                    P3=P(I,JJ)
                    P4=P(I,JJ+1)
                    P5=P(I,JJ-1)
                    D3=D1+D2+D4+D5
                    
                    dHdx=(term4*H(I1,JJ)+term3*H(I,JJ)+term2*H(I0,JJ)+term1*H(I2,JJ))/(2*DX)       
                    dPdx=(term4*P(I1,JJ)+term3*P(I,JJ)+term2*P(I0,JJ)+term1*P(I2,JJ))/(2*DX)
                    IF (xi(I,J) .GE. 0.95*xilim )Then               !.OR. xi(I,J) .LT. 0.01*xilim ej inkluderat för att xi har en derivata även när xi=0
                        dxidP=0.0
                        d2xidP2=0.0
                    ELSE
                        dxidP=taua+2*taua2*P(I,J)
                        d2xidP2=2*taua2
                    ENDIF
                    dROdP=RA1*PH/(RA2*P(I,JJ)*PH+1)**2
                    d2ROdP2=-2*RA1*PH*RA2*PH/(RA2*PH*P(I,JJ)+1)**3
                    
                    ! For previus timestep
                    D1p=D2p
                    D2p=0.5*(EPSxpast(I1,J)+EPSxpast(I,J))
                    D4p=0.5*(EPSypast(I,J0)+EPSypast(I,J))
                    D5p=0.5*(EPSypast(I,J1)+EPSypast(I,J))
                    P1p=Ppast(I0,JJ)
                    P2p=Ppast(I1,JJ)
                    P3p=Ppast(I,JJ)
                    P4p=Ppast(I,JJ+1)
                    P5p=Ppast(I,JJ-1)
                    D3p=D1p+D2p+D4p+D5p
                    
                    dHdxp=(term4*Hpast(I1,JJ)+term3*Hpast(I,JJ)+term2*Hpast(I0,JJ)+term1*Hpast(I2,JJ))/(2*DX)       
                    dPdxp=(term4*Ppast(I1,JJ)+term3*Ppast(I,JJ)+term2*Ppast(I0,JJ)+term1*Ppast(I2,JJ))/(2*DX)
                    IF (xipast(I,JJ) .GE. 0.95*xilim)Then
                        dxidPp=0.0
                    ELSE
                        dxidPp=taua+2*taua2*P(I,JJ)
                    ENDIF
                    dROdPp=RA1*PH/(RA2*Ppast(I,JJ)*PH+1)**2

                    ! Start calculations
                    
                    IF(J.EQ.NN.AND.ID(I).EQ.1)P(I,J)=P(I,J)-0.5*C2*D(I) ! If at the central node. Update the pressure before calculating
                    
                    ! If the lubrication height H is below zeroe. This is a built in fix
                    IF(H(I,J).LE.0.0)THEN
                        ID(I)=2
                        A(1,II)=0.0
                        A(2,II)=0.0
                        A(3,II)=1.0
                        A(4,II)=0.0
                        A(5,II)=1.0
                        A(1,II-1)=0.0
                        GOTO 50
                    ENDIF
                    
                    ! Calculate the residual of the node based on Reynolds equation
                    
                    ! d/dx ( eps * d/dx( P))
                    AII5dxx   = kvot*      (D1*P1+D2*P2+D4*P4+D5*P5-D3*P3)
                    AII5dxxp  = (1.0-kvot)*(D1p*P1p+D2p*P2p+D4p*P4p+D5p*P5p-D3p*P3p)
                    
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
                    
                    ! d(ROH)/dt 
                    AII5dt1   = 2*RO(I,J)*(H(I,J)-Hpast(I,J))
                    AII5dt2  =  2*(dROdP+dROdPp)/2*(P(I,J)-Ppast(I,J))*H(I,J)
                    
                    ! A(II+5) = Equation (14.11) rewritten to -L(i,j)=gamma in Equation (13.18) + time dep term d(RO*H)/dT
                    A5dxx(I,J)  =   -DX3*      ( AII5dxx  +   AII5dxxp)                                             
                    A5dx(I,J)   =             ( AII5dx   +   AII5dxp )                     
                    A5dt(I,J)   =   DT2*       ( AII5dt1   +   AII5dt2 )
                    A(5,II)     =   A5dxx(I,J)  +A5dx(I,J)  +A5dt(I,J)
                    
                    ! Depending on the parameters in the 2nd order derivative term, choose which method to use for the pressureuppdiating. 
                    IF(D1.GE.DX4)GOTO 30
                    IF(D2.GE.DX4)GOTO 30
                    IF(D4.GE.DX4)GOTO 30
                    IF(D5.GE.DX4)GOTO 30
                    ID(I)=1
                    IF(J.EQ.NN)P5=P4
                    IF(J.EQ.NN)P5p=P4p
                    
                    ! Methode 1 -----------------------------------
                    K00=BK00
                    K10=BK10
                    K20=BK20
                    K30=BK30
                    
                    ! d(kvot*(D1*P1+D2*P2+D4*P4+D5*P5-D3*P3))/dP(*,j)
                    AII2dxx =   D1+0.25*D3
                    AII3dxx =   -1.25*D3
                    AII4dxx =   D2+0.25*D3
                    
                    GOTO 40
                    ! Methode 2 -----------------------------------
                    
    30              ID(I)=0
                    P4=P(I,J0)
                    P4p=Ppast(I,J0)
                    IF(J.EQ.NN)P5=P4
                    IF(J.EQ.NN)P5p=P4p
                    
                    K00=AK00
                    K10=AK10
                    K20=AK20
                    K30=AK30
                    
                    ! d(kvot*(D1*P1+D2*P2+D4*P4+D5*P5-D3*P3))/dP(*,j)
                    AII2dxx =   D1
                    AII3dxx =   -D3
                    AII4dxx =   D2
            
                    !------------------------------------------------------------------------------------------------
                    ! Pressure derivatives of the residual
                    ! The 1/2 is added later in the derivatives of P and H
                    ! -d((dH/dx*(1-xi)*RO)/dP(*,J)                                          
40                  AII1dx1  =   -PAIAK*DX*(term4*K30 + term3*K20 + term2*K10 + term1*K00)*(1-xi(I,J))*RO(I,J)/(2*DX)
                    AII2dx1  =   -PAIAK*DX*(term4*K20 + term3*K10 + term2*K00 + term1*K10)*(1-xi(I,J))*RO(I,J)/(2*DX)
                    AII3dx1  =   -PAIAK*DX*(term4*K10 + term3*K00 + term2*K10 + term1*K20)*(1-xi(I,J))*RO(I,J)/(2*DX)       - dHdx*(-dxidP)*RO(I,J)     - (1-xi(I,J))*dROdP*dHdx
                    AII4dx1  =   -PAIAK*DX*(term4*K00 + term3*K10 + term2*K20 + term1*K30)*(1-xi(I,J))*RO(I,J)/(2*DX)   
                    
                    ! d(H(I,J)*RO(I,J)*dxidP*dPdx)/dP(*,j)
                    AII1dx2  =   PAIAK*DX*K20*RO(I,J)*dxidP*dPdx                                                                            + H(I,J)*RO(I,J)*dxidP*(term1)/(2*DX)
                    AII2dx2  =   PAIAK*DX*K10*RO(I,J)*dxidP*dPdx                                                                            + H(I,J)*RO(I,J)*dxidP*(term2)/(2*DX)
                    AII3dx2  =   PAIAK*DX*K00*RO(I,J)*dxidP*dPdx        + H(I,J)*dROdP*dxidP*dPdx       + H(I,J)*RO(I,J)*d2xidP2*dPdx       + H(I,J)*RO(I,J)*dxidP*(term3)/(2*DX)
                    AII4dx2  =   PAIAK*DX*K10*RO(I,J)*dxidP*dPdx                                                                            + H(I,J)*RO(I,J)*dxidP*(term4)/(2*DX)
                    
                    ! -d(H(I,J)*(1-xi(I,J))*dROdP*dPdx)/dP(*,j)
                    AII1dx3  =  - PAIAK*DX*(1-xi(I,J))*K20*dROdP*dPdx                                                                       - H(I,J)*(1-xi(I,J))*dROdP*(term1)/(2*DX)
                    AII2dx3  =  - PAIAK*DX*(1-xi(I,J))*K10*dROdP*dPdx                                                                       - H(I,J)*(1-xi(I,J))*dROdP*(term2)/(2*DX)
                    AII3dx3  =  - PAIAK*DX*(1-xi(I,J))*K00*dROdP*dPdx   + H(I,J)*dxidP*dROdP*dPdx       -H(I,J)*(1-xi(I,J))*d2ROdP2*dPdx    - H(I,J)*(1-xi(I,J))*dROdP*(term3)/(2*DX)
                    AII4dx3  =  - PAIAK*DX*(1-xi(I,J))*K10*dROdP*dPdx                                                                       - H(I,J)*(1-xi(I,J))*dROdP*(term4)/(2*DX)
                    
                    ! -d(d(2*RO(I,J)*H(I,J))/dt)/dP(*,j) = -d( 2*RO(I,J)*(H(I,J)-Hpast(I,J))*DT2 +  2*(dROdP+dROdPp)/2*(P(I,J)-Ppast(I,J))*DT2*H(I,J) )/dP(*,j)
                    AII1dt  =   -2* (RO(I,J)*(PAIAK*DX*K20)  +                                                                               (dROdP+dROdPp)/2*P(I,J)*(PAIAK*DX*K20)   )    *DT2 ! DT2 contains a 2
                    AII2dt  =   -2* (RO(I,J)*(PAIAK*DX*K10)  +                                                                               (dROdP+dROdPp)/2*P(I,J)*(PAIAK*DX*K10)   )    *DT2 
                    AII3dt  =   -2* (RO(I,J)*(PAIAK*DX*K00)  +  dROdP*H(I,J)  + d2ROdP2*P(I,J)*H(I,J)/2 + (dROdP+dROdPp)/2*1.0*H(I,J)    +   (dROdP+dROdPp)/2*P(I,J)*(PAIAK*DX*K00)   )    *DT2 
                    AII4dt  =   -2* (RO(I,J)*(PAIAK*DX*K10)  +                                                                               (dROdP+dROdPp)/2*P(I,J)*(PAIAK*DX*K10)   )    *DT2 
                    
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