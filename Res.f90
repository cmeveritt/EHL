 ! Subrutine for finding the residual of the Reynolds equation. Pressure increments are then calculated based on this residual and the evaluated derivatives due to pressure changes 
	SUBROUTINE Res(D2, D2p, MM, NN, J0, J1, JJ, J, kvot,Tcvot, DT2, SS, Pcenter)
    	implicit none
        include     'inc_COMAK2D.h'
        include     'inc_Contact_mat.h'
        include     'inc_CTRA4.h'
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_Grid.h'
        include     'inc_Grid2.h'
        include     'inc_Iterp.h'
        include     'inc_Method.h'
        include     'inc_NonNew.h'
        include     'inc_Outp.h'
        include     'inc_Past.h'
	    include     'inc_Past2.h'
	    include     'inc_PastRO.h'
        include     'inc_Ref.h'
	    include     'inc_Residual.h'
	    include     'inc_Rho.h'
        include     'inc_RLarsson.h'
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'
        real        Pcenter(NX,2)                           !D_IN: /Grid/ -> Nx
        Integer     NN, ss
        integer     MM, J,JJ,J1,J0, I, I0,I1, I2, II
        Real        PAI
        real        kvot, DT2
        real        dxidP, dxidPp, dROdP, dROdPp, dHdx, dHdxp, dPdx, dPdxp, d2ROdP2, d2xidP2 
        real        K00, K10, K20, K30
        real        D1, D2, D3,D4, D5,D1p, D2p, D3p,D4p, D5p
        real        P1,P2,P3,P4,P5,P1p,P2p,P3p,P4p,P5p, P6, P6p
        real        AII5dxx, AII5dxxp, AII5dx, AII5dxp, AII5dx1, AII5dx1p, AII5dx2, AII5dx2p, AII5dx3, AII5dx3p, AII5dt1, AII5dt2
        real        AII4dxx, AII4dx1, AII4dx2, AII4dx3,  AII4dt
        real        AII3dxx, AII3dx1, AII3dx2, AII3dx3,  AII3dt
        real        AII2dxx, AII2dx1, AII2dx2, AII2dx3,  AII2dt
        real        AII1dxx, AII1dx1, AII1dx2, AII1dx3,  AII1dt
        real        Tcvot
        real        test1, test2, test3, test4,  dp2dx, dp2dy
        ! output
        SAVE      /Current/                                     ! Current timestep
        SAVE      /Residual/                                    ! Parameters for the residual
        SAVE      /CTRA4/
        SAVE      /ITERp/
        DATA       PAI/3.14159265/
        
        Ub          =2*um-Ua             !D_IN: /Outp/ -> Um, Ua   !D_OUT: Ub -> /Outp/ (not saved)
                ! Uppdate the presure for NX=2 untill there is zero pressure near the exit. at NX=1 BCs restric this pressure to be 0
                DO 50 I=1+1*SS,MM,SS    ! !!! DO with a label here don't really know what this does
                    I0=I-1*SS
                    I2=I-2*SS
                    IF (I2 .LE. 0) I2=I0
                    I1=I+1*SS
                    II=I

                    ! For current time step -------------------------------------------------------------------------------                   
                    D1=D2
                    D2=0.5*(EPSx(I1,J)+EPSx(I,J))      !D_IN: /Current/ -> EPSx(I,J)
                    D4=0.5*(EPSy(I,J0)+EPSy(I,J))      !D_IN: /Current/ -> EPSy(I,J)
                    D5=0.5*(EPSy(I,J1)+EPSy(I,J))
                    P1=P(I0,JJ)                        !D_IN: /CurrentP/ -> P(I,J)
                    P2=P(I1,JJ)
                    P3=P(I,JJ)
                    P4=P(I,JJ+SS)
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
                    D2p=0.5*(EPSxpast(I1,J)+EPSxpast(I,J))       !D_IN: /Past/ -> EPSxpast(I,J)
                    D4p=0.5*(EPSypast(I,J0)+EPSypast(I,J))       !D_IN: /Past/ -> EPSypast(I,J)
                    D5p=0.5*(EPSypast(I,J1)+EPSypast(I,J))
                    P1p=Ppast(I0,JJ)                             !D_IN: /Past/ -> Ppast(I,J)
                    P2p=Ppast(I1,JJ)
                    P3p=Ppast(I,JJ)
                    P4p=Ppast(I,JJ+SS)
                    P5p=Ppast(I,JJ-SS)
                    P6p=Ppast(I2,JJ)
                    D3p=D1p+D2p+D4p+D5p
                    
                    ! Depending on the parameters in the 2nd order derivative term, choose which method to use for the pressureuppdiating. 
                    IF(D1.GE.DX4)GOTO 30                         !D_IN: /Grid2/ -> DX4            ! !!! Branch point to target 30
                    IF(D2.GE.DX4)GOTO 30
                    IF(D4.GE.DX4)GOTO 30
                    IF(D5.GE.DX4)GOTO 30
                    ID(I)=1                                      !D_OUT: /ITERp/ -> ID(I)
                    
                    ! Methode 1 -----------------------------------
                    ! d(kvot*(D1*P1+D2*P2+D4*P4+D5*P5-D3*P3))/dP(*,j)
                    Call Poiseuille_increment(i,j,i0,j0,i1,j1,i2, K00, K10, K20, K30, D1, D2, D3, 1, Tcvot, ss, AII1dxx, AII2dxx, AII3dxx, AII4dxx, dRodP, d2ROdP2) !CALL: AII1dxx, AII2dxx, AII3dxx, AII4dxx, dRodP, d2ROdP2 are not initialized yet and probably used for output, control = 1
                    GOTO 40        ! !!! Branch point to target 40
                    
                    ! Methode 2 -----------------------------------
30                  ID(I)=0        ! !!! Branch Target 30
                    ! Use updated values for P4. Not for P5 since these values has not been updated yet.  
                    P4=P(I,J0)
                    P4p=Ppast(I,J0)

                    ! d(kvot*(D1*P1+D2*P2+D4*P4+D5*P5-D3*P3))/dP(*,j)
                    Call Poiseuille_increment(i,j,i0,j0,i1,j1,i2, K00, K10, K20, K30, D1, D2, D3, 2, Tcvot, ss, AII1dxx, AII2dxx, AII3dxx, AII4dxx, dRodP, d2ROdP2) !CALL: No comments on arguments. Remember again: the last 6 arguments and KX0s are outputs
                    
                    
                    ! Derivatives --------------------------------------------------------------------------------------------------------------------
40                  dRodPp=dRodPpast(i,j)           ! !!! Branch target 40                                                     !D_IN: /PastRO/ -> dRodPpast
                                                                
                    !if(contact(I,JJ)==1 .or. contact(I1,JJ)==1 .or. contact(I0,JJ)==1)then
                    !    dhdx=(1*H(I1,JJ)+0*H(I,JJ)-1*H(I0,JJ)+0*H(I2,JJ))/(2*DX*SS) !central difference
                    !else
                        dHdx=(term4*H(I1,JJ)+term3*H(I,JJ)+term2*H(I0,JJ)+term1*H(I2,JJ))/(2*DX*SS)                            !D_IN: /Method/ -> term1, term2, term3, term4; /CurrentH/ -> H; /Grid/ -> DX 
                    !endif
                    
                    dPdx=(term4*P2+term3*P3+term2*P1+term1*P6)/(2*DX*SS) 
                    
                    ! xi is a sloidification parameter proposed by P Ehret 1997 On lubricant transport conditions in elastohydrodynamic conjunctions
                    IF(lub_param .EQ. 5 .or. lub_param .eq. 57 .or. lub_param .eq. 58 .or. lub_param .eq. 59) Then
                        IF (xi(I,J) .GE. 0.95*xilim )Then                                                                       !D_IN: /Current/ -> xi
                            dxidP=0.0
                            d2xidP2=0.0
                        ELSE
                            dxidP=taua+2*taua2*P(I,J)                                                                           !D_IN: /NonNew/ -> taua
                            d2xidP2=2*taua2                                                                                     !D_IN: /NonNew/ -> taua2
                        ENDIF
                    elseif(Lub_param .EQ. 51) Then
                        
                        dxidp       =abs(Ua-Ub)/(PAI*Um)    * (taua/((taua*P3)**2+1)) 
                        d2xidP2     =2*taua**3 *P3 / (taua**2*P3**2+1)**2  

                    elseif(Lub_param .EQ. 55) Then
                        
                        dxidp       =abs(Ua-Ub)/(PAI*Um)    * (taua/((taua*dPdx)**2+1)) * (term3)/(2*DX*SS)
                        d2xidP2     =2*taua**3*dPdx/(taua**2*dPdx**2+1)**2              * ((term3)/(2*DX*SS))**2
                    else
                        dxidp       =0.0
                        d2xidP2     =0.0
                    endif
                    
                    dHdxp=(term4*Hpast(I1,JJ)+term3*Hpast(I,JJ)+term2*Hpast(I0,JJ)+term1*Hpast(I2,JJ))/(2*DX*SS)                !D_IN: /Past/ -> Hpast(I,J)
                    dPdxp=(term4*P2p+term3*P3p+term2*P1p+term1*P6p)/(2*DX*SS)

                    IF(lub_param .EQ. 5 .or. lub_param .eq. 57 .or. lub_param .eq. 58 .or. lub_param .eq. 59) Then
                        IF (xipast(I,JJ) .GE. 0.95*xilim)Then                                                                   !D_IN: /Past/ -> xipast(I,J)
                            dxidPp=0.0
                        ELSE
                            dxidPp=taua+2*taua2*Ppast(I,JJ)
                        ENDIF
                    elseif(Lub_param .EQ. 51) Then
                        dxidpp       =abs(Ua-Ub)/(PAI*Um)    * (taua/((taua*P3p)**2+1)) 
                        
                    elseif(Lub_param .EQ. 55) Then
                        
                        dxidPp       =abs(Ua-Ub)/(PAI*Um)    * (taua/((taua*dPdxp)**2+1)) * (term3)/(2*DX*SS)
                    else
                        dxidPp       =0.0

                    endif

                    ! Start calculations ---------------------------------------------------------------------------------------------------------------------
                    ! If the lubrication height H is below zeroe. This is a built-in fix
                    IF(H(I,J).LE.0.0)THEN
                        ID(I)=2
                        A(1,II)=0.0                                                                                             !D_OUT: A -> /CTRA4/
                        A(2,II)=0.0
                        A(3,II)=1.0
                        A(4,II)=0.0
                        A(5,II)=1.0
                        A(1,II-1)=0.0
                        GOTO 50    ! Skip the rest of the derivative calculations
                    ENDIF
                    
                    ! Calculate the residual of the node based on Reynolds equation
                    
                    ! d/dx ( eps * d/dx( P))
                    AII5dxx   = kvot*      (D1*P1+D2*P2+D4*P4+D5*P5-D3*P3)
                    AII5dxxp  = (1.0-kvot)*(D1p*P1p+D2p*P2p+D4p*P4p+D5p*P5p-D3p*P3p)
                    
                    ! d((1-xi)*RO*H)/dx = dH/dx*(1-xi)*RO-H*RO*dxidP*dPdx+H*(1-xi)*dRdP*DPdx        at crrent time
                    AII5dx1  =  dHdx*       (1-xi(I,J))*    RO(I,J)                                                             !D_IN: /Current/ -> RO
                    AII5dx2  =  H(I,J)*    (-dxidP*dPdx)*   RO(I,J)
                    AII5dx3  =  H(I,J)*     (1-xi(I,J))*    dROdP*dPdx
   
                    ! d((1-xi)*ROH)/dx at previous time
                    AII5dx1p  =  dHdxp*         (1-xipast(I,J))*    ROpast(I,J)                                                 !D_IN: /Past/ -> ROpast
                    AII5dx2p  =  Hpast(I,J)*    (-dxidPp*dPdxp)*    ROpast(I,J)
                    AII5dx3p  =  Hpast(I,J)*    (1-xipast(I,J))*    dROdPp*dPdxp
                    AII5dx   =  kvot*      (AII5dx1     +AII5dx2     +AII5dx3)    
                    AII5dxp  =  (1.0-kvot)*(AII5dx1p    +AII5dx2p    +AII5dx3p )
                    
                    ! d(RO*H)/dt 
                    IF( tmeth .LE. 2) Then                                                                                      !D_IN: /Ref/ -> tmeth
                        AII5dt1   = 2*RO(I,J)*(H(I,J)-Hpast(I,J))
                        AII5dt2  =  2*(dROdP+dROdPp)/2*(P(I,JJ)-Ppast(I,JJ))*H(I,J)
                    elseif( tmeth .eq. 3) then
                        AII5dt1   = 2*RO(I,J)*  (1.5*H(I,J)  -2.0*Hpast(I,J)  +0.5*Hpast2(I,J))
                        AII5dt2  =  2*dROdP*    (1.5*P(I,JJ) -2.0*Ppast(I,JJ) +0.5*Ppast2(I,J)) *H(I,J)
                    elseIF( tmeth .EQ. 4) Then
                        AII5dt1   = 2*RO(I,J)*(H(I,J)-Hpast(I,J))
                        AII5dt2  =  2*dROdP*(P(I,JJ)-Ppast(I,JJ))*H(I,J)
                    endif
                    
                    ! A(II+5) = Equation (14.11) rewritten to -L(i,j)=gamma in Equation (13.18) + time dep term d(RO*H)/dT
                    A5dxx(I,J)  =   -DX3*      ( AII5dxx  +   AII5dxxp)                                                                                            !D_IN: /Grid2/ -> DX3                        !D_OUT: /Residual/ -> A5dxx
                    A5dx(I,J)   =             ( AII5dx   +   AII5dxp )                                                                                                                                          !D_OUT: /Residual/ -> A5dx
                    A5dt(I,J)   =   DT2*       ( AII5dt1   +   AII5dt2 )                                                                                                                                        !D_OUT: /Residual/ -> A5dt
                    A(5,II)     =   A5dxx(I,J)  +A5dx(I,J)  +A5dt(I,J)
                    dp2dx = P(i1,j) - 2*P(i,j) + P(i0,j)
                    dp2dy = P(i,j1) - 2*P(i,j) + P(i,j0)
                    

                    If( abs(A(5,II)) .gt. 0.5 ) then ! Added for investigations
                        A(5,II) = sign( 0.5, A(5,II)) 
                    endif

                    !------------------------------------------------------------------------------------------------
                    ! Pressure derivatives of the residual
                    ! The 1/2 is added later in the derivatives of P and H
                    ! -d((dH/dx*(1-xi)*RO)/dP(*,J)                                          
                    AII1dx1  =   -PAIAK*DX*SS*(term4*K30 + term3*K20 + term2*K10 + term1*K00)*(1-xi(I,J))*RO(I,J)/(2*DX*SS)                                        !D_IN: /Ref/ -> PAIAK
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
                    AII1dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K20)  +                                                                                                        (dROdP+dROdPp)/2*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K20)   )    *DT2 ! DT2 contains a 2
                    AII2dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10)  +                                                                                                        (dROdP+dROdPp)/2*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K10)   )    *DT2 
                    AII3dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K00)  +  dROdP*(H(I,J)-Hpast(I,J))  + d2ROdP2*(P(I,JJ)-Ppast(I,J))*H(I,J) + (dROdP+dROdPp)/2*1.0*H(I,J)    +   (dROdP+dROdPp)/2*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K00)   )    *DT2 
                    AII4dt  =   -2* (RO(I,J)*(PAIAK*DX*SS*K10)  +                                                                                                        (dROdP+dROdPp)/2*(P(I,JJ)-Ppast(I,J))*(PAIAK*DX*SS*K10)   )    *DT2 
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
                                        
                    ! A(II+1) = dL(i,j)/dP(k=i-2,l=j)
                    A(1,II)=    kvot*DX3*AII1dxx    +   kvot*(AII1dx1 + AII1dx2 +AII1dx3)  +   AII1dt                                
                    ! A(II+2)= dL(i,j)/dP(k=i-1,l=j) 
                    A(2,II) =   kvot*DX3*AII2dxx    +   kvot*(AII2dx1 + AII2dx2 +AII2dx3)  +   AII2dt                   
                    ! A(II+3) = dL(i,j)/dP(k=i,l=j)
                    A(3,II) =   kvot*DX3*AII3dxx    +   kvot*(AII3dx1 + AII3dx2 +AII3dx3)  +   AII3dt 
                    ! A(II+4) = dL(i,j)/dP(k=i+1,l=j)
                    A(4,II) =   kvot*DX3*AII4dxx    +   kvot*(AII4dx1 + AII4dx2 +AII4dx3)  +   AII4dt     
                    
                    if ( isnan(A(1,i) )) then
                        WRITE(4,*)'A(1,i)=', A(1,i)
                    endif
                    if ( isnan(A(2,i) )) then
                        WRITE(4,*)'A(2,i)=', A(2,i)
                    endif
                    if ( isnan(A(3,i) )) then
                        WRITE(4,*)'A(3,i)=', A(3,i)
                    endif
                    if ( isnan(A(4,i) )) then
                        WRITE(4,*)'A(4,i)=', A(4,i)
                    endif
                    if ( isnan(A(5,i) )) then
                        WRITE(4,*)'A(5,i)=', A(5,i)
                    endif
                    
                       
50              CONTINUE
        RETURN
    END