! Subrutine for itteration according to Reynolds equation of the timedependent solution
	SUBROUTINE ITER(t,H00,G0,DW,kvot,DT2)
    	implicit none
        COMMON /CTRA4/ D, A
        COMMON /ITERp/ID
        COMMON /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
        COMMON /NonNew/ tauc, taua, taua2,xilim                               ! Non Newtonian viscosity parameters
        COMMON /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag    ! Lubrication parameters acc Yasutomi
        COMMON /Rho/RA1,RA2                                             ! Density parameters
        COMMON /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
        COMMON /Itt/MK_stat, MK_time, ER_stat, ER_time, Ntime, KK   ! Itteration paramters
        COMMON /asp/asph,aspw                                           ! Apserity parameters
        COMMON /Method/term1,term2,term3,term4                          ! Controling the numerical method
        COMMON /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referenseCOMMON      /HREE/Emat, pg, Tg,YF
        COMMON  /COMAK2D/AK(0:300,0:300)                                       ! Parameters for HREE subrutine
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast  ! Past
        COMMON      /Residual/A5dxx,A5dx,A5dt                                   ! Parameters for the residual
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        real         A5dxx(1:300,1:300),A5dx(1:300,1:300),A5dt(1:300,1:300)
        real         Emat(1:300,1:300) ,pg(1:300,1:300), Tg(1:300,1:300), YF(1:300,1:300)
        real         ID(1:300), D(1:300),A(1:5,1:300)
        real        AK, RA1, RA2
        Integer NX, NY, NN, ref, t, tmeth
        integer Ntime, MK_stat, MK_time, KK, cyl 
        integer term1, term2, term3, term4
        integer MM,K, J,JJ,J1,J0, I, I0,I1, I2, IA, II
        Real        X0, XE, PAI
        Real        Z,  EDA0, Pref, alpha
        Real        Elast1, Elast2, EE
        Real        asph_real, aspw_real, ER_stat, ER_time
        Real        DX, width
        Real        asph, aspw,US, G0, DW
        real        ENDA, load, Sload
        real        A1, A2, A3, ALFA, G, AHM, HM0, UTL, HM0r, PH, H00
        real        tauc_real, tauc, tauS, taua, taua2,xilim
        real        Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
        real        kH, xH, gH, lH, L, M
        real        kvot, DT2
        real        dxidP, dxidPp, dROdP, dROdPp, dHdx, dHdxp, dPdx, dPdxp, d2ROdP2, d2xidP2 
        real        K00, K10, K20, K30, PAIAK
        real        D1, D2, D3,D4, D5,D1p, D2p, D3p,D4p, D5p
        real        P1,P2,P3,P4,P5,P1p,P2p,P3p,P4p,P5p
        real        KG1, C1, C2, DD
        SAVE      /Current/                                     ! Current timestep
        SAVE      /Residual/                                    ! Parameters for the residual
        SAVE      /ITERp/
        DATA    KG1,C1,C2/0,0.31,0.31/
    
        
        IF( t .GE. 0) C1=0.15
        IF( t .GE. 0) C2=0.15
        
        ! Initiating parameters
        KG1=1
    2   NN=(NY+1)/2
        MM=NX-1
        
        ! Itterate KK times untill checking if converged or not. 
        DO 100 K=1,KK
                    
            ! Start from NY=2 ang to to the senter node in Y direction and uppdate the pressure
            DO 70 J=2,NN
                
                ! Define nodenumbers for past and next nodes
                J0=J-1
                J1=J+1
                JJ=NY-J+1
                
                ! Check if its posible to save time by excluding nodes near the exit 
                IA=1
    8           MM=NX-IA
                
                IF(P(MM,J0).GT.1.E-6)GOTO 20
                IF(P(MM,J).GT.1.E-6)GOTO 20
                IF(P(MM,J1).GT.1.E-6)GOTO 20
                IA=IA+1
                IF(IA.LT.NX)GOTO 8
                GOTO 70
    20          IF(MM.LT.NX-1)MM=MM+1 
        
                D2=0.5*(EPSx(1,J)+EPSx(2,J))
                D2p=0.5*(EPSxpast(1,J)+EPSxpast(2,J))
                
                ! Uppdate the presure for NX=2 untill there is zero pressure near the exit. at NX=1 BCs restric this pressure to be 0
                CALL Res(D2, D2p, MM, J0,J1, JJ, J,C1,C2,kvot,DT2)
                       
                ! Subrutien for adjusting the pressureupdates according to nearby pressureupdates in X-direction
                CALL TRA4(MM)
                    
                DO 60 I=2,MM
                    IF(ID(I).EQ.2)GOTO 60
                    IF(ID(I).EQ.0)GOTO 52
                    DD=D(I+1)
                    IF(I.EQ.MM)DD=0
                    P(I,J)=P(I,J)+C2*(D(I)-0.25*(D(I-1)+DD))
                    IF(J0.NE.1)P(I,J0)=P(I,J0)-0.25*C2*D(I)
                    IF(P(I,J0).LT.0.)P(I,J0)=0.0
                    IF(J1.GE.NN)GOTO 54
                    P(I,J1)=P(I,J1)-0.25*C2*D(I)
                    GOTO 54
    52              P(I,J)=P(I,J)+C1*D(I)
    54              IF(P(I,J).LT.0.0)P(I,J)=0.0
    60          CONTINUE
70          CONTINUE
            
            NN=(NY+1)/2
                
            ! Mirror the pressure around the X-axis
            DO 80 J=2,NN
                JJ=NY+1-J
                DO 80 I=1,NX
80          P(I,JJ)=P(I,J)
            
            ! Update the borderpressures
            DO I=2,NX
                P(I,1)=P(I,2)
                P(I,NY)=P(I,NY-1)
            ENDDO
            
            ! Enforcing boundaryconditions
            DO J=1,NY
                P(1,J)=0.0
                P(NX,J)=0.0
            ENDDO
            
            ! Now when the pressure is updated, update the other parameters based on the new pressure
            CALL HREE(H00,G0,t)
100     CONTINUE
        
        load=sum(P)                     ! The load
        Sload=DX*DX*load/G0             ! Normalized load
        DW=Sload-1.0                    ! 1-scaled loadsum
             
        RETURN
    END