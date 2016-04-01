	PROGRAM POINTEHL
    implicit none                                                   ! För att kräva initiering av samtliga variabler http://www.obliquity.com/computer/fortran/common.html
    Integer NX, NY
    integer Ntime, MK_stat, MK_time, KK
    integer term1, term2, term3, term4, tmeth
    integer meth, ref, nr_N
    Real X0, XE, PAI, PAIAK, Y0
    Real RX, W0, Ua, Ub
    Real Z,  EDA0, Pref, alpha
    Real Elast1, Elast2, EE
    Real asph_real, aspw_real, ER_stat, ER_time
    Real DX, DT, width
    Real asph, aspw,US, G0
    real ENDA, HM0f, Uhd,Whd,Ghd
    real B, PH, H00, U, SRR, H00past, DWpast, Sloadpast 
    real A1, A2, A3, W, ALFA, G, AHM, HM0, UTL, HM0r
    real tauc_real, tauc, tauS, taua, taua2,xilim
    real Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
    real RA1,RA2, C, D
    real kH, xH, gH, lH, L, M, F
    COMMON /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
    COMMON /NonNew/ tauc, taua, taua2, xilim                        ! Non Newtonian viscosity parameters
    COMMON /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag    ! Lubrication parameters acc Yasutomi
    COMMON /Rho/RA1,RA2                                             ! Density parameters
    COMMON /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
    COMMON /Itt/MK_stat, MK_time, ER_stat, ER_time, Ntime, KK       ! Itteration paramters
    COMMON /asp/asph,aspw                                           ! Apserity parameters
    COMMON /outp/ W0,EE,RX,US, Ua, B                                ! Output parameters
    COMMON /Holmes/kH, xH, gH, lH                                   ! Viscosity and denisty param acc Holems et al
    COMMON /Dimless/L,M                                             ! Diminsionless params
    COMMON /Method/term1,term2,term3,term4                          ! Controling the numerical method
    COMMON /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referense
    COMMON /Higginson/Uhd,Ghd,Whd                                   ! Dimensionless params acc Higginsson
    COMMON      /H00/ H00past, DWpast, Sloadpast                    ! Params for uppdating H00
    COMMON /G0DT/G0,DT                                              ! G0 and DT
    
    ! Input data -----------------------------------------------------------------------------------------------------------
    open (unit = 7, file = "Input3.csv")
    read (7,*)
    read (7,*) meth, ref, tmeth
    
    ! Model geometry ------------------------------------------------------------------------------
    !NX=150
    !NY=70                          ! NY <= NX 
    !X0=-2.0
    !XE=1.5
    PAI=3.14159265                  ! Only first 7 digits that counts since single pressition
    !F=2.0 !0.0625                  ! Kvot between time and space steps
    !Ntime=450                      ! Number of steps in time direction. The code is based on that the asperity travels from X0 to XE
    read (7,*)
    read (7,*)
    read (7,*) NX, NY, X0, XE, F
    
    DX=(XE-X0)/(NX-1.)              ! Space increment
    Y0=-0.5*DX*NY+0.5*DX            ! Width in Y-direction
    width=DX*NY                     ! Width of model
    
    ! Choice of method? ---------------------------------------------------------------------------
    IF( meth .EQ. 1) THEN
        ! 1st dx
        term1=0
        term2=-2
        term3=2
        term4=0
    ELSE IF( meth .EQ. 2)THEN
        ! Backspace
        term1=1
        term2=-4
        term3=3
        term4=0
    ELSE            
        ! Central space
        term1=0
        term2=-1
        term3=0
        term4=1
    ENDIF
    

    ! EHL parameters -------------------------------------------------------------------------------------------------------------
    !RX=0.0136/2 ! 0.0136/2       ! Radius  acc Dave Hannes                                          [m]
    !W0=10.0E5  ![N]       ! Load parameter  to get ish Pherts=2GPa                                  [N/m]
    !Ua=23.0        ! =U2 The speed of the surface with the asperity    [m/s]
    !Ub=17.0        ! Velocity of lower surface                         [m/s]       ! Mean velocity is US=(Ua+Ub)/2
    read (7,*)
    read (7,*)
    read (7,*) RX, W0, Ua, Ub
    US=(Ua+Ub)/2                    ! Mean velocity of fluid                            [m/s]
    
    ! Asperity Data
    !asph_real= 1.5E-6   ! Real asperity height                             [m]
    !aspw_real= 50E-6     ! Real asperity radius                            [m]
    read (7,*)
    read (7,*)
    read (7,*)  asph_real, aspw_real
    asph_real= asph_real*1E-6
    aspw_real= aspw_real*1E-6
    
     !Elast1=206E9    ! Elastic modulus of first body                     [Pa]
    !Elast2=206E9    ! Elastic modulus of second body                    [Pa]
    !EE=2.26E11      ! Equvivalent E =E'                                 [N/m^2]
    read (7,*)
    read (7,*)
    read (7,*) Elast1, Elast2, EE
    Elast1=Elast1*1E9
    Elast2=Elast2*1E9
        
    ! Choice of input structure ------------------------------------------------------------------------------------------------ 
    IF( ref .LT. 3 .OR. ref .EQ. 6) THEN
        ! Newtonian acc X.Tans paper ref 30 31 or the asperity runs, Cylinder-----------------------------------------------------------------------
        EE=EE*1E9                           ! Equvivalent Elastic modulus [Pa]
        W0=W0*1E5
        nr_N = 2*NX*NY/(XE-X0)              ! Number of nodes in the contact
        PAIAK=0.15915489024
        
        !Lubrication
        !Z=0.68              ! Viscocity exponent                                [-]
        !EDA0=0.04           ! Ref viscocity                                     [kg/sm]
        !Pref=1.96E8         ! Ref pressure                                      [N/m^2]
        !alpha=2.2E-8        ! Pressure viscocity index.
        read (7,*)
        read (7,*)
        read (7,*) Z, EDA0, Pref, alpha
        alpha=alpha*1E-8
        Pref=Pref*1E8
        
        IF(ref .EQ. 2) EE=2/((1-0.3**2)/Elast1+(1-0.3**2)/Elast2)               ! If asperity case, ref2, Equvivalent Elastic modulus [Pa]
        
        Ntime=NX
        G0=width*PAI/2                  ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
        B=SQRT(8*W0*RX/(PAI*EE))        ! Contact halfwidth                     [m]
	    PH=2*W0/(PAI*B)                 ! Hertzian pressure                     [N/m^2]  
        
        IF(ref .EQ. 1) asph_real=asph_real*B**2/RX*1E6     ! IF acc X.Tan paper Rescaling the asperity since X.Tan gives asperity in relativ numbers.
        
        
        ! Reformulated density parameters
        C=5.9E8
        D=1.34
        RA1=(D-1)/C
        RA2=1/C
        
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 3
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 4
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 4        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 4
        
    ELSE IF( ref .EQ. 3)THEN
        ! Input parameters ac Holmes et al Transient EHL point contact analysis 2003, Ball ---------------------------------------------------------------------------------
        EE=EE*1E9                                                           ! Equvivalent Elastic modulus [Pa]
        
        PAIAK=0.2026423
        nr_N = 2*NX*NY/(XE-X0)*1/(-Y0)              ! Number of nodes in the contact
        
        !EDA0=0.0096
        !alpha=17E-9
        !kH = 63.15E-6
        !xH = 5.1E-9
        !lH = 1.683E-9
        !gH = 2.266E-9
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 1
        read (7,*)
        read (7,*)
        read (7,*) EDA0, alpha, kH, xH, lH, gH
        alpha=alpha*1E-9
        kH=kH*1E-6
        xH=xH*1E-9
        lH=lH*1E-9
        gH=gH*1E-9
        Z=alpha/(xH*LOG(EDA0/kH))
        
        Ntime=US*NX/(F*Ua)
        G0=2.0/3.0*PAI                                              ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan 
        PH=(3*W0*EE**2/(2*PAI**3*RX**2))**0.3333333              ! Hertzian pressure                     [N/m^2]
        B=PAI*PH*RX/EE                                              ! Contact halfwidth                     [m]
	    
        ! Reformulated density parameters
        RA1=gH-lH
        RA2=lH
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 4
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 4        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 4
    
    ELSE  
        ! NonNewtonian acc P.Eheret et all, ref 4 which is a Ball, or cylinder with asperity, ref 5------------------------------------------------------------------------------

        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 1
        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 3
        
        ! Lubrication
        ! Yasutomi viscosity parameters for 5P4E !HVI650 (ref P.Ehret D.Dowsin and C.M. Taylor)
        !Ta=92.0           ! Current avrage temperature deg Celcius
        !Tg0= -28.2 !-30.4      ! Reference temperature, deg Celcius
        !YA1= 134.2 ! 309.1       ! Deg Celcius
        !YA2= 1.929E-9 ! 0.3064E-9  ! 1/Pa
        !YB1= 4.815 ! 0.2186      ! [-]
        !YB2=0.160E-9 !29.99E-9    ! 1/Pa
        read (7,*)
        read (7,*)
        read (7,*) Ta, Tg0, YA1, YA2, YB1, YB2
        YA2=YA2*1E-9
        YB2=YB2*1E-9
    
        !YC1= 16.01 ! 10.264     ! [-]
        !YC2=20.69 ! 27.04       ! Deg Celcius
        !Yedag=1E12       ! Pa*s
        !Yalfag= 3.814697E-15!=40**(-9)    ! 1/Pa
        !EDA0=1.0!0.04 ! Resetting EDA0 since not scaled with this model.  
        !Z=0.68           ! Viscocity exponent for HMO                                [-]
        read (7,*)
        read (7,*)
        read (7,*) YC1, YC2, Yedag, Yalfag, EDA0, Z
        Yedag=Yedag*1E12
        Yalfag=Yalfag*1E-15
    
        ! Density parameters (ref P.Ehret D.Dowsin and C.M. Taylor)
        RA1=0.6E-9
        RA2=1.7E-9
    
        ! Non-Newtonian parameters
        !tauc_real=25E6           ! Limeting shear stress values [PA]
        !tauS=0.4                 ! The values of alpha*pherts 
        !taua2=0.0                ! The quadratic coefficiant alfa for glassy state
        !xilim=0.9                ! The maximum allowed value from xi
        read (7,*)
        read (7,*)
        read (7,*)  tauc_real, tauS, taua2, xilim
        tauc_real=tauc_real*1E6
        taua=tauS                ! The linear coefficiant alfa for glassy state of lubricant. At PH P=1. 
    
        Ntime=NX
        IF( ref .EQ. 4) THEN
            EE=EE*1E9                                                           ! Equvivalent Elastic modulus [Pa]
            PAIAK=0.2026423
            G0=2.0/3.0*PAI                                              ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan 
                    
            PH=(3*W0*EE**2/(2*PAI**3*RX**2))**0.3333333              ! Hertzian pressure                     [N/m^2]  
            B=PAI*PH*RX/EE                                              ! Contact halfwidth                     [m]
            nr_N = 2*NX*NY/(XE-X0)*1/(-Y0)              ! Number of nodes in the contact
        ENDIF
        
        
        
        IF( ref .EQ. 5) THEN ! Nonnewtonian fluid with Cylinder for asperity runs
            EE=2/((1-0.3**2)/Elast1+(1-0.3**2)/Elast2)  ! Equvivalent Elastic modulus [Pa]
            PAIAK=0.15915489024
            G0=width*PAI/2                              ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
            
            W0=W0*1E5
            B=SQRT(8*W0*RX/(PAI*EE))                    ! Contact halfwidth                     [m]
	        PH=2*W0/(PAI*B)                             ! Hertzian pressure                     [N/m^2]
            nr_N = 2*NX*NY/(XE-X0)              ! Number of nodes in the contact
    

        ENDIF

        
    ENDIF
    
    
    ! Alternative viscocity eq parameters
	A1=ALOG(EDA0)+9.67
	A2=5.1E-9*PH
	A3=0.59/(PH*1.E-9)
    
    
    
    
    ! Accurace paramters -----------------------------------------------------------------------------------------------------------
    !MK_stat = 200   ! Maximum numbers of itterations for static solution
    !MK_time = 200   ! Maximum numbers of itterations for timedep solution
    !ER_stat = 1E-4  ! Maximum numbers of itterations for static solution
    !ER_time = 2E-5  ! Maximum numbers of itterations for timedep solution
    !KK=15               ! Number of internal itterations
    
    read (7,*)
    read (7,*)
    read (7,*)     MK_stat, MK_time, ER_stat, ER_time, KK
    ER_stat=ER_stat*1E-8*nr_N               ! The error limit is set proportional to the number of nodes in the contact area
    ER_time=ER_time*1E-8*nr_N
    
    
    !H00=-0.1                   ! Initial offset for lubrication height     
    read (7,*)
    read (7,*)
    read (7,*)  H00, HM0f
    
    
    ! End Input data -------------------------------------------------------------------------------------------------------------
    
    ! Pre calculations to generate dimentionless paramters ------------------------------------------------------------------------
    
    ! Geometry setup
    DT=(US*(XE-X0))/(Ntime*Ua)      ! The time increment. T=t*Us/a. t_end=(XE-X0)*a/Ua => T_end=(XE-X0)Us/Ua
    SRR=(Ua-Ub)/US                  ! Slide to Roll Ratio

             
	U=EDA0*US/(2.*EE*RX)        ! Alternative way of defining U. ENDA ends up correct anyway.  Probably for HM0 that U is define this strange
    
	W=2.*PAI*PH/(3.*EE)*(B/RX)**2
	ALFA=Z*5.1E-9*A1
	G=ALFA*EE
	AHM=1.0-EXP(-0.68*1.03)
	HM0=1.0!    3.63*(RX/B)**2*G**0.49*U**0.68*W**(-0.073)*AHM      ! Parameter for uppdating H0, the fim thicknes ofset value
    HM0r=HM0*HM0f                                            ! Reduced lubrication uppdation. 
    ENDA=12 *EDA0 *Us * RX**2 / ( B**3*PH)                  ! Epsilon in Reynolds equation. = 12.*(2*U)*(EE/PH)*(RX/B)**3 
	UTL=EDA0*US*RX/(B*B*2.E7)
    H00past=H00
    
    ! Moes and Bosma non-dimentional paramters
    L = alpha * EE * (2*EDA0*Us/(EE*RX))**0.25
    M = (W0/(EE*RX**2))*(EE*RX/(2*EDA0*Us))**0.75
    
    ! Higginson and Dowsson non-dim param from Spikes review from 2006 for line load
    Uhd=Us*EDA0/(EE*RX)
    Ghd=alfa*EE             
    Whd=W/(EE*RX)

	
    ! Asperity  rescaling data
    asph=asph_real*RX/B**2      ! Normalized height of asperity            [-]   
    aspw=aspw_real/B            ! Normalized radius of asperity            [-]
    
    ! Non-Newtonian rescaling data
    tauc=ABS(tauc_real*B**2/(RX*EDA0*(Ua-Ub))) ! Rescaled to minimize calculations later
    if (Ua .Eq. Ub) tauc=1E18          !Big number to disable Non-Newtonian shear thinning. 

    ! Printout the data
    CALL Printout(H00,Tauc_real,asph_real, aspw_real,HM0,DT)

    ! START simulation ---------------------------------------------------------------------------------------------------------------------------
    ! Only passing dimentionless paramters and setup parameters
    CALL Time(H00)

    

	STOP
    END   
