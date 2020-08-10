! Subroutine to reed the imput data sheet where the parameters of the lubricant is defined. 
    ! The parameter Lubrication_param defines which set of pressure dependent newtonian equations should be used for the lubricant. 
    ! It also affects the parameters for the temperature simulations of the metals
    ! The parameter Shear_thin defines the non newtonian aspects of the lubracant based on the shear rates. 
    Subroutine Lubrication_def(HM0r1,  Z1, EDA01, alpha1)
    implicit none
    ! Input
    include     'inc_Grid.h'
    include     'inc_Itt.h'
    include     'inc_Outp.h'
    include     'inc_Ref.h'
    real    Z1 
    ! Calculations
    real        eda_1, eda_2, eda_3, eda_4, Delta_s
    real        YF, TG, PG
    real        C, D, taus
    ! Output
    real        EDA01, HM0r1, Alpha1
    include     'inc_Holmes.h'
    include     'inc_EDA_contact.h'
    include     'inc_Glass_const.h'
    include     'inc_NonNew.h'
	include     'inc_NonNew2.h'
    include     'inc_Rho.h'
    include     'inc_RLarsson.h'
    include     'inc_Shear_lim.h'
    include     'inc_Therm_param.h'
    include     'inc_Therm_cond.h'
    include     'inc_Temp_reduction.h'
    include     'inc_Visc.h'
    include     'inc_Yasutomi.h'
	include     'inc_Y_liu.h'  
    
    SAVE /RLarsson/                                                     ! Parameters for ref 10, R. Larssons formulation
    SAVE /Holmes/                                                       ! Viscosity and denisty param acc Holems et al
    SAVE /Rho/                                                          ! Density parameters
    SAVE /Visc/                                                         ! Lubrication parameters
    SAVE /NonNew/                                                       ! Non Newtonian viscosity parameters
    SAVE /NonNew2/
    SAVE /Yasutomi/                                                     ! Lubrication parameters acc Yasutomi
    SAVE /Y_Liu/
    SAVE /shear_lim/
    save /Therm_param/
    save /EDA_contact/
    save /Glass_const/
    save /Therm_cond/
    save /Temp_reduction/
    save /Itt/
    
    ! Read common variables for almost all ref
    read (7,*)
    read (7,*)
    read (7,*) Ta ,Z, EDA0, Pref, alpha, RA1, RA2           !D_IN: Input8.csv !D_OUT: alpha, EDA0, Pref, Z -> /Visc/; RA1, RA2 -> /Rho/; Ta -> /Yasutomi/
    alpha   = alpha*1E-8
    Pref    = Pref*1E8
    RA1     = RA1*1E-9             
    RA2     = RA2*1E-9
    
    ! Read values for material parameters defining the pressure dependent newtonian behaviour of the lubricant
    IF( Lub_param .EQ. 1 .or. Lub_param .EQ. 6) THEN !D_IN: /ref/ -> Lub_param ! Ends @ line: 534
        ! Newtonian acc X.Tans paper ref 30 31 and CH Venners paper from 1991  
        ! These parameters was already read above. 
        !Typical numbers are:
        !Z=0.68              ! Viscocity exponent                                [-]
        !EDA0=0.04           ! Ref viscocity                                     [kg/sm]
        !Pref=1.96E8         ! Ref pressure                                      [N/m^2]
        !alpha=2.2E-8        ! Pressure viscocity index.
        
        Ntime=NX/F                                          !D_IN: /Grid/ ->  NX; /ref/ -> F !D_OUT: Ntime -> /itt/ 

        ! Reformulated density parameters
        C=5.9E8
        D=1.34
        RA1=(D-1)/C
        RA2=1/C
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        
        
    ELSEIF( Lub_param .EQ. 2 ) THEN
        ! R Gohar 1988 Elastohydrodynamics p20-23, formulation of equations 
        ! The data was read above
        
        Z=alpha/(5.1*0.1*(log(EDA0)+9.67))     ! This eq I think should be updated according to lub_param 7          
        Ntime=NX/F
        
        RA1=0.6/10**(9)                         ! Gohar page 24. Scaled for GPa
        RA2=1.7/10**(9)
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        
    ELSE IF( Lub_param .EQ. 3)THEN
        ! Lub_param 3 is Input parameters ac Holmes et al Transient EHL point contact analysis 2003  
        
        ! Typical parameter values
        !EDA0=0.0096
        !alpha=17E-9
        !kH = 63.15E-6
        !xH = 5.1E-9
        !lH = 1.683E-9
        !gH = 2.266E-9
        
        read (7,*)
        read (7,*)
        read (7,*) kH, xH, lH, gH                           !D_in: input8.csv !D_OUT: kH, xH, lH, gH -> /Holmes/
        kH=kH*1E-6
        xH=xH*1E-9
        lH=lH*1E-9
        gH=gH*1E-9
        Z=alpha/(xH*LOG(EDA0/kH))
        
        Ntime=Um*NX/(F*Ua)                                  ! dT=F*um/um*dX             !D_IN: /outp/ -> Um, Ua 

        ! Reformulated density parameters
        RA1=gH-lH
        RA2=lH
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
    ELSE IF( Lub_param .EQ. 4)THEN
         ! R. Larsson et al. 2000 Lubricant properties for input to hydrodynamic and elastohydrodynamic lubricatio analyses formulation of equations 

        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original name of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New name of the inputs !Ta is the temperature EpsT0 is the density temperature coeff      !D_IN: Input8.csv !D_OUT: EpsT0, S0, RL_G0, Dz, Cz -> /RLarsson/
        
        !temperatures
        RL_T0= 40
        RL_Ta=Ta        ! To be on the safe side if used later ! !D_OUT: RL_T0, RL_Ta -> /RLarsson/
        Z=Dz+Cz*log10(1+Ta/135)
        S0=-S0
        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135)**S0)               !  EDA0 is updated in the newtonian subroutine,  EDA0=10**(-4.2+RL_G0*(1+Ta/135)**S0)           ! Eq (3)
        Ntime=NX/F
        
        ! Reformulated density parameter
        RL_c=1.5E-9                                             !D_OUT: RL_c -> /RLarsson/
        
        ! Viscosity coeff
        alpha=1.0/101325*(log(EDA0)+9.67)*(-1+(1+5.1/10**9*101325)**Z) ! 1 ATM = 101325 Pa
        
        read (7,*)
        read (7,*)
        read (7,*)  Temp_param, p_red, temp_red, temp_fac, tau_lim, tau_gamma, lim_factor     ! Different meaning of temp_red and temp_fac depending on temp param ==1 or == 2.   !D_OUT: Temp_param, p_red, temp_red, temp_fac, tau_lim, tau_gamma, lim_factor -> /Temp_reduction/
        read (7,*)
        read (7,*)
        read (7,*)  !This is for lub_param 5
        
        if ( temp_param==2) then
            temp_fac = temp_fac * PH                            ! Rescaling with Ph now instead of later
            !D_IN: /Visc/ -> PH
        endif
        temp_red = temp_red * PH
        tau_gamma = tau_gamma * PH                              ! Rescaling with Ph now instead of later
        tau_lim = tau_lim * 1e6
        p_red   = p_red * 1e9/PH
        
        ! These material parametres could be reformulated as input parameters. Now they are extracted from R. Larsson et al. 2000
        Cp_Oil  = 2080    !J/kgK=Nm/(kgK)                       ! Recalculated later
        CP_met  = 450
        K_met   = 47        !W/(mK)=Nm/(smK)    
        K_oil   = 0.14                                          ! Recalculated later
        rho_met = 7850
        rho_oil = 853   !PaoB
        !D_OUT: CP_oil, CP_met, K_met, K_oil, rho_met, rho_oil -> /Therm_param/
        
        ! These material parametres  could be reformulated as input parameters. Now they are extracted from R. Larsson et al. 2000
        Cp0     = 2080
        Cp1     = 0.41  *1e-9*PH
        Cp2     = 1.05  *1e-9*PH
        be0     = 6.5e-4
        be1     = 2.7   *1e-9*PH
        be2     = -1.5  *1e-9*PH*1e-9*PH
        ka0     = 0.154
        ka1     = 1.40  *1e-9*PH
        ka2     = 0.34  *1e-9*PH
        !D_OUT: Cp0, Cp1, Cp2, be0, be1, be2, ka0, ka1, ka2 -> /Therm_cond/
        
     ELSEIF( Lub_param .eq. 5) Then  
        ! NonNewtonian acc P.Eheret et all, Lub_param 5 -----------------------------------------------------------------------------
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        
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
        read (7,*) Tg0, YA1, YA2, YB1, YB2 !D_IN: Input8.csv !D_OUT: Tg0, YA1, YA2, YB1, YB2 -> /Yasutomi/
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
        read (7,*) YC1, YC2, Yedag, Yalfag !D_IN: Input8.csv !D_OUT: YC1, YC2, Yedag, Yalfag -> /Yasutomi/
        Yedag=Yedag*1E12
        Yalfag=1.0/(Yalfag**9)
        
        Tg=Tg0+YA1*LOG(1+YA2*0)
        pg=1/YA2*(EXP((1/YA1)*((Ta-Tg0)-1)))
        YF=1-YB1*LOG(1+YB2*0)
                
        EDA0=Yedag*10**(-(YC1*((Ta)-Tg)*YF) / (YC2+((Ta)-Tg)*YF)) 
        
        ! Density parameters (ref P.Ehret D.Dowsin and C.M. Taylor)
        RA1=0.6E-9
        RA2=1.7E-9
        
        ! From lub_param=4
        CP_oil  = 2000    !J/kgK=Nm/(kgK)
        CP_met  = 460
        K_met   = 47        !W/(mK)=Nm/(smK)
        K_oil   = 0.14  
        rho_met = 7850
        rho_oil = 853   !PaoB
        !D_OUT: CP_oil, CP_met, K_met, K_oil, rho_met, rho_oil -> /Therm_param/ ! !!! Duplicate
        
        ! Estimated for glass disk
        K_glass     = 0.84  ! From article
        rho_glass   = 4     ! (1.5-7.2)
        CP_glass    = 700   ! 670-840  
        !D_OUT: K_glass, CP_glass, rho_glass -> /Glass_const/
        
        ELSE IF( Lub_param .EQ. 7)THEN
         ! Wang 2017 formulation of equations 
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original namin of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New namin of the inputs !Ta is the temperature EpsT0 is the density temperature coeff
        
        ! Input in deg C, calculations in Kelvin
        RL_Ta=Ta+273
        RL_T0= 40+273
        
        Ta=RL_Ta
        ! Equations from Dong Zhu email conversation
        S0=0.042*(RL_T0-138)/(log(EDA0)+9.67)
        Z= alpha/(5.1e-9*(log(EDA0)+9.67))
        
        Ntime=NX/F
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5       
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
        CP_oil  = 2000      !J/kgK=Nm/(kgK)
        CP_met  = 460
        K_met   = 47        !W/(mK)=Nm/(smK)
        K_oil   = 0.14  
        rho_met = 7850
        rho_oil = 846       !Wang oil        ELSE IF( Lub_param .EQ. 7)THEN
             
    ELSE IF( Lub_param .EQ. 8)THEN
        ! Liu 2005 formulation of equations 
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original namin of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New namin of the inputs !Ta is the temperature EpsT0 is the density temperature coeff
        
        ! Input in deg C, calculations in Kelvin
        RL_Ta=Ta+273
        RL_T0= 30+273
        Ta=RL_Ta
        ! Equations from Dong Zhu email conversation
        S0=0.042*(RL_T0-138)/(log(EDA0)+9.67)
        Z= alpha/(5.1e-9*(log(EDA0)+9.67))
        
        Ntime=NX/F
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5       
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
        CP_oil  = 2000      !J/kgK=Nm/(kgK)
        CP_met  = 460
        K_met   = 47        !W/(mK)=Nm/(smK)
        K_oil   = 0.14  
        rho_met = 7850
        rho_oil = 870       !Wang oil
             
     ELSE IF( Lub_param .EQ. 9)THEN
        ! Bruyere 2010 formulation of equations 
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original namin of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New namin of the inputs !Ta is the temperature EpsT0 is the density temperature coeff called beta_DH
        
        ! Input in deg C, calculations in Kelvin
        RL_Ta=Ta+273
        RL_T0= 80+273
        Ta=RL_Ta
        ! Equations from Dong Zhu email conversation
        S0=0.0476*(RL_T0-138)/(log(EDA0)+9.67)    ! beta=0.0476
        Z= alpha*1.98*1E8/((log(EDA0)+9.67))     ! Almost equal to alpha/(5.1e-9*(log(EDA0)+9.67))
        
        RA1=RA1*1E9             
        RA2=RA2*1E9
        
        Ntime=NX/F
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5       
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
        CP_oil  = 2080      !J/kgK=Nm/(kgK)
        CP_met  = 450
        K_met   = 47        !W/(mK)=Nm/(smK)
        K_oil   = 0.15  
        rho_met = 7850
        rho_oil = 850       ! Bruyere oil
        
    ELSE IF( Lub_param .EQ. 10)THEN
        ! Hartinger 2008 formulation of equations 
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original namin of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New namin of the inputs !Ta is the temperature EpsT0 is the density temperature coeff called beta_DH
        
        ! Input in deg C, calculations in Kelvin
        RL_Ta=Ta+273
        RL_T0= 80+273
        Ta=RL_Ta
        ! Equations from Dong Zhu email conversation
        S0=0.0476*(RL_T0-138)/(log(EDA0)+9.67)    ! beta=0.0476
        Z= alpha/(5.1e-9*(log(EDA0)+9.67))
        
        RA1=RA1*1E9             
        RA2=RA2*1E9
        
        Ntime=NX/F
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5       
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
        CP_oil  = 2300      !J/kgK=Nm/(kgK)
        CP_met  = 450
        K_met   = 47        !W/(mK)=Nm/(smK)
        K_oil   = 0.15  
        rho_met = 7850
        rho_oil = 850       !
        
        Lub_param=9         ! Because same equations
        
    ELSE IF( Lub_param .EQ. 11 ) THEN
        ! Newtonian n=n_o^(alhpa*p^Z) Barue equation with an extra Z just to decrease pressure dependance if wanting to 

        Ntime=NX/F               
        ! Reformulated density parameters
        C=5.9E8  
        D=1.34
        RA1=(D-1)/C
        RA2=1/C
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 4
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
    ELSE IF( Lub_param .EQ. 12)THEN
        ! Kim 2001 formulation of equations 
        ! Case 1 - low slip

        read (7,*)
        read (7,*)
        read (7,*)  !This is for case 3
        read (7,*)
        read (7,*)
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original namin of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New namin of the inputs !Ta is the temperature EpsT0 is the density temperature coeff
        
        ! Input in deg C, calculations in Kelvin
        RL_Ta=Ta+273
        RL_T0= 40+273
        Ta=RL_Ta
        RA1=0.6E-9
        RA2=1.7E-9
        Ntime=NX/F
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5   
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
        CP_oil  = 2000      !J/kgK=Nm/(kgK)
        CP_met  = 460
        K_met   = 47        !W/(mK)=Nm/(smK)
        K_oil   = 0.14  
        rho_met = 7850
        rho_oil = 846  
        S0      = 0.042
        
    ELSE IF( Lub_param .EQ. 13)THEN
        ! Kim 2001 formulation of equations 
        ! Case 2 - High slip

        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 3
        read (7,*)
        read (7,*)
        !read (7,*) EDA0, alpha, kH, xH, lH, gH     Original namin of the inputs
        read (7,*) EpsT0,S0, RL_G0, Dz, Cz   !New namin of the inputs !Ta is the temperature EpsT0 is the density temperature coeff
        
        ! Input in deg C, calculations in Kelvin
        RL_Ta=Ta+273
        RL_T0= 40+273
        Ta=RL_Ta
        RA1=0.6E-9
        RA2=1.7E-9
        Ntime=NX/F
        
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5       
        read (7,*)
        read (7,*)
        read (7,*)  !This is for Lub_param 5
        
        CP_oil  = 2000      !J/kgK=Nm/(kgK)
        CP_met  = 460
        K_met   = 52        !W/(mK)=Nm/(smK)
        K_oil   = 0.124  
        rho_met = 7850
        rho_oil = 866   
        S0      = 0.0216
    Else
        WRITE(4,*)'Bad lubrication. lub_param ='
	    WRITE(4,*) lub_param
        stop 'Bad Lubrication'
    ENDIF

    ! Material defining the shear thinning behaviour. 
    IF( shear_thin .EQ. 0) Then !D_IN: /ref/ -> shear_thin
        ! No shear thinning applied
        tauS=0.0
        taua=0.0
        taua2=0.0
        xilim=0.0
        xi_param=0.0             !D_OUT: taua, taua2,  xilim,   xi_param -> /NonNew/
        
        read (7,*)
        read (7,*)
        read (7,*)
        read (7,*)
        read (7,*)
        read (7,*)
    ELSEIF( shear_thin .EQ. 1) Then
        ! NonNewtonian acc P.Eheret et al. 1997
        
        ! Non-Newtonian parameters
        !tauc_real=25E6           ! Limeting shear stress values [PA]
        !tauS=0.4                 ! The values of alpha*pherts 
        !taua2=0.0                ! The quadratic coefficiant alfa for glassy state
        !xilim=0.9                ! The maximum allowed value from xi
        read (7,*)
        read (7,*)
        read (7,*)  tauc_real, tauS, taua2, xilim, xi_param !D_IN: Input8.csv !D_OUT: tauc_real -> /NonNew2/
        tauc_real=tauc_real*1E6
        taua=tauS                ! The linear coefficiant alfa for glassy state of lubricant. At PH P=1. 
    
        tauc=ABS(tauc_real*B**2/(RX*EDA0*(Ua-Ub))) ! Rescaled to minimize calculations later   !D_OUT: Tauc -> /NonNew/
        if (Ua .Eq. Ub) tauc=1E18          !Big number to disable Non-Newtonian shear thinning is no shearing. 
        Ntime=NX/F
        Uslip=abs(Ua-Ub)                    !Slipspeed !D_OUT: Uslip -> /NonNew2/
        
        read (7,*)
        read (7,*)
        read (7,*) ! For shear thin = 2
    ELSE IF( shear_thin .EQ. 2)THEN
         ! Y Liu shear thinning 
        read (7,*)
        read (7,*)
        read (7,*) ! For shear thin = 1
        
        read (7,*)
        read (7,*)
        read (7,*)  L_n, L_G, L_h_limit, L_iter, L_stab !D_IN: Input8.csv !D_OUT: L_n, L_G, L_h_limit, L_iter, L_stab -> /Y_Liu/
        L_G=L_g*1e4
        L_h_limit=L_h_limit*1e-9
    ELSE IF( shear_thin .EQ. 3)THEN
         ! Y Liu shear thinning with basic temperature dependencd 
        read (7,*)
        read (7,*)
        read (7,*) ! For shear thin = §
        
        read (7,*)
        read (7,*)
        read (7,*)  L_n, L_G, L_h_limit, L_iter, L_stab
        L_G=L_g*1e4
        L_h_limit=L_h_limit*1e-9
    ELSE IF( shear_thin .EQ. 4)THEN
         ! Only temp increase for newtonian 4
        read (7,*)
        read (7,*)
        read (7,*) ! For shear thin = §
        
        read (7,*)
        read (7,*)
        read (7,*)  L_n, L_G, L_h_limit, L_iter, L_stab,shear_max,shear_min, temp_max !D_OUT: shear_max, shear_min, temp_max -> /shear_lim/
        L_G=L_g*1e4
        L_h_limit=L_h_limit*1e-9
        shear_min=shear_min*1e9
        shear_max=shear_max*1e9
    Else
        WRITE(4,*)'Bad shearthinning. Shear thin ='
	    WRITE(4,*) shear_thin
        stop 'Bad Shearthinning'
    ENDIF
        
    ! Alternative viscocity eq parameters
	A1=ALOG(EDA0)+9.67
	A2=5.1E-9*PH
	A3=0.59/(PH*1.E-9)
    
    ! To enable passing back to main. Could be formulated better
    Z1      = Z                       
    EDA01   = EDA0
    alpha1  = alpha
    HM0r1   = HM0r                     
    
        ! Define minimum viscosity if contact
        If( lub_param .EQ. 1)THEN    
            EDA_cont=EXP(alpha*Pref/Z*(-1+(1+1.0*PH/Pref)**Z)) !D_OUT: EDA_cont -> /EDA_contact/
        else If( lub_param .EQ. 11)THEN    
            EDA_cont=EXP(alpha*(1.0*PH)**Z)
        else If( lub_param .EQ. 2)THEN     
            ! Roelands equation acc Gohar
            EDA_1=5.1*PH*1.0/10**(9)
            EDA_2=(1.0+EDA_1)**Z
            EDA_3=(log(EDA0)+9.67)
            EDA_cont=EXP(EDA_3*(EDA_2-1.0))
        else If( lub_param .EQ. 3)THEN        
            EDA_cont=EXP(LOG(EDA0/kH)*((1+xH*PH*1.0)**Z-1)) 
        else If( lub_param .EQ. 4 .OR. lub_param .EQ. 60)THEN ! ??? Should lub_param here be 6? 
            
        
            eda_1=log(EDA0)+9.67
            eda_2=(1.0+5.1E-9*1.0*PH)**Z
            EDA_cont=EDA0*EXP(eda_1*(-1.0+eda_2))  !Roland Larsson pressure visc relation
            if( temp_param==2) then
                EDA_cont=tau_lim*1e-6 / 1           ! The limit shear strength at a highith of 1e-6m and a shear speed of 1 m/s
            endif
            
        elseif( Lub_param .eq.5) then
            EDA_cont=Yedag
            
        else If( lub_param .EQ. 6)THEN
            EDA_cont=EXP(log(EDA0+9.67)*((1+1.0*PH/Pref)**Z-1)) 
            
        else If( lub_param .EQ. 7 .or. Lub_param .EQ. 8)THEN
            eda_1=log(EDA0)+9.67
            eda_2=(1.0+5.1E-9*1.0*PH)**Z
            eda_3=((RL_Ta-138)/(RL_T0-138))**(-S0)
            EDA_cont=EDA0*exp(eda_1 * (-1.0 + eda_2 * eda_3))
            
        else if( lub_param == 9) then
            eda_1=log(EDA0)+9.67
            eda_2=(1.0+1.0*PH/(1.98*1E8))**Z
            eda_3=((RL_Ta-138)/(RL_T0-138))    
            delta_s=eda_1*eda_2*S0/(RL_T0-138)  
            eda_4=delta_s*(RL_Ta-RL_T0)
              
            EDA_cont=EDA0*exp(eda_1 * ((-1.0 + eda_2 )* eda_3)-eda_4)
        else
            WRITE(4,*)'STOP since NO defined EDA_Contact' 
	        WRITE(4,*) Ntime
            stop 'NO defined EDA_Contact'
        ENDIF
        
        ! Since the convergence routine stores if the time steps succseeded or not the number of time steps has an upper limit
        IF( ntime .GE. 1600) then 
            WRITE(4,*)'To large Ntime for the success vector'
	        WRITE(4,*) Ntime
            stop 'To large Ntime for the success vector'
        ENDIF
        
        
    return
    end
    