! Initiate output files and pring out some values of used parameters. These values could bee good to have to double check that the input has been read as intended. 
    SUBROUTINE Printout(H00,Tauc_real,asph_real, aspw_real,HM0,alpha_bar) !D_IN: H00, Tauc_real, asph_real, aspw_real, HM0, alpha_bar as arguments
    implicit none
    include     'inc_Asp.h'
    include     'inc_Dimless.h'
    include     'inc_DZ_com.h'
    include     'inc_Grid.h'
    include     'inc_Geom5.h'
    include     'inc_G0DT.h'
    include     'inc_Holmes.h'
    include     'inc_Higginson.h'
    include     'inc_Hoglund_ball.h'
    include     'inc_Itt.h'
    include     'inc_Method.h'
    include     'inc_NonNew.h'
    include     'inc_Outp.h'
    include     'inc_Rho.h'
    include     'inc_Ref.h'
    include     'inc_RLarsson.h'
    include     'inc_Visc.h'
   	include     'inc_Yasutomi.h'
	include     'inc_Y_liu.h'
    include     'inc_Tsize.h'
    include     'inc_Temp_param.h'
	include     'inc_Temp_reduction.h'
    include     'inc_Temp_conduction.h'
    include     'inc_Output_control.h'
    include     'inc_Single_step.h'

        Integer  i, Ntime_s
        Real    asph_real, aspw_real
        real    H00, U, SRR
        real    Tauc_real, HM0, alpha_bar
                
        ! Initaial data output -----------------------------------------------------------------------------------------------------------------
	    OPEN(8,FILE='FILM.DAT',STATUS='UNKNOWN') ! ??? Why make the files here if they aren't used in this subroutine
	    OPEN(10,FILE='PRESSURE.DAT',STATUS='UNKNOWN')
        !OPEN(2,FILE='EDA.DAT',STATUS='UNKNOWN')
        !OPEN(3,FILE='EPS.DAT',STATUS='UNKNOWN')
        !OPEN(5,FILE='Pside.DAT',STATUS='UNKNOWN')
        !OPEN(11,FILE='A5dxx.DAT',STATUS='UNKNOWN')
        !OPEN(12,FILE='A5dx.DAT',STATUS='UNKNOWN')
        !OPEN(13,FILE='A5dt.DAT',STATUS='UNKNOWN')
        OPEN(23,FILE='Temp.DAT',STATUS='UNKNOWN')

        OPEN(50,FILE='Temperature.DAT',STATUS='UNKNOWN')
        OPEN(51,FILE='MetalA_temp.DAT',STATUS='UNKNOWN')
        
        OPEN(52,FILE='Hoglund_param.DAT',STATUS='UNKNOWN')
        OPEN(53,FILE='DZ_vector.DAT',STATUS='UNKNOWN')
        
        Open(70, File = 'ExecutionTime.DAT',Status='Unknown')
        
        if (collect_full_data) then
            OPEN(60, FILE = 'Film_initial.DAT', STATUS = 'UNKNOWN')
            OPEN(54,FILE='MetalB_temp.DAT',STATUS='UNKNOWN')
            OPEN(55, FILE = 'EDAx.DAT', STATUS = 'UNKNOWN')
            OPEN(56, FILE = 'EDAy.DAT',STATUS = 'UNKNOWN')
            OPEN(57, FILE = 'EPSx.DAT', STATUS = 'UNKNOWN')
            OPEN(58, FILE = 'EPSy.DAT', STATUS = 'UNKNOWN')
            OPEN(59, FILE = 'RO.DAT', STATUS = 'UNKNOWN')
            OPEN(61, FILE = 'xi.DAT', STATUS = 'UNKNOWN')
            OPEN(62, FILE = 'w.DAT', STATUS = 'UNKNOWN')
            OPEN(63, FILE = 'drodP_mat.DAT', STATUS = 'UNKNOWN')
            !OPEN(71, File = 'Debug.DAT', Status = 'Unknown')
        else if (single_step_only) then
            OPEN(64, File = 'Pressure_single_step.DAT',Status = 'Unknown')
            OPEN(65, File = 'Film_single_step.DAT',Status = 'Unknown')
            OPEN(66, File = 'Temp_single_step.DAT',Status = 'Unknown')
            OPEN(67, File = 'MetalA_Temp_single_step.DAT',Status = 'Unknown')
            !Open(70, File = 'Debug_single_step.DAT',Status='Unknown')
            !OPEN(68, FILE = 'w_input_single_step.DAT', Status = 'Unknown')
            !OPEN(69, FILE = 'H_input_single_step.DAT', Status = 'Unknown')
        endif
        
        
        ! To Output file
        WRITE(4,*)'Simulation of EHL'
        WRITE(4,*)''
        WRITE(4,*)'The reference choosen:' 
        WRITE(4,*)'Meth,            Tmeth,      Geom,       PAIAK,      Asp shape,      contact_alg    and  Multi_grid_param'
        WRITE(4,*) Meth,            Tmeth,      Geom,       PAIAK,      asp_shape ,     contact_alg,        Multi_grid_param!D_IN: /Ref/
        WRITE(4,*)''
        
        WRITE(4,*)'Grid infromation'
        WRITE(4,*)'NX,              NY,         X0,         XE,     DZ_method'
	    WRITE(4,*) NX,              NY,         X0,         XE,     DZ_method                       !D_IN: /Grid/ -> NX, NY, X0, XE; /Dz_com/ -> DZ_method 
        WRITE(4,*)''
        WRITE(4,*)'DX,              DT,         kk_stat,     kk_time,    Ntime'
        WRITE(4,*) DX,              DT,         kk_stat,     kk_time,    Ntime                      !D_IN: /Itt/ -> Ntime, kk_time, kk_stat; /G0DT/ -> DT; /Grid/ -> Dx
        WRITE(4,*)''
        
        WRITE(4,*)'Geometry information'
        WRITE(4,*)'Rx,              W0'
        WRITE(4,*) Rx,              W0                                                              !D_IN: /outp/ -> Rx, W0
        WRITE(4,*)' Ua,             Um,         U2_over_Us'
        WRITE(4,*)  Ua,             Um,         U2_over_Us                                          !D_IN: /outp/ -> Ua, Um, U2_over_Us
        WRITE(4,*)'Lub param,       Shear thin  lub_temp'
        WRITE(4,*) Lub_param,       Shear_thin, lub_temp                                            !D_IN: /Ref/ -> Lub_param, shear_thin; /temp_param/ -> lub_temp
        Write(4,*)''
        WRITE(4,*)'B,               PH,         By,         Ry,             PH_new (For geom5)'
	    WRITE(4,*) B,               PH,         By,         Ry,             PH_new                  !D_IN: /Visc/ -> PH; /outp/ -> Ry, By, B; /Geom5/ -> PH_new
        WRITE(4,*)'EE'
        WRITE(4,*) EE
        WRITE(4,*)''
        
        WRITE(4,*)'Asperity'
        WRITE(4,*)'asph_real,       aspw_real,  asph,       aspw'
        WRITE(4,*) asph_real,       aspw_real,  asph,       aspw                                    !D_IN: /asp/ -> asph, aspw
        Write(4,*)'asp_shape,       asph2,      aspw2,      asp_ratio,      Hminimum '
        WRITE(4,*) asp_shape,       asph2,      aspw2,      asp_ratio,      Hminimum                !D_IN: /asp/ -> asph2, aspw2, asp_ratio, Hminimum; /Ref/ -> asp_shape
        WRITE(4,*)'Tsize'
        WRITE(4,*) Tsize
        WRITE(4,*)''
        
        WRITE(4,*)'Iteration'
        Write(4,*)'Ntime,           MK_stat,    MK_time,    ER_stat,        ER_time'
        Write(4,*) Ntime,           MK_stat,    MK_time,    ER_stat,        ER_time                 !D_IN: /Itt/ -> Ntime, MK_stat, MK_time, ER_stat, ER_time
        WRITE(4,*)'H00,             HM0,        HM0r'
	    WRITE(4,*) H00,             HM0,        HM0r
        WRITE(4,*)''
        
        WRITE(4,*)'Temp Ta='
        WRITE(4,*) Ta                                                                               !D_IN: /Yasutomi/ -> Ta
        
        WRITE(4,*)'Density'
        WRITE(4,*)'RA1,         RA2'
        WRITE(4,*) RA1,         RA2                                                                 !D_IN: /Rho/ -> Ra1, Ra2
        
        WRITE(4,*)'Viscosity'
        WRITE(4,*) 'Z,          EDA0,       Pref,       alpha,      ENDA'
        if(lub_param ==4 .or. lub_param .GE. 7) then! EDA0 is treated differently in these cases since affects the power generation
            WRITE(4,*) Z,           EDA0,       Pref,       alpha,      ENDA*EDA0                   !D_IN: /Visc/ -> EDA0, ENDA, alpha, Pref, Z            
        else
            WRITE(4,*) Z,           EDA0,       Pref,       alpha,      ENDA
        endif
        WRITE(4,*)'Heat partitioning '
        WRITE(4,*)'Shear_band,  Shear_xi'
        WRITE(4,*) Shear_band,  Shear_xi
        
        WRITE(4,*)''
        
        If( Lub_param .EQ. 3) Then
            WRITE(4,*)'Holmes 2003 parameters'
            WRITE(4,*)'lH,      gH,         xH,         kH'
            WRITE(4,*) lH,      gH,         xH,         kH                                          !D_IN: /Holmes/ -> kH, xH, gH, lH
            WRITE(4,*)''
            
        ELSEIF( Lub_param .eq. 4) THEN
            WRITE(4,*)'R.Larssons 2000 parameters'
            WRITE(4,*)'EpsT0,   S0,         RL_G0,      Dz' 
            WRITE(4,*) EpsT0,   S0,         RL_G0,      Dz                                          !D_IN: /RLarsson/ -> EpsT0, S0, RL_G0, Dz
            WRITE(4,*)'Cz,      RL_T0,      RL_c'
            WRITE(4,*) Cz,      RL_T0,      RL_c                                                    !D_IN: /RLarsson/ -> Cz, RL_T0, RL_c
            WRITE(4,*)''
            
        ELSEIF( Lub_param .eq. 5) THEN
            WRITE(4,*)'Non Nwetonian Ehret Yasutomi equation'
            WRITE(4,*)'Ta,      Tg0,        YA1,        YA2,        YB1,        YB2'
            WRITE(4,*) Ta,      Tg0,        YA1,        YA2,        YB1,        YB2                 !D_IN: /Yasutomi/ -> Ta, Tg0, YA1, YB1, YB2
            WRITE(4,*)'YC1,     YC2,        Yedag,      Yalfag'
            WRITE(4,*)YC1,      YC2,        Yedag,      Yalfag                                      !D_IN: /Yasutomi/ -> YC1, YC2, Yedag, Yalfag
            WRITE(4,*)''
        ENDIF
        
        
        IF(shear_thin .EQ. 1)THEN
            WRITE(4,*)'Shearthinning and Solidification'
            WRITE(4,*)'Tauc_real, Taua,     Taua2,      xilim,      xi_param '                      ! ??? Is Tauc_real assigned anywhere?
            WRITE(4,*) Tauc_real, Taua,     Taua2,      xilim,      xi_param                        !D_IN: /NonNew/ -> Taua, Taua2, xilim, xi_param 
            WRITE(4,*)''
        ELSEIF(shear_thin .EQ. 2)THEN
            WRITE(4,*)'Shear thinning acc Y. Liu'
            WRITE(4,*)'L_n,     L_G,	    L_h_limit,	L_iter,	    L_stab '
            WRITE(4,*) L_n,	    L_G,	    L_h_limit,	L_iter,	    L_stab                          !D_IN: /Y_Liu/ -> L_N, L_G, L_H_limit, L_iter, L_stab
            WRITE(4,*)''
        ELSEIF(shear_thin .EQ. 3)THEN
            WRITE(4,*)'Shear thinning acc Y. Liu with basic temperature dependence'
            WRITE(4,*)'L_n,     L_G,	    L_h_limit,	L_iter,	    L_stab '
            WRITE(4,*) L_n,	    L_G,	    L_h_limit,	L_iter,	    L_stab 
            WRITE(4,*)''
        ENDIF
        
        IF(temp_param .GE. 1)THEN
            WRITE(4,*)'Shear limit included in simulations'
            WRITE(4,*)'Temp_param, p_red,   temp_red, temp_fac, tau_lim, tau_gamma, lim_factor'
            WRITE(4,*) Temp_param, p_red,   temp_red, temp_fac, tau_lim, tau_gamma, lim_factor      !D_IN: /Temp_redutcion/ -> Temp_param, p_red,   temp_red, temp_fac, tau_lim, tau_gamma, lim_factor
            WRITE(4,*)''
        ENDIF
        
        IF( geom == 6) then
            WRITE(4,*)'Höglund ball simulation'
            WRITE(4,*)'V_ball, Vv_ball,Vh_ball,theta_ball, M_ball, I_ball, omega_ball, DT0'
            WRITE(4,*) V_ball, Vv_ball,Vh_ball,theta_ball, M_ball, I_ball, omega_ball, DT0          !D_IN: /Hoglund_ball/ -> V_ball, Vv_ball,Vh_ball,theta_ball, M_ball, I_ball, omega_ball, DT0
            WRITE(4,*) ''
        endif
        
        if( asp_shape .ge. 200 .and. asp_shape .le. 209) then
            Ntime_s = Ntime + (135-15)*Ntime/Nx          ! Ntime/Nx=1/F  ! added time steps because some the shape is slightly bigger due to the added frame
            WRITE(4,*)'Rough surface simulations'
            WRITE(4,*)'The updated number of timesteps are, Ntime=', Ntime_s
            WRITE(4,*) ''
        ENDIF
        

    
        WRITE(4,*)'Dimensionless params'
        WRITE(4,*)'L ,          M and       alpha_bar'
        WRITE(4,*) L,           M,          alpha_bar                                               !D_IN: /Dimless/ -> L,M
        WRITE(4,*)'U,           G,          W '
        WRITE(4,*) Uhd,         Ghd,        Whd                                                     !D_IN: /Higginson/ -> Uhd, Ghd, Whd 
        WRITE(4,*)''
        
        WRITE(4,*)'Choise of method'
        WRITE(4,*)'The four numbers with the witght of the nodes'
        WRITE(4,*) term1,       term2,      term3,      term4                                       !D_IN: /Method/ -> term1, term2, term3, term4        

        WRITE(*,*)'               Wait please'
        WRITE(4,*)''
        WRITE(4,*)'               Wait please'
        WRITE(4,*)''
        
        !IF(DZ_method .NE. 0) then
            WRITE(53,*)(DZ_vec(I),I=1,40)                                                           !D_IN: /Dz_com/ -> Dz_vec
        !ENDIF
        Write(70,*)'Starting the simulation: '
            
            
        RETURN
    END
!