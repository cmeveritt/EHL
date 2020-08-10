! A subroutine to calculate the temperature insde the fluid film at the node i,j
! For more information about the energy equation see Eq (17) at https://www.comsol.com/multiphysics/heat-transfer-conservation-of-energy
subroutine temp_calc_ij(i,j, SS, dT_real, Pi, Hi_real, dhdt, ROi, EDAxi, NYs, umcxi,umcxi0,umcxi1, umcyj, umcyj0, umcyj1, u_slip , &
     tempi0, tempi, tempi1, tempj0, tempj1, temp0, ROi0,Hi0_real,Hi1_real,Hj0_real,Hj1_real, dROdTemp, dpdx,dpdy,t, dpdt,  temp_incre)
    implicit none
    include     'inc_CurrentT_met.h' 
    include     'inc_Grid.h'
    include     'inc_Glass_const.h'
    include     'inc_Outp.h'
    include     'inc_Rho.h'
    include     'inc_Ref.h'
    include     'inc_Shear_lim.h'
    include     'inc_Therm_param.h'
    include     'inc_Therm_matr.h'
	include     'inc_Therm_matr_past.h'
    include     'inc_Temp_reduction.h'
    include     'inc_Visc.h'
    ! Input
    integer      t
    real        Pi, EDAxi
    real        RO_r(5),    RO0, ROi
    real        P_press, dpdx, dpdy
    real        tempi0, tempi, tempi1, tempj0, tempj1, temp0
    integer     cont
    real        Vi, Vi0, ROi0, Hi_real, Hi0_real,Hi1_real,Hj0_real,Hj1_real, partition
    real        dROdTemp, dpdt, dROdp, p_oil_y, P_oil_x
    real        dhdt, umcxi0, umcxi1, umcyj0, umcyj1,    umcxi, umcyj, dRodP_real
    ! Calculations
    integer     I, J, I0, I1, J0, J1
    integer     ii, jj, ii1, jj1, ii0, jj0
    integer     h_oil_c
    real        P_tau, P_met, P_oil_f, P_oil_b, P_oil_r, P_oil_l, ptot, P_vol
    real        P_trans_i, P_trans_o, P_trans_iy, P_trans_oy, P_tau_press, P_shear
    real        dt_real, u_slip, tau_x,  temp_incre
    integer     ss, NYs, im, jm, NN, NX_m, NY_m, NN_m, im1, jm1
    real        temp_mai, temp_mbi, Area
    real        T_sa, T_sb, dz, P_comp, DXSSB
    real        K_oil_i1, K_oil_i0, K_oil_j0, K_oil_j1, CP_oil_i1, CP_oil_i0, CP_oil_j0, CP_oil_j1, dCp_dt
    real        E_dens, E_cp, E_temp, P_term, P_trans, P_cond
    real        tau_a, tau_b, tau_gamma_r, tau_lim_c, z_tau
    !other
    real        test, test2
    
    NN=(NYs+1)/2
    NY_m=(nys+1)/2
    NN_m=(NY_m+1)/2
    NX_m=(NX+1)/2                                                           !D_IN: /Grid/ -> Nx
    
    J0=J-SS
    IF( J0 .LE. 0)  J0=J+SS
    J1=J+1*SS
    IF( J1 .GT. NN) J1=NN-SS
                
    I0=I-1*SS
    IF( I0 .LE. 0)  I0=1
    I1=I+1*SS
    IF( I1 .GT. NX) I1=NX
    
    DXSSB=(SS*DX*B)                                                         !D_IN: /Grid/ -> DX, /outp/ -> B
    Area=DXSSB*DXSSB
    dz=2*DXSSB
    h_oil_c=4
    
    ! To calculate dRodP
    CALL Epsilon_derivative(i,j,i,j, 1.0, 1.0, 1.0, 1.0, 1.0, ss, test, dRodP, test2)             !CALL: No comment
    
    dRodP_real=ROi*dRodP/Ph ! Since epsilon_derivative calculates in dimensionless form and here all dimensions are added   !D_IN: /Visc/ -> PH
    
    ! The metal has twiz the stepsize between nodes as the oil. This makes it more stable in time dep simulations
    IF (MOD(i-1+ss,2*SS)==0 .AND. MOD(j-1+ss,2*SS)==0 ) THEN
            !THE NUMBER IS AN EVEN NUMBER
            ! For nodes of temp_metal is needed
            im=i/2
            jm=j/2
            
            im1=im+ss
            if( im1 .GT. NX_m) Im1=NX_m-ss
            jm1=jm+ss
            if( jm1 .GT. NN_m) jm1=NN_m-ss
            
            temp_mai    = (temp_ma(im,jm,1) +temp_mb(im1,jm,1) +temp_mb(im1,jm1,1)+ temp_ma(im,jm1,1))/4                    !D_IN: /CurrentT_met/ -> temp_ma, temp_mb
            temp_mbi    = (temp_mb(im,jm,1) +temp_mb(im1,jm,1) +temp_mb(im1,jm1,1)+ temp_mb(im,jm1,1))/4

    ELSEIF(MOD(i-1+ss,2*SS)==0 )THEN
            ! Two nodes of temp_metal is needed in x-direction
            im=i/2
            jm=(j-1)/2+1   
            im1=im+ss
            if( im1 .GT. NX_m) Im1=NX_m-ss
            
            temp_mai    = (temp_ma(im,jm,1)+ temp_ma(im1,jm,1))/2
            temp_mbi    = (temp_mb(im,jm,1)+ temp_mb(im1,jm,1))/2

        
    ELSEIF(MOD(j-1+ss,2*SS)==0 )THEN
            ! Two nodes of temp_metal is needed in y-direction
            im=(i-1)/2+1
            jm=(j)/2
            jm1=jm+ss
            if( jm1 .GT. NN_m) jm1=NN_m-ss
            
            temp_mai    = (temp_ma(im,jm,1)+ temp_ma(im,jm1,1))/2
            temp_mbi    = (temp_mb(im,jm,1)+ temp_mb(im,jm1,1))/2
            
            
    ELSE
            ! Only one temp_metal node is needed
            im=(i-1)/2+1
            jm=(j-1)/2+1
            
            temp_mai    = temp_ma(im,jm,1)
            temp_mbi    = temp_mb(im,jm,1)
 
    END IF
    
    If(lub_param==4) Then                                                                                                   !D_IN: /Ref/ -> lub_param
        ii=i*2-1
        jj=j*2-1
        
        ii1=i*2
        jj1=j*2
        
        ii0=i*2-2*ss
        if(ii0 .LT. 1) ii0=1
        jj0=j*2-2*ss
        if(jj0 .LT. 1) jj0=1
        
        K_oil   = K_O(ii,jj)                                                                                                !D_OUT: K_oil -> /Therm_param/ (Not saved)                                       !D_IN: /Therm_matr/ -> K_o
        K_oil_i1 = K_O(ii1,jj)
        K_oil_i0 = K_O(ii0,jj)
        K_oil_j0 = K_O(ii,jj0)
        K_oil_j1 = K_O(ii,jj1)
        
        Cp_oil   = Cp_O(ii,jj)                                                                                              !D_OUT: Cp_oil -> /Therm_param/ (Not saved)                                      !D_IN: /Therm_matr/ -> Cp_oil
        Cp_oil_i1 = Cp_O(ii1,jj)
        Cp_oil_i0 = Cp_O(ii0,jj)
        Cp_oil_j0 = Cp_O(ii,jj0)
        Cp_oil_j1 = Cp_O(ii,jj1)
        
        if (t .GT. -11) then
            dCp_dt   = (Cp_O(ii,jj)-Cp_O_past(ii,jj))/dt_real                                                               !D_IN: /Therm_matr_past/ -> CP_O_past 
        else
            dCp_dt   = 0.0
        endif
        
    else
        K_oil   = K_oil
        K_oil_i1 = K_oil
        K_oil_i0 = K_oil
        K_oil_j0 = K_oil
        K_oil_j1 = K_oil
        
        Cp_oil   = Cp_oil
        Cp_oil_i1 = Cp_oil
        Cp_oil_i0 = Cp_oil
        Cp_oil_j0 = Cp_oil
        Cp_oil_j1 = Cp_oil
        
        dCp_dt = 0.0
    endif
    if(tempi .gt. 400) then
        test=tempi
    endif
    
                    ! Power from pressure gradient ---------------------------------------------------
                    P_term=hi_real**2 /(12.0*EDAxi)                                     ! [m^2/(N/m^2)*s = m^4/Ns]
                    IF(P_term .LT. 0) then
                        WRITE(4,*)' Waringing, Bad P_term value. P_term was =', P_term
                        P_term=abs(P_term)
                    endif
                    P_press=(dpdx**2+dpdy**2)*P_term                                    ! [N^2/(m^4*m^2) * m^4/Ns = N/(m^2*s)] 
                    if(temp_param==4) then                                              !D_IN: /Temp_reduction/ -> temp_param
                        P_press = P_press * temp_fac                                    !D_IN: /Temp_reduction/ -> temp_fac
                    endif
                    
                    ! Power from slip ----------------------------------------------------------------------
                    tau_x   = EDAXi*u_slip/hi_real                                      !  [N/m^2*s *m/s/m = N/m^2]
                    P_tau    =tau_x*u_slip/hi_real                                      ! [N/m^2 *m/s/m = N/(m^2*s)]
                    If( P_tau .LT. 0) then
                        WRITE(4,*)'Negative power from tau'
                        WRITE(4,*)'P_tau = ', P_tau, ', P_tau = ', P_tau, ', tau_x = ', tau_x, ', u_slip = ', u_slip
                        WRITE(4,*)'i,j = ', i,j, ',  hi_real = ', hi_real,', tempi = ', tempi
                        STOP
                    endif
                    
                    ! The combination term from shear and pressure gradient ----------------------
                    P_tau_press = abs(dpdx * u_slip / 2)                                ! [N/m^2/m * m/s = N/(m^2*s)] Gives contributions no mather what direction slip is defined in
                    
                    ! The total power from shear
                    P_shear = P_tau + P_press + P_tau_press
                    
                    If( P_shear .LT. 0) then
                        WRITE(4,*)'Negative total power from shearing'
                        WRITE(4,*)'P_tau = ', P_tau, ', P_press = ', P_press, ', P_tau_press = ', P_tau_press
                        WRITE(4,*)'i,j = ', i,j, ',  hi_real = ', hi_real,', tempi = ', tempi, 'dpdx', dpdx, 'u_slip',u_slip
                        P_shear = abs( P_shear ) ! Power is still generated 
                    endif
                    
                    
                    ! Power from pressure change on adiabatic compression  -------------------------------------
                    P_comp = umcxi*dpdx+umcyj*dpdy                ! [1 * m/s * N/(m^2*m) = N/(m^2*s)]
                    
                    if( t .GT. -4) then             ! Add the temperature from time dependent compressuon
                        P_comp=P_comp + (dpdt)
                    elseif( t .GT. -6) then
                        P_comp=P_comp + (dpdt)*0.5  ! Add it a bit earlier. dpdt shold anyway be zero at time indep parts
                    endif
                    P_comp= P_comp * tempi/ROi*dROdTemp

                    ! Power from heat flow due to material conductivity -----------------------------------------
                    if( solid_material==2)then                                                                                                          !D_IN: /Glass_const/ -> solid_material
                        T_sb         = (tempi*K_oil/(hi_real/h_oil_c)+temp_mbi*K_glass/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_glass/dz)                 !D_IN: /Glass_const/ -> K_glass
                        T_sa         = (tempi*K_oil/(hi_real/h_oil_c)+temp_mai*K_met/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_met/dz)
                    else if( solid_material==3)then
                        T_sa         = (tempi*K_oil/(hi_real/h_oil_c)+temp_mbi*K_glass/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_glass/dz)
                        T_sb         = (tempi*K_oil/(hi_real/h_oil_c)+temp_mbi*K_met/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_met/dz)
                    else
                        T_sb         = (tempi*K_oil/(hi_real/h_oil_c)+temp_mbi*K_met/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_met/dz)
                        T_sa         = (tempi*K_oil/(hi_real/h_oil_c)+temp_mai*K_met/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_met/dz)
                    endif
                    
                    ! Power flow down into metal
                    P_met       = (T_sa-2*tempi+T_sb)*K_oil/(hi_real/h_oil_c)**2
                    
                    ! Power flow in oil
                    if( i == 1 .or. i == Nx) then                                                                                                       
                        P_oil_x=0
                    else
                        P_oil_x     = ( (K_oil_i1+K_oil)*tempi1 - (2*K_oil+K_oil_i1+K_oil_i0)*tempi + (K_oil_i0+K_oil)*tempi0 ) / (2 *DXSSB**2) 
                    endif
                    
                    IF( J == 1 )then
                        P_oil_y= 0
                    elseif( J == NN ) then
                        P_oil_y     = ( (K_oil_j0+K_oil)*tempj0 - (2*K_oil+K_oil_j0+K_oil_j0)*tempi + (K_oil_j0+K_oil)*tempj0 ) / (2 *DXSSB**2)  ! due to symmetry
                    else
                        P_oil_y     = ( (K_oil_j1+K_oil)*tempj1 - (2*K_oil+K_oil_j1+K_oil_j0)*tempi + (K_oil_j0+K_oil)*tempj0 ) / (2 *DXSSB**2) 
                    endif
                    
                    P_cond      =  P_oil_x + P_oil_y + P_met
                    
                    ! Power transportation -------------------------------------------------------------------------------------------------
                    
                    ! Flow in x-direction
                    P_trans_i     =  ROi *CP_oil*0.5*(tempi1-tempi0)/DXSSB*umcxi ! Energy transported  along with the flow
                    P_trans_o     =  0 ! DXSSB*0.5*(hi_real+hi1_real) *ROi *CP_oil*0.5*(tempi+tempi1)*0.5*(umcxi+umcxi1)
                    
                    ! The central method has the drawback that if the current node is warmer than both the previus and the next, it will still heat up if the derivative is negative. 
                    IF( tempi1 .GE. tempi .and. tempi .GE. tempi0) then
                        continue ! Since then central difference works
                    elseIF( tempi1 .LE. tempi .and. tempi .LE. tempi0) then   
                        continue ! Since then central difference works
                    else
                        ! Central difference does not work
                            if(umcxi .ge. 0)then
                                P_trans_i = ROi *CP_oil*(tempi-tempi0)/DXSSB*umcxi
                            else
                                P_trans_i = ROi *CP_oil*(tempi1-tempi)/DXSSB*umcxi
                            endif
                    endif
    
                    
                    IF( I==NX) then
                        P_trans_o =  ROi *CP_oil*(tempi-tempi0)/DXSSB*umcxi  ! DXSSB*(hi_real) *ROi *CP_oil*(tempi)*(umcxi)
                    endif

                    ! Flow in y-direction 
                    P_trans_iy     =  ROi *CP_oil*0.5*(tempj1-tempj0)/DXSSB*umcyj
                    P_trans_oy     =  0 ! DXSSB*0.5*(hi_real+hj1_real) *ROi *CP_oil*0.5*(tempi+tempj1)*0.5*(umcyj+umcyj1)
                    
                    !IF( tempi .GT. max(tempj1,tempj0)) then
                    IF( tempj1 .GE. tempi .and. tempi .GE. tempj0) then
                        continue ! Since then central difference works
                    elseIF( tempj1 .LE. tempi .and. tempi .LE. tempj0) then   
                        continue ! Since then central difference works
                    else
                            if(umcyj .ge. 0)then
                                P_trans_iy = ROi *CP_oil*(tempi-tempj0)/DXSSB*umcyj
                            else
                                P_trans_iy = ROi *CP_oil*(tempj1-tempi)/DXSSB*umcyj
                            endif
                    endif
                    
                    IF( J == 1 .or. j == NN) then
                        P_trans_iy =  0 ! Due to symmetry. No transportation of energy except via volume changes which is not looked at directly. Conduction still occurs
                        P_trans_oy =  0 ! -P_trans_iy ! Due to symetry and a minus sign in the Ptot equation
                    endif
                    
                    P_trans          = P_trans_i + P_trans_o + P_trans_iy + P_trans_oy
                    !------------------------------------------------------------
                    ! Equation changing how much energy is stored per degree-------------------------------------------------------------
                    !----------------------------------------------------------------
                    
                    ! Density change
                    E_dens  = CP_oil*tempi*dROdTemp
                    
                    ! Cp hange
                    E_cp    = ROi*tempi*DCP_Dtemp(i,j)                         !D_IN: /Therm_matr/ -> DCP_Dtemp(I,J)
                    
                    ! Temp change
                    E_temp  = ROi*Cp_oil

                    Ptot        = P_shear - P_comp - P_trans + P_cond
                    
                    IF( P_tau .LT. 0 .or. P_press .LT. 0) then
                        WRITE(4,*)'Negative power from shearing'
                        WRITE(4,*)'P_tau = ', P_tau, ', P_press = ', P_press
                        WRITE(4,*)'i,j = ', i,j, ',  hi_real = ', hi_real,', tempi = ', tempi
                        STOP
                    endif

                    temp_incre= Ptot/(E_temp + E_Cp + E_dens)
                    
                    if( abs(temp_incre) .gt. 1e14) then
                        WRITE(4,*)'High temp incre in oil'
                        if( abs(P_shear) .gt. 1e13)then
                            WRITE(4,*) 'P_shear=', P_shear
                            WRITE(4,*) 'P_tau', P_tau,'P_press',P_press,' P_tau_press', P_tau_press
                            WRITE(4,*) 'P_term ', P_term , 'dpdx ', dpdx, 'dpdy ', dpdy
                            WRITE(4,*) 'temp_fac ', temp_fac
                        endif
                        
                        if( abs(P_comp) .gt. 1e13) WRITE(4,*) 'P_comp=', P_comp
                        
                        if( abs(P_trans) .gt. 1e13) Then
                            WRITE(4,*) 'P_trans=', P_trans
                            WRITE(4,*) 'umcxi', umcxi, 'umcyj', umcyj
                        endif
                        
                        if( abs(P_cond) .gt. 1e13) WRITE(4,*) 'P_cond=', P_cond
                        
                    endif
    
                IF( isnan(temp_incre )) then
                    WRITE(4,*)'Nan temperature increment in Oil'
                    Write(4,*)'For I J = ', I ,J , 'temp(i,j) = ',tempi
                    STOP
                endif

        return
    end
    