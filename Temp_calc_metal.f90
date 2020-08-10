    ! A subroutine to calculate the temperature increments insde metal domains. The temperature fields are then updated in the subroutine Metal_temp_upd
    subroutine temp_calc_metal(dt_lim, SS, t, NYs, temp_convm)
    implicit none
    ! Input
    include     'inc_DZ_com.h'
    include     'inc_Contact_mat.h'
	include     'inc_CTRA4.h'
	include     'inc_Current.h'
	include     'inc_CurrentH.h'
	include     'inc_CurrentP.h'
	include     'inc_CurrentT.h'
    include     'inc_CurrentT_met.h' 
    include     'inc_Glass_const.h'
    include     'inc_Grid.h'
    include     'inc_Met_Temp_inc.h'
    include     'inc_Outp.h'
    include     'inc_PastT_met.h'
    include     'inc_Rho.h'
	include     'inc_Ref.h'
    include     'inc_Shear_lim.h'
    include     'inc_Therm_param.h'
    include     'inc_Therm_matr.h'
    include     'inc_Temp_conduction.h'
    include     'inc_Visc.h'  
    include     'inc_Cores_used.h'
    integer      t, ss, NYs
    real        dt_lim
    ! Calculations
    integer     t_lim, n_met
    integer     I, J, I0, I1, J0, J1, k0, k1
    integer     Im, Jm, I0m, I1m, J0m, J1m
    integer     NX_m,NY_m, NN_m, i2, j2
    integer     k, l, start, fini
    real        P_met_f, P_met_b, P_met_r, P_met_l, ptot, P_trans_i, P_trans_o, p_met_out, P_met_in
    real, pointer ::    T_met(:,:,:)
    real        dt_real, u_slip, hi_real, new_temp, test
    integer     NN
    real        speed, umax, ddt, t_diff
    real        DX_met, DZ, frack, P_tau, tau_x, T_s
    real        time
    integer     fin, h_oil_c
    real        test1, test2, test3, test4, test5, test6
    real        K_m, CP_m, rho_m    
    real        DT_lim_c, DT_nodes
    real        P_trans, P_cond_x, P_cond_y, P_cond_z
    real        temp_convm, temp_conv_v
    ! Output
    save  /Met_Temp_inc/
    save  /CurrentT_met/
 
    u_slip  = abs(ua-ub)      !D_IN: /Outp/ -> Ua, Ub
    
    NN      = (NYs+1)/2
    
    DX_met      = 2*DX*SS*B ! A 2 since these grid in the metals are 2 times courser than the grid in the lubricant   !D_IN: /Grid/ -> Dx; /Outp/ -> B
    DZ=DZ_vec(1)            ! Vertical distance dependes on DZ_method    !D_IN: /DZ_com/ -> Dz_vec(I)
    
    ! Material data
    n_met   = 39      ! Number of nodelayers of metal
    
    ! The metal has twize the stepsize between nodes as the oil. To make it more stable in time dep simulations
    NX_m    = (nx+1)/2    !D_IN: /Grid/ -> Nx 
    NY_m    = (nys+1)/2   !D_IN: /Grid/ -> Ny
    NN_m    = (NY_m+1)/2
    
    !!$OMP PARALLEL DO PRIVATE(l,speed,h_oil_c,k,j,T_met,K_m,CP_m,rho_m,start,fini) DEFAULT(SHARED)
    do l= 1,2                    ! Do twice, once fore each contacting body 
        K_m     = K_met          !D_IN: /Therm_param/ -> K_met
        CP_m    = CP_met         !D_IN: /Therm_param/ -> CP_met
        rho_m   = rho_met        !D_IN: /Therm_param/ -> rho_met

        if( l==1)Then ! Surface with surface iregularities
            speed   = ua
            IF( Shear_band == 0) then
                h_oil_c = 4 ! Old setup untill 7001c
            else  ! Opposit surface
                IF(ua == ub) then   ! Shear_xi has to be 0.5 since both surfaces moves with the same speed
                    h_oil_c = 1/(0.5 / Shear_band) 
                else if(ua .gt. ub) then ! we're near the faster surface
                    h_oil_c = 1/( Shear_xi / Shear_band) ! Shear band =1 means that the energy is transported from the shear band. Shear band= 2 means that the energy is spread out in the lubricant
                else
                    h_oil_c = 1/ ((1-Shear_xi) / Shear_band)
                endif
            endif

            do k=1,n_met,ss
                do j=1,NN_m,ss
                    if(speed .GE. 0) then
                        Temp_ma(1,j,k)=temp(1,1)                                                                                                                !D_IN: /CurrentT/ -> temp(I,J) !D_OUT: Temp_ma -> /CurrentT_met/
                        if(Temp_ma(NX_m,j,k)  .GT. Temp_ma(NX_m-ss,j,k))    Temp_ma(NX_m,j,k)=Temp_ma(NX_m-ss,j,k)  !At the boundary the heat has to disipate.
                        if(Temp_ma(NX_m,j,k)  .LT. temp(1,1))               Temp_ma(NX_m,j,k)=temp(1,1)
                    else
                        Temp_ma(NX_m,j,k)=temp(1,1)
                        if(Temp_ma(1,j,k)  .GT. Temp_ma(1+ss,j,k))          Temp_ma(1,j,k)=Temp_ma(1+ss,j,k)  !At the boundary the heat has to disipate.
                        if(Temp_ma(1,j,k)  .LT. temp(1,1))                  Temp_ma(1,j,k)=temp(1,1)
                    endif
                enddo
            enddo
            
            T_met   => Temp_ma
            
            if(solid_material==3)then                                   !D_IN: /Glass_const/ -> solid_material
                K_m     = K_glass                                       !D_IN: /Glass_const/ -> K_glass
                CP_m    = CP_glass                                      !D_IN: /Glass_const/ -> CP_glass
                rho_m   = rho_glass                                     !D_IN: /Glass_const/ -> Rho_glass
            endif

        else
            speed   = ub
            IF( Shear_band == 0) then
                h_oil_c = 4 ! Old setup untill 7001c
            else
                IF(ub == ua) then   ! Shear_xi has to be 0.5 since both surfaces moves with the same speed
                    h_oil_c = 1/(0.5 / Shear_band) 
                else if(ub .gt. ua) then ! we're near the faster surface
                    h_oil_c = 1/(Shear_xi / Shear_band) ! Shear band =1 means that the energy is transported from the shear band. Shear band= 2 means that the energy is spread out in the lubricant
                else
                    h_oil_c = 1/((1-Shear_xi) / Shear_band)
                endif
            endif
            
            ! Enforce boundary conditions
            do k=1,n_met,ss
                do j=1,NN_m,ss
                    if(speed .GE. 0) then
                        Temp_mb(1,j,k)=temp(1,1)
                        if(Temp_mb(NX_m,j,k)  .GT. Temp_mb(NX_m-ss,j,k))    Temp_mb(NX_m,j,k)=Temp_mb(NX_m-ss,j,k)  !At the boundary the heat has to disipate.          !D_OUT: Temp_mb -> /CurrentT_met/
                        if(Temp_mb(NX_m,j,k)  .LT. temp(1,1))               Temp_mb(NX_m,j,k)=temp(1,1)
                    else
                        Temp_mb(NX_m,j,k)=temp(1,1)
                        if(Temp_mb(1,j,k)  .GT. Temp_mb(1+ss,j,k))          Temp_mb(1,j,k)=Temp_mb(1+ss,j,k)  !At the boundary the heat has to disipate.
                        if(Temp_mb(1,j,k)  .LT. temp(1,1))                  Temp_mb(1,j,k)=temp(1,1)
                    endif
                    
                enddo
            enddo
                
            T_met   => Temp_mb
            
            if(solid_material==2)then
                K_m     = K_glass
                CP_m    = CP_glass
                rho_m   = rho_glass
            endif
            
        endif
        
        ! Different boundary conditions depending on which direction the metal is moving
        if(speed .GT. 0) then
            start   = 1+ss
            fini    = NX_m
        elseif(speed == 0) then
            start   = 1
            fini    = NX_m
        else
            start   = 1
            fini    = NX_m-SS
        endif
        
        time=0                  ! Internal time
        fin=0                ! If run for enough time
        !dt_lim=dt_real      ! The value of 2 times the timestep
        
    ! ---------------- k=1 ---------------------------------------------------------------------------
        k=1
        !$OMP PARALLEL DO IF(use_multiple_cores)                                                                                    &
        !$OMP&            PRIVATE(Jm,J0m,J1m,Im,I0m,I1m,P_met_out,P_tau,j2,j,i2,frack,i,j0,j1,I0,I1,hi_real,tau_x,                  &
        !$OMP&                      t_diff,T_s,P_met_in,P_cond_z,P_cond_x,P_cond_y,P_trans,Ptot,dT_nodes,dt_lim_c)                  &
        !$OMP&            FIRSTPRIVATE(K_oil)                                                                                       &
        !$OMP&            SHARED(NN_m,SS,start,fini,NX_m,NX,H,Rx,B,contact,P,ph,u_slip,DZ,Temp_ma,Temp_mb,k,K_m,NN,                 & 
        !$OMP&                      temp,T_met,lub_param,K_O,h_oil_c,DX_met,t,geom,speed,rho_m,CP_m,l,temp_incre,EDA0)              &
        !$OMP&            REDUCTION(min:dt_lim)
        DO Jm=1,NN_m,SS
            
            ! Define nodenumbers for past and next nodes
            ! For the metal
            J0m=Jm-SS
            J1m=Jm+1*SS
            IF( J1m .GT. NN_m)  J1m=NN_m-SS
            
            
                            
            DO Im=start,fini,SS
                    
                    ! Define nodenumbers for past and next nodes
                    ! For the metal
                    I0m=Im-SS
                    I1m=Im+SS
                    
                    IF(I0m .LT. 1) I0m = 1
                    IF(I1m .GT. NX_m) I1m=NX_m
                    
                    ! Conduction ------------------------------------------------------------------------------
                    ! z direction -----------------------------------------------------------------------------
                    ! For the oil
                    P_met_out   =0                             
                    P_tau       =0
  
                    do j2=-1,1
                        j=(jm-1)*2+1+(j2)*SS
                        if( j .Lt. 1) j=1
                        if( j .GT. NN) j=NN-ss
                        
                        do i2=-1,1
                            frack=0.25
                        
                            i=(im-1)*2+1+(i2)*SS
                            if( i .Lt. 1)  i=1
                            if( i .GT. NX) i=NX
                            
                            J0=J-SS
                            IF( J0 .LE. 0)  J0=J+SS
                            J1=J+1*SS
                            IF( J1 .GT. NN) J1=NN-SS
                            
                            I0=I-1*SS
                        
                            IF( I0 .LE. 0)  I0=1
                            I1=I+1*SS
                            IF( I1 .GE. NX) I1=NX

                            hi_real  =H(I,J)/Rx*B**2          !D_IN: /CurrentH/ -> H(i,j); /Outp/ -> Rx
                        
                            ! If larger connection area
                            if( i2==0) frack=frack*2 
                            if( j2==0) frack=frack*2
                        
                           ! If metal contact, the energy is generated in the metal and transported without the oil inbetween
                            if( contact(i,j)==2 .or. contact(i0,j)==2 .or. contact(i1,j)==2  .or. contact(i,j1)==2 .or. contact(i,j0)==2 ) then   !D_IN: /Contact_mat/ -> contact(i,j)    
                                tau_x       = P(i,j)*ph*0.3 ! Shear at centre of film     !D_IN: /CurrentP/ -> P(I,J)
                                P_tau       = P_tau+ tau_x*u_slip*frack*0.5 /(DZ)  ! Energy from shear. The energy is divided equally between the two surfaces
                                
                                if (l==1) then
                                    t_diff=  Temp_ma(im,jm,k) -Temp_mb(im,jm,k)
                                else
                                    t_diff= -Temp_ma(im,jm,k) +Temp_mb(im,jm,k)
                                endif
                                
                                P_met_out       = P_met_out+ (T_diff)*K_m*frack/Dz !*B**2
  
                            ! If not metal contact
                            else
                                if( temp(i,j) .NE. T_met(im,jm,k))then
                                    if(lub_param .EQ. 4)then             !D_IN: /Ref/ -> Lub_param
                                        ! The material parameters are already normalized with Ph in this IF statement
                                        k_oil   = K_O(i*2-1,j*2-1)       !D_IN: /Term_matr/ -> K_O(I,J) !D_OUT: k_oil -> /Therm_param/ (Not saved)
                                    else
                                        k_oil   = K_oil
                                    endif
                                                        
                                    T_s=(temp(i,j)*K_oil/(hi_real/h_oil_c)+T_met(im,jm,k)*K_m/(dz/2)) / (h_oil_c*K_oil/(hi_real)+2*K_m/dz)
                                    P_met_out       = P_met_out+ (T_s-temp(i,j))*K_oil/(hi_real/h_oil_c)*frack!*B**2         ! Energy cunducted outwards towards oil , A=DX^2, L=H/4
                                    
                                endif
                            
                            endif
                    
                        enddo
                    enddo
 
                    
                    P_met_in     = (T_met(im,jm,k)-T_met(im,jm,k+SS))*K_m/DZ ! Energy cunducted into metal
                    
                    P_cond_z     = (-P_met_in-P_met_out)/DZ ! The change in first derivative = the second derivative
                    
                    ! x- direction ---------------------------------------------------------------------
  
                    if( im==NX_m .or. im==1 ) Then
                        P_cond_x     = 0
                    else
                        P_cond_x  = (T_met(i1m,jm,k)-2*T_met(im,jm,k) +T_met(i0m,jm,k)  )*K_m / (DX_met**2) ! Energy cunducted outwards, towards the oil
                    endif
                    
                    ! y- direction --------------------------------------------------------------------
                    
                    if( t .LT. 0 .and. geom .eq. 1) then            !D_IN: /Ref/ -> geom
                        P_cond_y  = 0.0                                                                 ! Since no y-derivatives
                    elseif( jm==NN_m ) then
                        P_cond_y  = (T_met(im,j0m,k)-2*T_met(im,jm,k) +T_met(im,j0m,k)  )*K_m / (DX_met**2) ! Due to symmetry
                    elseif( jm==1 ) Then
                        P_cond_y  = (T_met(im,j1m,k)-2*T_met(im,jm,k) +T_met(im,j1m,k)  )*K_m / (DX_met**2) ! Due to symmetry
                    else
                         P_cond_y  = (T_met(im,j1m,k)-2*T_met(im,jm,k) +T_met(im,j0m,k)  )*K_m / (DX_met**2) ! Energy cunducted outwards, towards the oil
                    endif
                    
                    ! Transportation --------------------------------------------------------------------------------
                    if( speed==0)then
                        P_trans=0
                    else
                        if(im==1 ) then
                            P_trans  = speed*rho_m*CP_m*( (T_met(i1m,jm,k)-T_met(im,jm,k)) / (DX_met) + (T_met(im+2*ss,jm,k)-T_met(im,jm,k)) / (2*DX_met) ) *0.5
                            !if(P_trans .GT. 0) P_trans=0  ! No input of energy at outlet
                            
                        elseif( im==Nx_m)then
                            P_trans  = speed*rho_m*CP_m*( (T_met(im,jm,k)-T_met(i0m,jm,k)) / (DX_met) + (T_met(im,jm,k)-T_met(im -2*ss,jm,k)) / (2*DX_met) ) *0.5
                            !if(P_trans .GT. 0) P_trans=0  ! No input of energy at outlet
                            
                        else
                            P_trans    = speed*rho_m*CP_m*(T_met(i1m,jm,k)-T_met(i0m,jm,k)) / (2*DX_met)! Energy transported in along with the flow
                        endif
                        
                    endif

                    Ptot =  P_tau - P_trans + P_cond_x + P_cond_y + P_cond_z
                    temp_incre(im,jm,k,l)    =  Ptot/(CP_m*rho_m )  !D_OUT: temp_incre(I,J) -> /Met_Temp_inc/ 

                    ! Checking for limit on time increment
                    IF( Im==NX_m) then
                        dT_nodes=abs(T_met(im,jm,k)-T_met(i0m,jm,k))
                    else
                        dT_nodes=abs(T_met(im,jm,k)-T_met(i1m,jm,k))
                    endif
                    
                    if( dT_nodes .GT. 1 ) then
                        dt_lim_c = abs(dT_nodes/temp_incre(im,jm,k,l)*0.15)
                        if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                    else
                        dt_lim_c = abs(1/temp_incre(im,jm,k,l))
                        if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                    endif
                    
                    IF( isnan(temp_incre(im,jm,k,l) )) then
                        !$OMP CRITICAL (C_1)
                        WRITE(4,*)'Nan temperature increment in Metal'
                        Write(4,*)'For Im Jm = ', Im ,Jm , ' k= ', k, ' P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp(i,j) = ',temp(i,j), 'T_met(i,j,k) = ',T_met(im,jm,k)   !D_IN: /Visc/ -> EDA0
                        !$OMP END CRITICAL (C_1)
                        STOP
                    endif
 
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
    
    !------- k=2:N_met -----------------------------------------------------------------------------------------------------------
        !$OMP PARALLEL DO IF(use_multiple_cores)                                                                      &
        !$OMP&            PRIVATE(K,K0,K1,J,J0,J1,I,I0,I1,P_cond_x,P_cond_y,P_cond_z,P_trans,Ptot,dT_nodes,dt_lim_c)  &
        !$OMP&            SHARED(n_met,SS,NN_m,start,fini,NX_m,T_met,K_m,DX_met,t,geom,temp,DZ_vec,speed,rho_m,CP_m,  &
        !$OMP&                   temp_incre,l,P,H,EDA0)                                                               &
        !$OMP&            REDUCTION(min:dt_lim)
        do K= 1+SS,n_met,SS      ! Number of nodelayers of metal
            K0=K-SS
            K1=K+SS
            IF (K1 .gt. N_met) K1=N_met
            
            DO J=1,NN_m,SS
                ! Define nodenumbers for past and next nodes
                J0=J-SS
                J1=J+1*SS
                IF( J1 .GT. NN_m) J1=NN_m-SS
            
                DO I=start,fini,SS
                    I0=I-1*SS
                    IF(I0 .LT. 1) I0 = 1
                    I1=I+1*SS
                    IF(I1 .GT. NX_m) I1 = NX_m

                    ! x- direction ---------------------------------------------------------------------
                        
                    if( i==Nx_m .or. i==1 ) Then
                        P_cond_x     = 0
                    else
                        P_cond_x  = (T_met(i1,j,k)-2*T_met(i,j,k) +T_met(i0,j,k)  )*K_m / (DX_met**2) ! Energy cunducted outwards, towards the oil
                    endif
                    
                    ! y- direction ---------------------------------------------------------------------
                        
                    if( t .LT. 0 .and. geom .eq. 1) then
                        P_cond_y  = 0.0                                                                 ! Since no y-derivatives
                    elseif( j==NN_m)then
                        P_cond_y  = (T_met(i,j0,k)-2*T_met(i,j,k) +T_met(i,j0,k)  )*K_m / (DX_met**2)       ! Due to symmetry
                    elseif(j==1 ) Then
                        P_cond_y  = (T_met(i,j1,k)-2*T_met(i,j,k) +T_met(i,j1,k)  )*K_m / (DX_met**2)
                    else
                        P_cond_y  = (T_met(i,j1,k)-2*T_met(i,j,k) +T_met(i,j0,k)  )*K_m / (DX_met**2)       ! Energy cunducted outwards, towards the oil
                    endif
                    
                    ! z- direction ---------------------------------------------------------------------
                        
                    if( k==n_met) Then
                        P_cond_z  = ( (temp(1,1)    -T_met(i,j,k))/DZ_vec(k1) - (T_met(i,j,k) - T_met(i,j,k0))/DZ_vec(K) )*K_m / ( (DZ_vec(K1)+DZ_vec(K))*0.5 ) ! The temperature outside in z-direction is set to be fixed. Since othervise the corner nodes will be overconstrained. 
                    else
                        P_cond_z  = ( (T_met(i,j,k1)-T_met(i,j,k))/DZ_vec(k1) - (T_met(i,j,k) - T_met(i,j,k0))/DZ_vec(K) )*K_m / ( (DZ_vec(K1)+DZ_vec(K))*0.5 ) ! Energy cunducted outwards, towards the oil
                    endif
                    
                    
                    ! Transportation ----------------------------------------------------------------------
                    if( speed   == 0 )then
                        P_trans = 0

                    else 
                        if(i==1 ) then
                            if(t_met(i,j,k) .GT. T_met(i1,j,k)) T_met(i,j,k)=T_met(i1,j,k)  !At the boundary the heat has to disipate. 
                            if(t_met(i,j,k) .LT. temp(1,1)) T_met(i,j,k)=temp(1,1)
                            P_trans  = speed*rho_m*CP_m* ( (T_met(i1,j,k)-T_met(i,j,k)) / (DX_met) + (T_met(i+2*ss,j,k)-T_met(i,j,k)) / (2*DX_met) ) *0.5
                            !if(P_trans .GT. 0) P_trans=0  ! No input of energy at outlet
                            
                        elseif( i==Nx_m)then
                            if(t_met(i,j,k) .GT. T_met(i0,j,k)) T_met(i,j,k)=T_met(i0,j,k)  !At the boundary the heat has to disipate.
                            if(t_met(i,j,k) .LT. temp(1,1)) T_met(i,j,k)=temp(1,1)
                            P_trans  = speed*rho_m*CP_m* ( (T_met(i,j,k)-T_met(i0,j,k)) / (DX_met) + (T_met(i,j,k)-T_met(i-2*ss,j,k)) / (2*DX_met) ) *0.5
                            !if(P_trans .GT. 0) P_trans=0  ! No input of energy at outlet
                            
                        else
                            P_trans  = speed*rho_m*CP_m*(T_met(i1,j,k)-T_met(i0,j,k)) / (2*DX_met)
                        endif
                        
                    endif
                    
                    Ptot = - P_trans + (P_cond_x + P_cond_y + P_cond_z ) 
                    temp_incre(i,j,k,l)    =  Ptot/(CP_m*rho_m )
                    
                    ! Checking for limit on time increment
                    IF( I==NX_m) then
                        dT_nodes=abs(T_met(i,j,k)-T_met(i0,j,k))
                    else
                        dT_nodes=abs(T_met(i,j,k)-T_met(i1,j,k))
                    endif
                    
                    
                    if( dT_nodes .GT. 1 ) then
                        dt_lim_c = abs(dT_nodes/temp_incre(i,j,k,l)*0.15)
                        if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                    else
                        dt_lim_c = abs(1/temp_incre(i,j,k,l))
                        if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                    endif
                    
                    IF( isnan(temp_incre(i,j,k,l) )) then
                        !$OMP CRITICAL (C_2)
                        WRITE(4,*)'Nan temperature increment in Metal'
                        Write(4,*)'For I J = ', I ,J , ' k= ', k, ' P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp(i,j) = ',temp(i,j), 'T_met(i,j,k) = ',T_met(im,jm,k)
                        !$OMP END CRITICAL (C_2)
                        STOP
                    endif
                    
                ENDDO
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
    enddo
    !!$OMP END PARALLEL DO    
    
    temp_convm=0.0
    
    Do i=1,5
        temp_conv_v     = sum(abs(temp_incre(:,NN_m,i,1)))/NX_m
        temp_convm      = max(temp_convm,temp_conv_v)                          !??? Do those go anywhere? !D_OUT: temp_convm -> returned
    Enddo
 
        
        return
    end
    