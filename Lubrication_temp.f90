! A set of subroitines are used to calculate the temperature in the lubricant and the metal bodies
    ! The subroutine temp_calc_ij calculates the temperature increments inside the lubricant one node at each call 
    ! The Subroutine temp_calc_metal calculates the temperature increments in the metal bodies
    ! The subroutine Metal_temp_upd updates the temperature fields inside the metal based on the increments generated by temp_calc_metal
    Subroutine Lubrication_temp(SS, NYs, C1, t, sim_end, sim_pass, iter_lim, iter_cvot, temp_conv, temp_convm)
        implicit none 
        include     'inc_Average_past.h'
        include     'inc_Contact_mat.h'
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_CurrentT_met.h'
        include     'inc_CurrentRO.h'
        include     'inc_C_method.h'
        include     'inc_EDA_average.h'
        include     'inc_Grid.h'
        include     'inc_G0DT.h'
        include     'inc_Holmes.h'
        include     'inc_Method.h'
        include     'inc_Outp.h'
        include     'inc_Past.h'
        include     'inc_PastT.h'
        include     'inc_PastT_met.h'
        include     'inc_Rho.h'
        include     'inc_RLarsson.h'
        include     'inc_Ref.h'
        include     'inc_Setup.h'
	    include     'inc_Shear_lim.h'
	    include     'inc_Therm_param.h'
        include     'inc_Therm_matr.h'
        include     'inc_Therm_matr_past.h'
        include     'inc_Glass_const.h'
        include     'inc_Temp_reduction.h'
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'
	    include     'inc_Y_liu.h'
        include     'inc_Cores_used.h'                                                           
        ! Input
        real        C1
        integer     NN, SS, t, sim_end, sim_pass, iter_lim, iter_cvot, M
        real        T40, Tcvot
        ! Calculations
        integer     err, iter
        real        u_slip
        real        hi_real
        real        umc
        real        b2Rx, b2RxbUm, DXSSPH, ddt, dt_real
        integer     I, J, I0, I1, I3, J0, J1, JJ, I2, J2, J3, im, jm
        real        umax, temp_1, temp_2
        real        tempi0, tempi, tempi1, tempj0, tempj1, time, temp_incre, dt_lim, temp00
        integer     fin
        real        DXSSB, temp_conv, temp_convm
        real        dt_nodes, dt_lim_c 
        integer     NX_m, NY_m, NN_m, im1, jm1
        real        umcxi, umcxi0, umcxi1, umcyj, umcyj0, umcyj1
        real        Hi0_real,Hi1_real,Hj0_real,Hj1_real
        real        Length, cvot
        real        P_lim1, P_lim2, diffP
        real        epst, u_c_max, temp_count_max, temp_ave
        integer     u_c_number, temp_count
        integer, allocatable :: Cavitaion(:,:), cont(:,:)
        real, allocatable :: dPdx(:,:), dPdy(:,:), h_real_m(:,:), temp_i(:,:), umcx(:,:), umcy(:,:)  
        real, allocatable :: dROdTemp(:,:), dpdt(:,:), dpdt_ave(:,:), h_real_m_ave(:,:)
        real, allocatable :: P_ave(:,:), EDAx_ave(:,:), dpdx_ave(:,:), dpdy_ave(:,:), dhdt(:,:)
        real, allocatable :: U_ref(:), H_ref(:)
        ! Other
        Integer     NYs
        ! Output
        SAVE        /Current/                    ! Current timestep
        save        /CurrentT/
        save        /CurrentT_met/
        
        ! A good description of the heat equation and its origin is availible at: https://www.comsol.se/multiphysics/heat-transfer-conservation-of-energy 
        
        NN=(NYs+1)/2
        u_slip=ua-ub        ! Includes direction                                   !D_IN: /Outp/ -> Ua, Ub.
        b2Rx=B**2/Rx                                                               !D_IN: /Outp/ -> B, Rx;
        dt_real=dt*b/um                                                            !D_IN: /G0DT/ -> dt; /Outp/ -> Um
        b2RxbUm=b2Rx/dt_real
        DXSSB=DX*SS*B                                                              !D_IN: /Grid/ -> Dx
        DXSSPH=PH/(2*DXSSB)                                                        !D_IN: /Visc/ -> Ph
        temp00=temp(1,1)                                                           !D_IN: /CurrentT/ -> temp(I,J)
        
        P_lim1=1e-4
        P_lim2=1e-4
        
        umax=max(ua,ub)
        
        !Allocating the arrays based on grid parameters. Deallocation happens automatically when the subroutine is completed. 
        allocate(dPdx(1:NX,1:NN), dPdy(1:NX,1:NN), h_real_m(1:NX,1:NN), temp_i(1:NX,1:NN), umcy(1:NX,1:NN), umcx(1:NX,1:NN))
        allocate(dROdTemp(1:NX,1:NN), dpdt(1:NX,1:NN), dpdt_ave(1:NX,1:NN), h_real_m_ave(1:NX,1:NN), P_ave(1:NX,1:NN))
        allocate(EDAx_ave(1:NX,1:NN), dpdx_ave(1:NX,1:NN), dpdy_ave(1:NX,1:NN), dhdt(1:NX,1:NN), U_ref(1:NX), H_ref(1:NX))
        allocate(Cavitaion(1:NX,1:NN), cont(1:NX,1:NN))
        
        Cavitaion=Cavitaion*0 ! Matrix storing the locations of where the lubricant cavitates
        h_real_m = H*b2Rx
        
        If( t .LT. 0) call pastupd(0, NX, NN, SS)                                  !CALL: NO comment on the arguments

        ! Reset the temperature if timedependent. othervise it's not needed since past and present are the same
        if( t .GE. 0) then
            temp=tempp                                                             !D_OUT: temp -> /CurrentT/ !D_IN: /PastT/ -> Tempp
            temp_ma=temp_map                                                       !D_OUT: temp_ma -> /CurrentT_met/; !D_IN: /PastT_met/ -> temp_map
            temp_mb=temp_mbp                                                       !D_OUT: temp_mb -> /CurrentT_met/; !D_in: /PastT_met/ -> temp_mbp
        endif
        
        time=0
        fin=0
        iter=0
        
        !$OMP  PARALLEL DO                                                                 &
        !$OMP& IF(use_multiple_cores)                                                      &
        !$OMP& PRIVATE(J,J0,J2,J1,J3,JJ,I,RL_Ta,err,I0,I2,I1,I3,diffp,EpsT)                &
        !$OMP& SHARED(P_ave,temp,dPdx,dPdy,dPdx_ave,dPdy_ave,Cavitaion,P,P_lim1,           &
        !$OMP&        h_real_m_ave,h_real_m,dhdt,H,Hpast,b2RXbUm,lub_param,dROdTemp,EpsT0, &
        !$OMP&        RL_c,PH,RA1,RA2,rho_oil,DXSSPH,contact,cont,t,Ppast,dpdt,dpdt_ave,   &
        !$OMP&        dt_real,NN,NY,NYs,NX)                                                &
        !$OMP& SCHEDULE (STATIC)
        !!$OMP& SCHEDULE (GUIDED, NN/(2*cores) + 1)
        DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J2=J-2*SS
            IF( J2 .LE. 0)  J2=J0
            J1=J+1*SS
            IF( J1 .GT. NN) J1=J-SS
            J3=J1+2*SS
            IF( J3 .GT. NN) J3=NY-J3
            JJ=NYs+1-J
            IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
                err=0
                RL_Ta=temp(I,J)
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I2=I-2*SS
                    IF( I2 .LE. 0)  I2=I0
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                    I3=I+2*SS
                    IF( I3 .GE. NX) I3=NX

                    P_ave(i,j)      = (P(i0,j) + P(i,j0) + 4*P(i,j) + P(i1,j) + P(i,j1)) /8.0                    !D_IN: /CurrentP/ -> P(I,J)
                    IF( J == NN)then
                        P_ave(i,j)  = (P(i0,j) + P(i,j0) + 4*P(i,j) + P(i1,j) + P(i,j0)) /8.0
                    endif
                    

                    IF( i == NX) then
                        dPdx(i,j)       = 0
                        dPdy(i,j)       = 0
                        dPdx_ave(i,j)   = 0
                        dPdy_ave(i,j)   = 0
                        Cavitaion(i,j)  = 1  ! The fluid cavitates at the edge
                        
                    elseIF(P(i,j) .GT. P_lim1 )then     ! Do not allowe a gradient higher than the pressure itself. 
                        diffP           = (P(I1,J) - P(I0,J))
                        if( abs(diffP) .gt. P(i,j) ) diffP = sign(P(i,j),diffP) 
                        if( abs(diffP) .gt. 0.5) diffP=sign(0.5,diffP)
                        dPdx(i,j)       = diffP * DXSSPH
                        
                        diffP           = (P(I,J1) - P(I,J0))
                        if( abs(diffP) .gt. P(i,j) )  diffP = sign(P(i,j),diffP) 
                        if( abs(diffP) .gt. 0.5) diffP=sign(0.5,diffP)
                        dPdy(i,j)       = diffP * DXSSPH
                        
                        diffP           = (P(I3,J) + 2.0 * ( P(I1,J) - P(I0,J) ) - P(I2,J)) / 4.0
                        if( abs(diffP) .gt. P(i,j) )  diffP = sign(P(i,j),diffP)  
                        if( abs(diffP) .gt. 0.5) diffP=sign(0.5,diffP)
                        dPdx_ave(i,j)       = diffP * DXSSPH
                        
                        diffP           = (P(I,J3) + 2.0 * ( P(I,J1) - P(I,J0) ) - P(I,J2)) / 4.0
                        if( abs(diffP) .gt. P(i,j) )  diffP = sign(P(i,j),diffP) 
                        if( abs(diffP) .gt. 0.5) diffP=sign(0.5,diffP)
                        dPdy_ave(i,j)       = diffP * DXSSPH
                        
                        Cavitaion(i,j)  = 0  
                        
                    elseif( x(i) .GT. 0 ) then ! In the outlet region where we're likely to have cavitations. 
                        dPdx(i,j)      = 0 !(P(I1,J)-P(I0,J))*DXSSPH
                        dPdy(i,j)      = 0 !P(I,J1)-P(I,J0))*DXSSPH
                        dPdx_ave(i,j)  = 0
                        dPdy_ave(i,j)  = 0
                        Cavitaion(i,j) = 1

                    else !if(x(i) .LE. 0 .and. max(P(I0,j), P(i0,j1), P(i0,j0)) .LT. 1e-4 )then ! In the inlet region with low presure. 
                        dPdx(i,j)      = 0.0 !(P(I1,J) - P(1,J))*2*DXSSPH / ((i1-1)/SS)
                        dPdy(i,j)      = 0.0
                        dPdx_ave(i,j)  = 0.0 !dPdx(i,j)
                        dPdy_ave(i,j)  = 0.0
                        Cavitaion(i,j) = 0
                        
                    endif
                    
                    IF(j==NN) dpdy(i,j) = 0.0 ! Enforcing symmetry
                    
                     !h_real_m(i,j)      = H(i,j)*b2Rx                                                            !D_IN: /CurrentH/ -> H(I,J)
                     h_real_m_ave(i,j)  = (h_real_m(i0,j)  + h_real_m(i,j0)  + 4*h_real_m(i,j)  + h_real_m(i1,j)  + h_real_m(i,j1) ) /8.0
                     dhdt(i,j)          = (H(i,j) - Hpast(i,j))*b2RXbUm                                          !D_IN: /Past/ -> Hpast(I,J)
                     
                    if( lub_param == 4) then                                                                     !D_IN: /Ref/ -> lub_param
                        EpsT            = EpsT0*exp(-RL_c*P(I,J)*PH)                                             !D_IN: /RLarsson/ -> EpsT0, RL_c;
                        dROdTemp(I,J)   = (1+RA1*PH*P(I,J) / (1+RA2*PH*P(I,J))) * (-EpsT)*rho_oil                !D_IN: /Rho/ -> RA1, RA2; /Therm_param/ -> rho_oil
                    else if( Lub_param == 5) then
                        dROdTemp(I,J)   = 0
                    else if( lub_param == 7 .or. Lub_param == 8) then
                        dROdTemp(I,J)   = -EpsT0*rho_oil
                    else if( lub_param == 9) then
                        dROdTemp(I,J)   = -EpsT0*rho_oil
                    else if( lub_param == 12 .or. Lub_param == 13) then
                        dROdTemp(I,J)   = (1+RA1*PH*P(I,J) / (1+RA2*PH*P(I,J))) * (-EpsT0)*rho_oil 
                    else
                        WRITE(4,*)'Bad lubrication in subroutine Lubrication_temp. lub_param ='
                        WRITE(4,*) lub_param
                        stop 'Bad Lubrication'
                    endif
                    
                    ! If node or nearby nodes are in contact
                    cont(i,j)       = max(contact(i,j), contact(i0,j), contact(i1,j), contact(i,j1), contact(i,j0)) ! IF in contact or next to contact                !D_IN: /Contact_mat/ -> contact(I,J)

                    ! Time derivative of pressure for power due to compaction
                    if( t .GT. -5) then
                        diffP       = (P(i,j) - Ppast(i,j))                                                      !D_IN: /Past/ -> Ppast
                        if( abs(diffP) .gt. 0.5) diffP=sign(0.5,diffP)
                        dpdt(i,j)       = diffP *PH /dt_real
                    
                        diffP   = (P_ave(i,j) - P_ave_past(i,j))
                        if( abs(diffP) .gt. 0.5) diffP=sign(0.5,diffP)
                        dpdt_ave(i,j)   = diffP*PH /dt_real
                        IF( t .LE. -3) dpdt_ave(i,j) = dpdt_ave(i,j)*0.5 ! A slower introduction of this term which idealy shold be zero for time indep. 
                    else
                        dpdt(i,j)     = 0.0
                        dpdt_ave(i,j) = 0.0
                    endif

            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
        u_c_max     = 0            
        u_c_number  = 0
        temp_count  = 0
        temp_count_max    = 0
        ! Calculate the temperature
        ! temp at current time step = temp at old time step + temp increment ( EDA at current time step)
222     dt_lim=dt_real                                                       
        
        if (t .LE. -10) then
            dt_lim=dt_lim*C1
        endif
        
        if( Lub_param .EQ. 4) Call EDA_calc(SS, NYs, 40.0, 1.0, NN, t, 0, dpdx_ave, dpdy_ave, P_ave)                                                                        !CALL: No comment (This writes global data)
        if( dt_lim .gt. DXSSB/4.0*1.0/umax) dt_lim=DXSSB/4.0*1.0/umax                                   ! The maximum timestep is set so that the metal will not move more than half a node step 
        if(lub_param==4) Call cp_calc(NX, NN,SS)                                                        ! The material parameters depend on the pressure and temperature    !CALL: No comment (This also writes global data)
            
        !$OMP  PARALLEL DO                                                                 &
        !$OMP& IF(use_multiple_cores)                                                      &
        !$OMP& PRIVATE(J,J0,J1,I,I0,I1,umc)                                                &
        !$OMP& SHARED(NN,NX,SS,Lub_param,EDAx_ave,EDAx_a,EDAx,dpdy_ave,dpdx_ave,umcx,umcy, &
        !$OMP&        Cavitaion,h_real_m_ave,um,umax,umax_p,P,P_lim2,t,geom,U_ref,H_ref,H) &
        !$OMP& REDUCTION(max:u_c_max) REDUCTION(+:u_c_number)                              &
        !$OMP& SCHEDULE (STATIC)
        DO J=1,NN,SS
            umcy(1,J) = 0
            umcx(1,J) = 0
            temp(1,J)=temp(1,1)
            J0=J-SS
            IF( J0 .LE. 0)  J0=J
            J1=J+1*SS
            IF( J1 .GT. NN) J1=J-SS
            
            DO I=1+ss,NX,SS
            
                I0=I-1*SS
                I1=I+1*SS
                IF( I0 .LE. 0)  I0=1
                IF( I1 .GE. NX) I1=NX
                
                if( Lub_param .EQ. 4) then ! Simly incoded to 6072 to test it this makes the code more stable. If works this should be incorporated in a nicer manner. 
                    EDAx_ave(i,j)   = (EDAx_a(i0,j) + EDAx_a(i,j0)  + 4*EDAx_a(i,j)  + EDAx_a(i1,j) + EDAx_a(i,j1)) /8.0
                else  
                    EDAx_ave(i,j)   = (EDAx(i0,j) + EDAx(i,j0)  + 4*EDAx(i,j)  + EDAx(i1,j) + EDAx(i,j1)) /8.0                                          !D_IN: /Current/ -> EDAx
                endif
                
                !Enforcing symetry 
                if(J == NN ) then
                    dpdy_ave(i,j)=0  
                endif
                
                If( Cavitaion(i,j) == 1) then
                    umcx(i,j)  = um !Since umc is limited by umc in the previous node this might have an effect. 
                else
                    umc     = um-dpdx_ave(i,j)*(h_real_m_ave(i,j))**2/(12.0*EDAx_ave(i,j))
                
                    IF( abs(umc) .GT. 5*umax) then
                        umc=sign(5*umax,umc)
                    endif
                
                    IF( i .gt. 1+SS .and. abs(umc-umcx(i0,j)) .GT. umax_p*um*SS ) then                                                                  !D_IN: /C_method/ -> umax_p
                        u_c_max     = max(u_c_max,abs(umc))
                        u_c_number  = u_c_number+1
                        umc         = umcx(i0,j)+sign(umax_p*um*SS , umc-umcx(i0,j))
                   
                    endif
                    
                    !Have a small value here so that the lubricant does not stand still
                    if( P(i,j) .LT. P_lim2 .and. umc .LT. 0.2) umc=0.2  
                    
                    ! If low pressure, no backflow is allowed. 
                    !if( P(i,j) .LT. 1e-2 .and. umc .LT. umcx(i0,j) .and. i .LT. NX/2) umc=umcx(i0,j)*H(i,j)/H(i0,j) 
                    umcx(i,j)  =umc
                endif
                
                if( (t .le. 0 .and. (geom == 1 .or. geom == 4 .or. geom == 5)) .or. Cavitaion(i,j) == 1) then                                           !D_IN: /Ref/ -> geom
                    umcy(i,j) = 0       ! Since then flat surfaces with zero transverse flow
                    
                else
                    umc         =-dpdy_ave(i,j)*h_real_m_ave(i,j)**2/(12.0*EDAx_ave(i,j))
                    
                    IF( abs(umc) .GT. 5*umax) then
                        umc=sign(5*umax,umc)
                    endif
                    
                    IF( abs(umc-umcy(i0,j)) .GT. 0.6*um ) then !0.6 is based on a test of code 6076c where the time indep speed change was up to 4 m/s and um was 8.5m/s
                        umc=umcy(i0,j)+sign(0.6*um , umc-umcy(i0,j))
                        u_c_number  = u_c_number+1
                    endif

                    umcy(i,j)  =umc
                endif
                
                if(t .lt. 10 .and. j .eq. 1)then
                    U_ref(i) = umcx(i,1)
                    H_ref(i) = H(i,1)
                endif
                
                ! At low presures, mostly at the inlet, small presure canges have a tendency to cause occilations due to thier affect on the flow. 
                !if( P(i,j) .LT. P_lim1) then
                !    umcx(i,j) = U_ref(i) * min(1.5 , H(i,j) / H_ref(i))
                !    
                !    if( t .gt. 0 .and. (geom == 1 .or. geom == 4 .or. geom == 5))  then ! Just to be conservative
                !        umcy(i,j)   =   umcy(i,j) * 0.5 
                !        
                !        if(  abs(umcy(i,j)) .gt. umax) then
                !            umcy(i,j) = sign(umax , umcy(i,j)) ! Just to be more conservative
                !        endif
                !    endif
                !endif
                !if ( P(i,j) .LT. P_lim2 .and. umcx(i,j) .lt. 0) then
                !    umcx(i,j) = 0
                !endif
                !if( x(i) .LT. -1.2 .and. umcx(i,j) .lt. 0) then
                !    umcx(i,j) = 0
                !endif
                
            enddo           
        enddo
        !$OMP END PARALLEL DO
        
        iter=iter+1
            !$OMP PARALLEL DO                                                                   &
            !$OMP&  IF(use_multiple_cores)                                                      &
            !$OMP& PRIVATE(J,J0,J1,JJ,I,I0,I1,tempi0,tempi,tempi1,tempj0,tempj1,umcxi0,         &
            !$OMP&         umcxi,umcxi1,umcyj0,umcyj,umcyj1,Hi0_real,Hi1_real,Hj0_real,         &
            !$OMP&         Hj1_real,temp_incre,dT_nodes,dt_lim_c)                               &
            !$OMP& SHARED(NN,SS,NYs,NX,temp,umcx,umcy,h_real_m_ave,cont,Cavitaion,dt_real,      &
            !$OMP&        P_ave,dhdt,RO,RA1,RA2,PH,rho_oil,EDAx_ave,u_slip,dROdTemp,dpdx_ave,   & 
            !$OMP&        dpdy_ave,EpsT0,RL_c,dRodp_mat,EDAx,t,dpdt_ave,temp_i,temp00,lub_param,&
            !$OMP&        DX,B,H,P,alpha,pref,ENDA,xH,kH,RL_G0,S0,Dz,Cz,RL_T0,Tg0,YA1,YA2,YB1,  &
            !$OMP&        YB2,YC1,YC2,Yedag,Yalfag,temp_ma,temp_mb,CP_O,K_O,CP_O_past,K_met,    &
            !$OMP&        K_glass,solid_material,DCP_Dtemp,temp_fac,temp_param)                 &
            !$OMP& FIRSTPRIVATE(Z,EDA0,RL_Ta,K_oil,Cp_oil)                                      &
            !$OMP& REDUCTION(max:temp_count_max) REDUCTION(+:temp_count) REDUCTION(min:dt_lim)  &
            !!$OMP& SCHEDULE (STATIC)
            !$OMP& SCHEDULE (GUIDED, NN/(2*cores) + 1)
            DO J=1,NN,SS
                J0=J-SS
                IF( J0 .LE. 0)  J0=J+SS
                J1=J+1*SS
                IF( J1 .GE. NN) J1=NN-SS
                JJ=NYs+1-J
                IF( JJ .LE. 1)  JJ=1+SS
                    
                DO I=1+ss,NX,SS
            
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                    
                    tempi0  = temp(i0,j)
                    tempi   = temp(i,j)
                    tempi1  = temp(i1,j)
    
                    tempj0  = temp(i,j0)
                    tempj1  = temp(i,j1) 

                    umcxi0  = umcx(i0,j)
                    umcxi   = umcx(i,j)
                    umcxi1  = umcx(i1,j)
                    
                    umcyj0  = umcy(i,j0)
                    umcyj   = umcy(i,j)
                    umcyj1  = umcy(i,j1)
                    
                    Hi0_real  = h_real_m_ave(i0,j)
                    Hi1_real  = h_real_m_ave(i1,j)
    
                    Hj0_real  = h_real_m_ave(i,j0)
                    Hj1_real  = h_real_m_ave(i,j1) 
                    
                    !if( x(i) .GT. 0 ) press(i,j)=max(P(i-ss,j),P(i-2*ss,j),P(i-3*ss,j),P(i,j),P(i,j1),P(i,j0))
                    
                    if( cont(i,j)==2 .or. Cavitaion(i,j)==1 ) then 
                        temp_incre=0! If contact, the energy is added to the metals instead
                    else

                        call temp_calc_ij(i,j, SS, dt_real, P_ave(i,j), h_real_m_ave(i,j), dhdt(i,j), RO(i,j)*rho_oil, EDAx_ave(i,j), NYs, umcxi,umcxi0,umcxi1, umcyj, umcyj0, umcyj1, u_slip, tempi0, &                    
                        tempi, tempi1, tempj0, tempj1,temp00,RO(i0,j)*rho_oil ,Hi0_real,Hi1_real,Hj0_real,Hj1_real, dROdTemp(I,J),dpdx_ave(i,j),dpdy_ave(i,j),t,dpdt_ave(i,j), temp_incre) !(i,j, SS, DT, pi, Hi_real, ROi, EDAxi, NYs, umcxi, umcyi,u_slip, new_temp)   !CALL: temp_incre is written here        !D_IN: /Current/ -> RO

                        dT_nodes=max(abs(temp(i,j)-temp(i0,j)),abs(temp(i1,j)-temp(i,j)),abs(temp(i,j)-temp(i,j0)),abs(temp(i,j1)-temp(i0,j)))

                        if( abs(temp_incre) .ge. 1e12)then
                            temp_count_max = max(temp_count_max, temp_incre)
                            temp_incre = sign(1e12,temp_incre)
                            temp_count=temp_count+1
                        endif
                        
                               
                        ! Scale the time depending on how big the increment in temp is and the difference in temp between the neibuoring nodes. 
                            if( dT_nodes .GT. 40 ) then     ! The higest temp difference allowed during one step is 10 deg
                                dt_lim_c = abs(40/temp_incre*0.25)
                                if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                            elseif( dT_nodes .GT. 4 ) then  ! The normal calse
                                dt_lim_c = abs(dT_nodes/temp_incre*0.25)
                                if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                            else                            ! Add a limit of a 1 degree temp differene between metal nodes which is always allowed
                                dt_lim_c = abs(4/(temp_incre*0.25))
                                if(dt_lim_c .LT. dt_lim) dt_lim=dt_lim_c
                            endif
                    endif
                    
                    temp_i(I,J) = temp_incre
                ENDDO
            ENDDO
            !$OMP END PARALLEL DO
        
        CALL temp_calc_metal(dt_lim, SS, t, NYs, temp_convm)                                                !CALL: temp_convm is assigned here also 
        
        time=time+dt_lim
        ! If more than enough time has elapsed
        if(time .GE. dt_real)then
            dt_lim=(dt_real-(time-dt_lim)) !dt_lim2 becomes 2 times the time we have left. 
            time=dt_real
            fin=1
        endif
        
        ! Update the temperatures
        ! In the metal
        Call Met_temp_upd(dt_lim, NYs, SS, temp00, 0.0)                                                     !CALL: No comment
        NY_m=(nys+1)/2
        NX_m=(NX+1)/2
        NN_m=(NY_m+1)/2
        
        ! In the oil
        !$OMP  PARALLEL DO                                                                 &
        !$OMP& IF(use_multiple_cores)                                                      &
        !$OMP& PRIVATE(J,J0,J1,I,I0,I1,JJ,im,jm,im1,jm1,temp_ave)                          &
        !$OMP& SHARED(NN,NX,NN_m,NX_m,SS,temp_ma,temp_mb,temp,cont,Cavitaion,temp_i,       &
        !$OMP&        dt_lim)                                                              &
        !$OMP& SCHEDULE (STATIC)
        DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J1=J+1*SS
            IF( J1 .GT. NN) J1=NN-SS
            JJ=NYs-J+1
            IF( JJ .LE. 1)  JJ=1+SS
            DO I=1+ss,NX,SS
                I0=I-1*SS
                IF( I0 .LE. 0)  I0=1
                I1=I+1*SS
                IF( I1 .GE. NX) I1=NX

                ! The metal has twize the stepsize between nodes as the oil. To make it more stable in time dep simulations
                IF (MOD(i-1+ss,2*SS)==0 .AND. MOD(j-1+ss,2*SS)==0 ) THEN
                    !THE NUMBER IS AN EVEN NUMBER! For nodes of temp_metal is needed
                    im  = i/2
                    jm  = j/2
                    im1 = im+ss
                    if( im1 .GT. NX_m) Im1=NX_m-ss
                    jm1 = jm+ss
                    if( jm1 .GT. NN_m) jm1=NN_m-ss
                        
                    temp_ave   = (temp_ma(im,jm,1) +temp_ma(im1,jm,1) +temp_ma(im1,jm1,1)+ temp_ma(im,jm1,1))/8 + (temp_mb(im,jm,1) +temp_mb(im1,jm,1) +temp_mb(im1,jm1,1)+ temp_mb(im,jm1,1))/8
                    
                        
                ELSEIF(MOD(i-1+ss,2*SS)==0 )THEN
                    ! Two nodes of temp_metal is needed in x-direction
                    im  = i/2
                    jm  = (j-1)/2+1  
                    im1 = im+ss
                    if( im1 .GT. NX_m) Im1=NX_m-ss
                        
                    temp_ave    = (temp_ma(im,jm,1)+ temp_ma(im1,jm,1))/4+   (temp_mb(im,jm,1)+ temp_mb(im1,jm,1))/4
                ELSEIF(MOD(j-1+ss,2*SS)==0 )THEN
                    ! Two nodes of temp_metal is needed in y-direction
                    im  = (i-1)/2+1
                    jm  = (j)/2
                    jm1 = jm+ss
                    if( jm1 .GT. NN_m) jm1=NN_m-ss
                        
                    temp_ave    = (temp_ma(im,jm,1)+ temp_ma(im,jm1,1))/4+   (temp_mb(im,jm,1)+ temp_mb(im,jm1,1))/4
                ELSE
                    ! Only one temp_metal node is needed
                    im  = (i-1)/2+1
                    jm  = (j-1)/2+1
                    temp_ave    = temp_ma(im,jm,1)/2 + temp_mb(im,jm,1)/2
                END IF
                
                ! If in contact or in the outlet region where cavitations are formed. 
                ! Give the librication the average temperature of the two metal surfaces. 
                if( cont(i,j)==2  .or. Cavitaion(i,j)==1 ) then
                    temp(i,j)  = temp_ave
                else
                    temp(i,j)  = temp(i,j) + temp_i(i,j)*dt_lim ! Storing the updated temperatures
                endif
                
                if( temp(i,j) .LT. temp(1,1)) then
                    temp(i,j)=temp(1,1)
                elseif( temp(i,j) .GT. 1000) then
                    WRITE(4,*) 'Warning, the temperature has hit the limit at i = ', i, ' j = ', j 
                    temp(i,j)=1000
                endif
                
                if( temp(i,j) .gt. 200 .and. P(i,j) .lt. 0.5 ) then
                    temp(i,j)=200 ! Experimental test in 6078h to se if stabilizes. If works should be used to not call temp_calc_ij as ofthen. 
                    IF( temp(i,j) .LT. temp_ave ) temp(i,j) = temp_ave  !This is needed in the outlet. Othervise the oil cools of to fast here. 
                endif
                
                IF( temp(i,j) .GT. temp_ave+200 ) temp(i,j) = temp_ave+200
                
                temp(i,jj) = temp(i,j) 

                ENDDO
        ENDDO
        !$OMP END PARALLEL DO

        temp_conv=temp_conv + sum(abs(temp_i(1+SS:NX:SS,NN)))/NX * dt_lim
        !if ( t .LE. 0) WRITE(4,*)'temp_conv = ',temp_conv
        
        CALL Newtonian(SS, NYs, 40.0, 1.0, NN, t, 0, 0)                                                     !CALL: No comment
        ! Update the metal temperature in every itteration if time dependent
        if( t .GE. 0)  Call Met_temp_upd(dt_lim, NYs, SS, temp00, dt_lim)                                   !CALL: No comment
        
        if (fin==0 .and. ((t .GE. 0 .and. iter .LT. iter_lim) .or. ( iter .LT. 10) .or. (sim_end==1 .and. iter .LT. iter_cvot*iter_lim) )) goto 222              ! !!! Branch point to target 222
        sim_pass=1
        if (iter .GE. iter_lim .and. sim_end .NE. 1) then
            WRITE(4,*)'warning, not finished with the temperature timestep'
            WRITE(4,*) 'iter =',iter, ' time = ',time,' dt_real = ', dt_real
            sim_pass=0
        endif
        
        if (iter .GE. iter_cvot*iter_lim .and. sim_end == 1) then
            WRITE(4,*)'WARNING, not finished with the END temperature timestep'
            WRITE(4,*) 'iter =',iter, ' time = ',time,' dt_real = ', dt_real
            sim_pass=0
        endif

        if( t .LT. 0 ) then     !.and. ss == 1
            ! if time independent, the time steps are more restricted so this will save some computational time
            Call Met_temp_upd(dt_lim, NYs, SS, temp00, time)                                                !CALL: No comment
            M=1
            Do while( temp_convm .gt. 0.05 .and. M .LE. 20)
            CALL pastupd(2, NX, NN, SS)                                                                     !CALL: No comment
            ! Since time indep the temperature in the metal should only change is not converged
            time=time*5
            if( time .gt. DXSSB/2*1/umax) time=DXSSB/2*1/umax  ! The maximum timestep is set so that the metal will not move more than half a node step 
            CALL temp_calc_metal(time, SS, t, NYs, temp_convm)                                              !CALL: No comment (remmember temp_conv is returned here)
            Call Met_temp_upd(time, NYs, SS, temp00, 0.0)                                                   !CALL: No comment 
            Call Met_temp_upd(time, NYs, SS, temp00, time)                                                  !CALL: No comment

            ! The metal temperature diffusion is much more stable than the oil simulation
            if (t .LT. -9) then
                CALL pastupd(2, NX, NN, SS)                                                                 !CALL: No comment
                time=100*time           ! Temp_calc_met rescales the time step anyway
                if( time .gt. DXSSB/2*1/umax) time=DXSSB/2*1/umax  ! The maximum timestep is set so that the metal will not move more than half a node step 
                CALL temp_calc_metal(time, SS, t, NYs, temp_convm)                                          !CALL: No comment (remmember temp_conv is returned here)
                Call Met_temp_upd(time, NYs, SS, temp00, 0.0)                                               !CALL: No comment
                Call Met_temp_upd(time, NYs, SS, temp00, time)                                              !CALL: No comment
            endif
            M=M+1
            temp_convm = temp_convm * time
            enddo
        endif
        
        If( t .LT. 0) call pastupd(0, NX, NN, SS)                                                           !CALL: No comment
        if( u_c_number .ge. 1) WRITE(4,*) 'WARNING. u_c_number =', u_c_number, ' u_c_mac =', u_c_max
        if( temp_count .ge. 1) WRITE(4,*) 'WARNING. temp_count =', temp_count, 'temp_count_max =', temp_count_max
        
        Return
    end
    