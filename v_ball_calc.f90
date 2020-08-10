! Subrutine for calculating the velocities of Höglunds ball
	SUBROUTINE v_ball_calc(t,H00)
        implicit none 
        include     'inc_Hoglund_ball.h'
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
        include     'inc_Grid.h'
        include     'inc_G0DT.h'
        include     'inc_Gamma_ball_vec.h'
        include     'inc_Method.h'
        include     'inc_Outp.h'
        include     'inc_Visc.h'
        
        ! Input
        integer     SS, NN, t
        real        H00
        ! Calcluations
        integer     i,j, i0, j0, i1, j1, i2, j2
        real        tau(1:601,1:601), h_real(1:601,1:601)
        real        U_slip, T_tot, P_tot, Ph2DXss, dpdx, DT2, hmin
        
        save  /Hoglund_ball/
        save  /outp/
        save  /gamma_ball_vec/
        save  /G0DT/
        
        SS      = 1
    
        U_slip  = ub - ua                                                      !D_IN: /outp/ -> Ua, Ub
        NN      = (NY+1)/2                                                     !D_IN: /Grid/ -> Ny
        
        h_real  = H*b*b/Rx                                                     !D_IN: /CurrentH/ -> H; /Outp/ -> Rx, B;
        Ph2DXss = Ph /(2*DX*SS)                                                !D_IN: /Visc/ -> Ph; /Grid/ -> Dx
        ! Calculate the shear stress
        
        do j=1,NN
            do i=1,NX
                ! Only relavant with the shear in x-direction
                i2=i-2
                i0=i-1
                i1=i+1
                
                if(i .LE. 1) then
                    i0=1
                    i2=1
                elseif( i==2) then
                    i2=1
                elseif(i==NX)then
                    i1=NX
                endif
                
                    
                dpdx=(term4*P(i1,j) + term3*P(i,j) + term2*P(i0,j) +term1*P(i2,j))*Ph2DXss           !D_IN: /Method/ -> term1, term2, term3, term4
                if( i==0 .or.  i==NX .or. j==0) then
                    dpdx=0
                endif
                
                tau(i,j) = - h_real(i,j) / 2 * dpdx  + EDAx(i,j) * u_slip / h_real(i,j)              !D_IN: /Current/ -> EDAx
                
            enddo
        enddo
        
        P_tot=sum(P)*Ph*DX**2*b**2
        T_tot=sum(tau)*DX**2*b**2
        gamma_ball(t+1)=T_tot/P_tot                                                                  !D_OUT:  gamma_ball(I) -> /gamma_ball_vec/ 
        
        if( P_tot .GE. T_tot) then                                                                    
            DT2 = M_ball/ P_tot *0.5                                                                 !D_IN: /Hoglund_ball/ -> M_ball
        else
            DT2 = M_ball/ T_tot *0.5
        endif
        
        j=NN
        i=(NX+1)/2
        hmin = H_real(I,J)
        if( hmin .LT. 0) then
            WRITE(4,*)'Warning bad Hmin. Hmin= ', hmin 
            hmin=0
        endif
        IF( DT2 .GT. DT ) then                                                                       !D_IN: /G0DT/ -> DT
                 ! Only look for Hmin at used gridpoint for current level. 

            IF( DT2 .GT. abs(hmin / Vv_ball * 0.05))  then                                           !D_IN: /Hoglnd_ball/ -> Vv_ball
                DT2 = abs(hmin / Vv_ball * 0.05)
            endif
            
            H00 = H00 - (Vv_ball * DT2)*Rx/(B*b)
        ELSE
            DT2=DT0     ! Original limit for DT
            H00 = H00 - (Vv_ball * DT)*Rx/(B*b)
        ENDIF
        DT=DT2          !Rescaling DT for the rest of the simulation
        
        Vv_ball     = Vv_ball - P_tot / M_ball * DT                                                  !D_OUT: Vv_ball -> /Hoglund_ball/
        Vh_ball     = Vh_ball - T_tot / M_ball * DT                                                  !D_OUT: Vh_ball -> /Hoglund_ball/
        omega_ball  = omega_ball + T_tot * Rx / I_ball *DT                                           !D_IN: /Hoglund_ball/ -> I_ball   !D_OUT: omega_ball -> /Hoglund_ball/
        
        Ub = Vh_ball            ! Bottom surfce moves in positive x-direction                        !D_OUT: Ub -> /outp/
        ua = omega_ball * Rx    ! Top surface                                                        !D_OUT: Ua -> /outp/
        
         WRITE(52,*) Vv_ball,    Vh_ball,    omega_ball,  gamma_ball(t+1)
         WRITE(52,*) ua,  ub, P_tot, T_tot
         WRITE(52,*) H00, DT, hmin
         
         If ( Vv_ball .LT. 0 .and. H00 .GT. 1.5)then
                     WRITE(4,*)'Finished :) Enough height reached ' 
                    stop 'Finished :) Enough height reached ' 
         endif
         
    return
    end
    
                