! Subroutine for uppdating the pressure outside the simulated domain. This is only used if simulating a cylinder. 
    ! The elastic deformation is dependent on the pressure outside of the simulated area. Therefore is the pressure extrapolated outside the simulated area in this subroutine. 
	! In the normal cylinder case the surface iregularities do not affect the pressure outside the simulated domain. Thus is this pressure not updated for the regions where the pressure could be afected by the time dependent parts
    ! However, some simulations looks inte line defects and thus has the side pressures to be updated for the whole area. This is added as the extra ittrations at the end. 
    SUBROUTINE Pressure_outside(t,SS, NYs, k, k_use)
        implicit none 
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_Grid.h'
        include     'inc_Itt.h'
        include     'inc_Outp.h'
        include     'inc_Ref.h'
        include     'inc_Pres_ave_param.h'
        include     'inc_P_line_side.h'
        ! Input
        integer     NYs, SS, k, k_use, M_conv
        Integer     t
        ! Calculations
        integer     NN, I, J, JJ, number_n
        real        av_p, av_T,  SRR
        integer     position
        ! Output
        save        /P_line_side/
        save        /CurrentT/

        NN=(NYs+1)/2    ! NYs is not the same as NY in all cases. NYs is redced for the initial time independent solutions in some cases
        
        IF (t .LE. -10  )then                   ! if time indep (Asp_shape=6 is a line asperity  
            position=2                          ! The position parameter is here used only for defining which areas could have been affected by the time dependent solution and which are not. 
                                                ! The position is the node number which starts at 1 in the inlet and ends at NX at the outlet
        else if( (t .LE. ntime/2 .and. p_ave_param ==1 ) .and. ( mod(k,3) == 1 ) .and. (k .ne. k_use) ) then    ! Update the later part of the pressure which is not yet affected by the asperity. 7074b shows that not good to smooth at last k iteration  !D_IN: /itt/ -> ntime; /pres_ave_param/ -> p_ave_param
            SRR=2*(ua-ub)/um                                                                                    ! D_IN: /outp/ -> Ua, Ub, Um
            position=max(1, floor(t*(1+abs(SRR))*f)) + floor(NX/10.0)                                           ! Assuming that the effect of the asperity travels with the fastest surface + a safty distance of 10%  !D_IN: /Grid/ -> Nx; /Ref/ -> F                
            WRITE(4,*)'P_smoothening used from position ', position
            IF( position .gt. NX) position = NX                                                                 
        else
            position = NX   
        endif   
            
        if( position .LT. NX) then
            Do I=position,NX,SS ! !!! Rowwise access to arrays
                number_n=0
                av_P=0
                av_T=0
                DO J=1,NN,SS
                    number_n=number_n+1
                    av_P=av_P+P(I,J)                !D_IN: /CurrentP/ -> P
                    av_T=av_T+Temp(I,J)             !D_IN: /CurrentT/ -> Temp
                enddo
                av_P=av_P/number_n
                av_T=av_T/number_n
                    
                DO J=1,NN,SS  
                    Temp(i,j)=av_T
                enddo
                    
                IF( t .LT. -10) then
                    P_line(I)   = av_P !* 0.75 + P_line(I) * 0.25                     !D_OUT: P_line -> /P_line_side/
                else
                    P_line(I)   = av_P !* 0.5 + P_line(I) * 0.5
                endif
                    
            enddo
        endif
        
        ! Update the initial parts of the side pressure as well in the following cases. Since these cases have geometries affecting the pressure outside of the simulated domain
        IF( asp_shape   == 1  .or. &     ! The folowing asp_shapes has line asperities that reaches outside of the analysed domain. Thus affecting the pressure outside the simulated domain 
            asp_shape   == 3  .or. &
            asp_shape   == 6  .or. &
            asp_shape   == 7  .or. &
            asp_shape   == 201.or. &
            asp_shape   == 203  )Then  
            Do I=1,position,ss                ! !!! Rowwise access
                number_n=0
                av_P=0
                       
                DO J=1,NN,SS
                    number_n=number_n+1
                    av_P=av_P+P(I,J)               
                enddo
                av_P=av_P/number_n
                P_line(I)   = av_P !* 0.5 + P_line(I) * 0.5
                        
            enddo
        endif  
            

    RETURN
    END
