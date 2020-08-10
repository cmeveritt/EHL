! A subroutine to update the pressure based on the increments calculated in Single grid. 
    ! This subroutine locks att all pressure increments in order to detrmine the relaxation factor C_true
    ! The further away from the solution the lower C becomes. 
    Subroutine P_update(NYs, SS, k, C11, t, k_use, M_conv, Di_max)
        implicit none
        include     'inc_Grid.h'
        include     'inc_CurrentP.h'
        include     'inc_C_method.h'
        include     'inc_Contact_mat.h'
        include     'inc_Error.h'
        include     'inc_Limit_pressure.h' 
        include     'inc_P_new_com.h'
        include     'inc_Ref.h'
        include     'inc_Temp_reduction.h'
        include     'inc_Single_step.h'
        include     'inc_Output_control.h'
        Common      /Error_conv/ Old_error, old_C_true
        integer      I, J, jj, NN, NYs,SS, k, t, k_use, M_conv
        integer     i0, i1, j0, j1
        integer     i01, i2, j01, j2
        real        ER, ABS, sum, c1, c2, C_true, C11
        real        diff, max_diff, incr, Old_error, old_C_true, C_true2
        real        test1, test2, cx_lim, Di_max
        real        dp2dx, dp2dy, dp2dx1, dp2dy1,dp2dx0, dp2dy0, dp2_max, dpdx0, dpdy0
        real        incr_limit
        save        /CurrentP/
        save        /Error_conv/
        save        /P_new_com/

        ! Similar as in Erp
        ER          = 0.0 !This error is evaluated before the pressure incements are decrease with the relaxation factor C1
        SUM         = 0.0
        max_diff    = Di_max ! Read the maximum increment as befor it was limited
        NN          = (NYs+1)/2
        
        ! A limit on the maximum pressure increment allowed. Depends on if in the initial 10 itterations of the timestep. If taking to big pressure increments straight away then the risk increases that the solution diverges
        if( M_conv .lt. 10 )then
            incr_limit = 0.001 * (M_conv+1)**2
        else
            incr_limit = 0.1
        endif

        DO 10 J=1+ss,NN, SS
            DO 10 I=1,NX, SS                                                                                                                                       !D_IN: /Grid/ -> NX
                If( P_upd(I,J) .LT. Plim .and. P_upd(I,J) .GE. 0.0 ) Then																						   !D_IN: /P_new_com/ -> P_upd; /Limit_pressure/ -> Plim		
                    diff        = ABS(P_upd(I,J)-Pold(I,J)) !Pold is updated in erp, so corresponds to the pressure used for calculating H                         !D_IN: /Error/ -> Pold
                    ER          = ER + diff
                    SUM         = SUM+P(I,J)                                                                                                                       !D_IN: /CurrentP/ -> P
                else if(P_upd(I,J) .LT. 0.0) then
                    P_upd(I,J)  = 0
                else
                    P_upd(I,J)  = Plim                      
                endif
                
    10	CONTINUE
        ER = ER/SUM
        
        IF( ER .GT. 15)then
            write(4,*) 'WARNING!!! Too high er, er= ',er
        endif

        ! Relax the pressure increments based upon the choosen C_meth defined in the input file
        If(c_meth .le. 4) then                                                                                                                                     !D_IN: /C_method/ -> c_meth
            C1=C_loc_* (1.0/k) *  1.0/max_diff                                                                                                                      !D_IN: /C_method/ -> c_loc
        
            C2=C_glob* (1.0/k) * 1.0/er                                                                                                                            !D_IN: /C_method/ -> c_glob
        else
            if(c_meth==7) then
                cx_lim = 1e-6
            else if(c_meth==8) then
                cx_lim = 1e-5
            else if(c_meth==9) then
                cx_lim = 2e-5
            else 
                cx_lim = 5e-5
            endif
                
            if(max_diff .lt. cx_lim) then
                C1 = C_max                                                                                                                                         !D_IN: /C_method/ -> c_max
            elseif(max_diff .lt. 1e-1) then
                C1 = (C_min-C_max)/(-log10(cx_lim)+log10(1e-1)) * (-log10(cx_lim) + log10(max_diff)) +C_max                                                        !D_IN: /C_method/ -> c_min
            else
                C1 = C_min
            endif
            
            if(er .lt. cx_lim) then
                C2 = C_max
            elseif( er .lt. 1e-2) then
                C2 = (C_min-C_max)/(-log10(cx_lim)+log10(1e-2))  * (-log10(cx_lim) + log10(er)) +C_max
            else
                C2 = C_min
            endif
        endif
        
        C_true=min(C1,C2)
        
        if( C_true .GT. C_max) C_true = C_max
        if( C_true .LT. C_min) C_true = C_min
        if(C_true .GT. 2*old_C_true) C_true=2*old_C_true                          ! This is initialized in Main.f90 since it uses the /Error_conv/ common block
        
        If( K .GE. 2) Then
            If( ER .LE. Old_error) then
               C_true = C_true * 0.5                                                ! Decreasing C if Er is not increasing through every itteration
            endif
        endif
        
        
        Old_error = ER
        old_C_true=C_true
        
        If( c_meth .GE. 1 .and. C_meth .LE. 9) then
            if( C_meth == 1) C_true=C11                 ! Old restriction method
            IF( (C_meth == 2 .or. C_meth == 5) .and. t .LE. -9) C_true=C11 ! Old restriction method
            IF( (C_meth == 3 .or. C_meth == 6) .and. t .LE. 0) C_true=C11  ! Old restriction method
            ! c_meth=4 and c_meth=7 means no old restrictions
            ! c_meth = 8 means higher error for when cmax is applied. 
        else
            WRITE(4,*) 'Bad C_meth input '
            stop
        endif
                    
        DO J=1+ss,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J01=J0-SS
            IF( J01 .LE. 0)  J01=J+SS
            J1=J+1*SS
            !IF( J1 .GT. NN) J1=J-SS removed since P is populated untill NY
            J2=J+2*SS
            !IF( J2 .GT. NN) J2=J-SS
            DO i=1+ss,NX-1,SS
                I0=I-1*SS
                IF( I0 .LE. 0)  I0=1
                I01=I0-1*SS
                IF( I01 .LE. 0)  I01=1
                I1=I+1*SS
                IF( I1 .GE. NX) I1=NX
                I2=I+2*SS
                IF( I2 .GE. NX) I2=NX

                incr    = C_true*(P_upd(i,j)-P(i,j))
                
                if( abs(incr) .gt. incr_limit .and. contact(i,j) == 0 .and. t .GE. 0) then                                                                         !D_IN: /Contact_mat/ -> contact
                    if (abs(incr) .gt. 0.1 .and. temp_param  .GE. 1) then                                                                                          !D_IN: /Temp_reduction/ -> temp_param
                        WRITE(4,*) 'Pressure incre limit at i,j = ', i, j
                        K    = K_use         ! To go back and update the temperature
                    endif
                    incr = sign( incr_limit , incr)
                elseif (abs(incr) .gt. 0.1) then
                    incr = sign( 0.1, incr) 
                endif
                
                P(i,j)  = P(i,j) + incr
                    
                ! Making sure that P is within the alowable range
                if(P(i,j)  .LE. 0.0)   P(i,j)  = 0.0 
                if(P(i,j)  .GT. Plim)  P(i,j)  = Plim
                
                ! The pressure should not be above twice the nerby pressures. A limit on 0.1 is added for the regions without presures. 
                ! An extra limit of maximum of 1Ph higher than nerby pressures is added which will affect the high presure region
                P(i,j) = min( P(i,j) , max( P(i0,j), P(i1,j), P(i,j0), P(i,j1), 0.01) * 1.5,  max( P(i0,j), P(i1,j), P(i,j0), P(i,j1), 0.1) + 0.5)
                
                if( P(i,j) .gt. 0.01) then ! Only do this check if the pressure has some value of significance. 
                    dp2dx = P(i1,j) - 2*P(i,j) + P(i0,j)
                    dp2dy = P(i,j1) - 2*P(i,j) + P(i,j0)
                    
                    dp2dx0 = P(i,j) - 2*P(i0,j) + P(i01,j)
                    dp2dy0 = P(i,j) - 2*P(i,j0) + P(i,j01)
                    dpdx0  = P(i,j) - P(i01,j)
                    dpdy0  = P(i,j) - P(i,j01)
                    
                    dp2dx1 = P(i2,j) - 2*P(i1,j) + P(i,j)
                    dp2dy1 = P(i,j2) - 2*P(i1,j) + P(i,j)
                    dp2_max= max( abs(dp2dx0), abs(dp2dy0), abs(dp2dx1), abs(dp2dy1) )
                    
                    if( dp2dx0 .gt. 0.3 .and. dp2dx .lt. 0.0 .and. dpdx0 .gt. 0.3) then
                        P(i,j) = ( P(i1,j) + P(i0,j) )*0.5 !This sets dp2dx to zero
                    endif
                    
                    if( -dp2dx0 .gt. 0.3 .and. -dp2dx .lt. 0.0 .and. -dpdx0 .gt. 0.3) then
                        P(i,j) = ( P(i1,j) + P(i0,j) )*0.5 !This sets dp2dx to zero
                    endif
                
                    if( dp2dy0 .gt. 0.3 .and. dp2dy .lt. 0.0 .and. dpdy0 .gt. 0.3) then
                        P(i,j) = ( P(i,j1) + P(i,j0) )*0.5 !This sets dp2dy to zero
                    endif
                    
                    if( -dp2dy0 .gt. 0.3 .and. -dp2dy .lt. 0.0 .and. -dpdy0 .gt. 0.3) then
                        P(i,j) = ( P(i,j1) + P(i,j0) )*0.5 !This sets dp2dy to zero
                    endif
                    
                    if( max( abs(dp2dx), abs(dp2dy) ) .gt. 1.0  ) then ! If the second derivative of P is to high it indicates that the pressure is diverging
                        if( P(i,j) .gt. max( P(i1,j), P(i0,j), P(i,j1), P(i,j0))) then ! Only adjust the highest pressure node
                            P(i,j) = min( abs( P(i1,j) + P(i0,j) -1.0)*0.5, abs( P(i,j1) + P(i,j0) -1.0)*0.5 ) !Adjust the pressure so that the second derivative is not higher than the allowed value
                        endif
                    endif
                endif
                
                if ( isnan(P(i,j) )) then
                    WRITE(4,*)'check 1111 P = Nan at i,j,=', i,j
                    P(i,j) = P(i-ss,j)
                endif
                
            enddo
        enddo
        
        ! Mirror the pressure around the X-axis
        DO J=1+SS,NN,SS
            JJ=NYs+1-J
            DO I=1+SS,NX,SS
                P(I,JJ)=P(I,J)
            ENDDO
        ENDDO
            
        IF(Geom .EQ. 2 .or. geom .EQ. 3 .OR. Geom .EQ. 6) THEN !Ball                                                                                               !D_IN: /Ref/ -> Geom
            CONTINUE
        ELSE    ! Cylinder
            ! Update the borderpressures
            DO I=1+SS,NX,SS
                P(I,1)=P(I,1+2*SS)
                P(I,NYs)=P(I,NYs-2*SS)
            ENDDO
        ENDIF
    
        ! Enforcing boundaryconditions
        DO J=1,NYs,SS
            P(1,J)=0.0
            P(NX,J)=0.0
        ENDDO
        
        IF(max_diff*C_true .gt. 0.1) then
            write(4,*) 'To high max_diff = ',max_diff
            k = k_use
        endif
        
        WRITE(4,*)'C1=', C1, ' C2=', C2,' C_true=', C_true, ' ER=', ER
        
    return
    end
    
    