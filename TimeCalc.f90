! Subrutine keeping track of the time step. This subroutine also ensures taht everyting is set-up in the begining
    SUBROUTINE TimeCalc(H00, aspw, executionTime) !D_IN: H00, aspw as arguments 
#ifdef _OPENMP        
        Use omp_lib
#endif
        implicit none 
        include     'inc_CurrentT.h'
        include     'inc_G0DT.h'
        include     'inc_Grid.h'
        include     'inc_Hoglund_ball.h'
        include     'inc_Itt.h'
        include     'inc_Outp.h'
        include     'inc_PastT.h'
        include     'inc_Ref.h'
        include     'inc_Temp_param.h'
        include     'inc_Y_liu.h'
        include     'inc_Single_step.h'
        integer     t,      Mm,     M,      iteration
        real        H00,    DDw,    DTT
        real        kvot,   DT2
        real        aspw 
        integer     i, Ntime_s
        save        /PastT/
        real        stepTime
        real        executionTime
        Iteration = 1
        
        ! Starting calculations	
        ! Calculate the deformation matrix AK
        CALL SUBAK  !CALL no data passed
        
        ! Setup grid and start guess. 
        CALL INITI(H00)  !CALL D_OUT: H00
        
        !-----------------------------------------------------------------------------------
        !-----------------Calculate a single step only if requested-------------------------
        !-----------------------------------------------------------------------------------
        if(single_step_only) then
            CALL Get_single_step_data(t,H00) !Get current t from here
#ifdef _OPENMP
            executionTime = OMP_GET_WTIME()              !Start measuring the execution time. Not needed if one isn't interested in the performance of the code.
#endif
            IF( tmeth .EQ. 2) kvot=0.5                                                                                                 
            IF( tmeth .EQ. 4) kvot=0.5
            DT2     = 1.0/(2*DT)
        
            IF( lub_temp .NE. 0) Then                                                                                                 
                tempp=temp                                           
            endif
        
            IF( asp_shape .eq. 0) GOTO 333                                                                                             

        
            Ntime_s=Ntime                                           
            if( asp_shape .GT. 30 .and. asp_shape .LT. 100)     Ntime_s = Ntime * (1+5*aspw/(XE-X0))
            if( asp_shape .ge. 200 .and. asp_shape .le. 209)    Ntime_s = Ntime + (135-15)*Ntime/Nx
             
            IF( Geom .EQ. 6) Then  ! Höglund 1999 Influence of lubricant properties on elastohydrodynamic lubrication
                WRITE(4,*)'Warning! Geom = 6 is for simulating the experiments by Höglund. This is not fully developed yet'  
                Call v_ball_calc(t,H00)                                                                                                                                                               
            endif
            
            if (Ntime .GT. 100 .and. t == 10) CAll P0save ! Update the time indep pressure                                                                                                  
            
            IF(asp_shape .LE. 30 .or. asp_shape .GE. 100)iteration=iteration+1
        
            ! Store solution from previous timestep
            CALL Pastupd(1,NX, (NY+1)/2, 1)                                                                                                                                                   
            CALL initial_guess                                                                                                                                                             

            WRITE(4,*)'Going to, t=',t 
            
            ! Solve the current timestep
            Call MultiGrid(t, kvot, DT2, H00, iteration)
        
            print *, 'Completed single step'
            Goto 333
        endif
        ! ----------------------------------------------------------------------------------
        ! --------------Calculate static soloution------------------------------------------
        !-----------------------------------------------------------------------------------
#ifdef _OPENMP
            executionTime = OMP_GET_WTIME()            !Start measuring execution time. Not needed if one isn't interested in the performance of the code.
#endif
        t       = -11                                                    ! Controles asperity position
        kvot    = 1.0                                                   ! How much weigth is put on the current time step. For static, kvot =1, for timedep acc C-N kvot=0.5      
        DT2     = 0.0                                                   ! =1/(2.*DT). For timeindependent DT=inf so DT2=0 
        
        ! WRITE(*,*)'Starting simulation, t=',t
        WRITE(4,*)'Starting simulation, t=',t
#ifdef _OPENMP
        stepTime = OMP_GET_WTIME()
#endif        
        Call MultiGrid(t, kvot, DT2, H00, iteration)                !CALL: t, kvot (weight of step), DT2 (=1/2*dt), H00, iteration
#ifdef _OPENMP
        stepTime = OMP_GET_WTIME() - stepTime
        Write(70,*) t , ' : ', stepTime
#endif        
        ! ----------------------------------------------------------------------------------
        ! --------------Gently add the timedependence ------------------------------------------
        !-----------------------------------------------------------------------------------
        kvot    = 1.0
        DO iteration = 1,10                                             ! Add the timedependance over a certain number of steps. If changing here, remember to switch length oc success
#ifdef _OPENMP
            stepTime = OMP_GET_WTIME()
#endif  
            DTT = 1000.0/(iteration**3)*DT                              ! Gradually increasing the timeaffects. 
            DT2 = 1./(2*DTT)
            t=-11+iteration                                             ! !!! and this is where the -10 to -1 steps come from
            CALL Pastupd(0,NX, (NY+1)/2, 1)                                                ! Update parameters of past timestep    !CALL: No comment
            ! WRITE(*,*)'Going to, t=',t ,'iteration =', iteration
            WRITE(4,*)'Going to, t=',t ,'iteration =', iteration
            Call MultiGrid(t, kvot, DT2, H00, iteration+1)                                                                         !CALL: No comment
#ifdef _OPENMP
        stepTime = OMP_GET_WTIME() - stepTime
        Write(70,*) t , ' : ', stepTime
#endif 
        ENDDO 
        
        iteration=iteration+1
        
        ! ----------------------------------------------------------------------------------
        ! --------------Calculate time dependent  sulution---------------------------------
        !-----------------------------------------------------------------------------------
        IF( tmeth .EQ. 2) kvot=0.5                                                                                                 !D_IN: /ref/ -> tmeth
        IF( tmeth .EQ. 4) kvot=0.5
        DT2     = 1.0/(2*DT)
        
        IF( lub_temp .NE. 0) Then                                                                                                  !D_IN: /temp_param/ -> lub_temp
            tempp=temp                                      !Update the past temperature                                           !D_OUT: tempp -> /PastT/  
        endif
        
        ! Skip timedepending if no asperity to save time. Since non time dep reference case
        IF( asp_shape .eq. 0) GOTO 333                                                                                             !D_IN: /Ref/ -> asp_shape     ! !!! Branching point to target 333 at the end of the file
  
        CAll P0save                                                     ! Saves the time independent pressure                       !CALL:  No comment
        CALL Pastupd(0,NX, (NY+1)/2, 1)                                                                                             !CALL: No comment
        
        ! For all timesteps, calculate the timedepending solution
        Ntime_s=Ntime                                                   ! The actualnumber of time steps taken. I'm increasing this instead of changing earlier in the code for simplicity.   !D_IN: /Itt/ -> Ntime
        if( asp_shape .GT. 30 .and. asp_shape .LT. 100)     Ntime_s = Ntime * (1+5*aspw/(XE-X0))
        if( asp_shape .ge. 200 .and. asp_shape .le. 209)    Ntime_s = Ntime + (135-15)*Ntime/Nx          ! Ntime/Nx=1/F         ! add (135-15) movements of the pressure since it's 135 is the width of the rough surface and this will give it time to move almost out of the simulated area.
        
        DO t = 0,Ntime_s                               !This is where the rest of the time dependent solution comes from
#ifdef _OPENMP
            stepTime = OMP_GET_WTIME()
#endif 
            IF( Geom .EQ. 6) Then                                                                                                                                                             !D_IN: /Ref/ -> Geom
                Call v_ball_calc(t,H00)                                                                                                                                                       !CALL: No comment
                    
            endif
            
            if (Ntime .GT. 100 .and. t == 10) CAll P0save ! Update the time indep pressure                                                                                                    !CALL: No comment
            
            IF( asp_shape .LE. 30 .or. asp_shape .GE. 100)iteration=iteration+1
        
            ! Store solution from previous timestep
            CALL Pastupd(1,NX, (NY+1)/2, 1)                                                                                                                                                   !CALL: No comment
            CALL initial_guess                                                                                                                                                                !CALL: No comment
            
            ! WRITE(*,*)'Going to, t=',t 
            WRITE(4,*)'Going to, t=',t 
            
            ! Solve the current timestep
            Call MultiGrid(t, kvot, DT2, H00, iteration)                                                                                                                                      !CALL: No comment
            

#ifdef _OPENMP
        stepTime = OMP_GET_WTIME() - stepTime
        Write(70,*) t , ' : ', stepTime
#endif       
        ENDDO
               
333 RETURN                                                                                                                                                                                    ! !!! Branch target 333
    END
! **********************************************************!************************************************************************************************************************************************