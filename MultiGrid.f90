! Subroutine for the multigrid setup
    ! The subroutine defines which level to start at and controles which level we are currently at
    Subroutine MultiGrid(t, kvot, DT2, H00, iteration)
    implicit none
    include     'inc_Grid.h'
    include     'inc_Ref.h'
    include     'inc_Yasutomi.h'
    include     'inc_Itt.h'
    ! Input
    integer     t,      iteration
    real        kvot,   DT2,        H00
    ! Calculations
    integer     SS,     step(4),    level,      Sstart
    integer     NYs
    real        Tcvot, T40

    ! This subroutine does not produce any internal output. Instead the results are printed out in the output subroutine
    
    ! Defining the different grid step sizes for each grid level
    step(1) = 8 
    step(2) = 4
    step(3) = 2
    step(4) = 1 
    
    T40     = 40                                    ! Defining a reference temperature for which a solution for the Non Newtonian lubrication is found easy enough
    Tcvot   = 0                                     ! The cvot between using the temperature of the problem and theT40 temperature
    IF(Ta .LE. T40) Tcvot=1.0                               !D_IN: /Yasutomi/ -> Ta
    IF(lub_param .NE. 5 .and. lub_param .NE. 58) Tcvot=1.0     !D_IN: /Ref/ -> lub_param
    
    NYs     = 89                                    ! Unnessesery to calculate with more nodes in y-direction at the start if cylinder
    IF(NY .LT. NYs) NYs = NY                        !D_IN: /Grid/ -> Ny
    If(Geom .EQ. 2 .OR. Geom .EQ. 3 .OR. Geom .EQ. 6) NYs = NY                        ! Since then its a ball !D_IN: /Ref/ -> Geom
    
    IF(t .LT. -10)THEN

        IF(Multi_grid_param == 1) Then ! Use multi grid to find h00   !D_IN: /Ref/ -> Multi_grid_param
            IF(NX .GT. 200 .AND. NY .GT. 160)THEN                     !D_IN: /Grid/ -> Nx
                Sstart  = 1                             ! Start level of simulation
                NYs     = 161                           ! Since using 4 levels instead of three
            
                IF(NY .LT. NYs) NYs = NY
                If(Geom .EQ. 2 .OR. Geom .EQ. 3 .OR. Geom .EQ. 6) NYs = NY                ! Since then its a ball
            ELSE
                Sstart  = 2                             ! The coarsest solution should not be too coarse. 
            ENDIF
        ELSE
            Sstart=4
        ENDIF
        
            
        DO level=Sstart,4    ! !!! this is where the 4 or 3 -11's come from 
            SS      = step(level)
            CALL Grid_p(SS)                         ! SS=Step size at current levenl   !CALL: SS
            
             WRITE(4,*)'Going to level =', level     ! Write it uot both to file and comand window
             WRITE(4,*)
                                                                            ! Update gridparameters
            CAll HREE(H00,t,SS, NYs, T40, Tcvot, 0, 1, 0)                   ! Update lubrication params !CALL: H00, t, SS, NYs, T40, Tcvot, k = 0, k_use = 1, M_conv = 0
            CALL Conv(t, kvot, DT2, H00, iteration, SS, NYs, T40, Tcvot)    ! Solve the problem         !CALL: t, kvot, DT2, H00, iteration, SS, NYs, T40, Tcvot
                        
            ! When converged or ran out of tries. Output whatever we have 
            IF(level .LT. 4)THEN
                CALL OUTPUT                         ! For lvl 4 Output is called at end of this file. 
                SS  = step(level+1)
                CALL Grid_p(SS)                                        !CALL: No comment
                CALL Fine(H00,SS, kvot, DT2, t, NYs, T40, Tcvot)       !CALL: No comment
            ENDIF

        ENDDO
        
        ! Increase the temperature slowly if Non Newtonian at high temperature
        IF(Lub_param .EQ. 5 .AND. Ta .GT. T40)THEN        
            DO WHILE ((Ta-T40)*tcvot .LT. Ta-T40)
                IF( (Ta-T40)*tcvot .LT. 20) THEN
                    tcvot   = tcvot+4/(Ta-T40)      ! Increasing 4 degreas at a time
                    
                elseif ( (Ta-T40)*tcvot .LT. 40) THEN
                    tcvot   = tcvot+2/(Ta-T40)      ! Increasing 2 degreas at a time
                    
                else
                    tcvot   = tcvot+1/(Ta-T40)      ! Increasing 1 degrea at a time
                    
                ENDIF
                
                WRITE(4,*)'Going to tcvot =', tcvot                             
                WRITE(4,*)
                
                ! Update material parameters
                CALL pastupd(0,NX, (NYS+1)/2, SS)                                                                           !CALL: No comment
                CAll HREE(H00,t,SS, NYs, T40, Tcvot, 0, 1, 0)                            ! Update lubrication params        !CALL: No comment
                CALL Conv(t, kvot, DT2, H00, iteration, SS, NYs, T40, Tcvot)                                                !CALL: No comment
                CALL OUTPUT                                                                                                 !CALL: No comment
            enddo
            
                tcvot=1.0
                WRITE(4,*)'Going to tcvot =', tcvot                         
                WRITE(4,*)
                CAll HREE(H00,t,SS, NYs, T40, Tcvot, 0, 1, 0)                            ! Update lubrication params        !CALL: No comment
                CALL Conv(t, kvot, DT2, H00, iteration, SS, NYs, T40, Tcvot)                                                !CALL: No comment
                CALL OUTPUT                                                                                                 !CALL: No comment
        ENDIF
        
        ! If did not solve for all nodes in Y-direction. Add solutions to these nodes
        IF (NYs .NE. NY) CALL NYinterp(H00, t, NYs, T40, Tcvot)  
        

    ! When time dependent. Solves only on finest grid. Had trouble using courser gids here due to time derivatives and different grids might yield different H00
    ELSE                        
        
        Tcvot = 1.0                                         ! Temperature cvote
        
            level = 4                                       ! !!! Level 4 for all time steps that aren't -11
            SS    = step(level)                             ! SS=Step size at current levenl
            WRITE(4,*)'Going to level =', level             ! Write it out both to file and command window
            WRITE(4,*)

            CALL Grid_p(SS)                                 ! Update gridparameters    !CALL: No comment
            CAll HREE(H00,t,SS, NY, T40, Tcvot, 0, 1, 0)             ! Update lubrication params  !CALL: No comment
            
            ! Converged or ran out of tries. Output whatever we have 
            CALL Conv(t, kvot, DT2, H00, iteration, SS, NY, T40, Tcvot)                           !CALL: No comment
                                              
    ENDIF
    
    if( Ntime .LT. 800) then                                                                      !D_IN: /Itt/ -> Ntime
        CALL OUTPUT                                             ! Converged or ran out of tries. Output whatever we have    !CALL: No comment
    elseif( Ntime .lt. 1600 .and. MOD(t,2)==0) then
        CALL OUTPUT                                                                                                         !CALL: No comment
    endif
    
    
    Return
END
    