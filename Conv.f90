! This subroutine generates a converged solution with help of multigrid methods
    ! The subroutine varies the relaxation factor C1 based on the difference between the current and the past solution. 
    Subroutine Conv(t, kvot, DT2, H00, iteration, SS, NYs, T40, Tcvot)
    implicit none
	include     'inc_Grid.h'
	include     'inc_Itt.h'
 	include     'inc_Numbers_of_bad.h'
	include     'inc_Ref.h'
    include     'inc_Shear_nbr.h'
    include     'inc_Success_vec.h'
	include     'inc_Temp_param.h'
	include     'inc_Visc.h'
    include     'inc_Yasutomi.h'
    include     'inc_Output_control.h'
    include     'inc_Single_step.h'
    
    ! Input
    integer     NYS, NN
    ! Calculation param
    integer     MK, MMmin, MM, M_conv, SS, lim, sim_pass, iter_lim, iter_cvot, ntime2
    real        H00ink, DWold, DW, DDW, ER
    real        C1, ER_ss, nr_DW, temp_conv, temp_convm
    ! Other
    integer     t, iteration, M_init
    real        kvot, DT2, H00,  T40, Tcvot, max_diff 
    ! Output
    save        /Numbers_of_bad/
    save        /shear_nbr/
    save        /success_vec/
    save        /outcontrol/
    
    ! Initiate counter for successfull iterations. Assume all iterations will be successfull   
    if( t == -11 .or. single_step_only) then 
        if( ntime .lt. 800) then              !D_IN: /itt/ -> ntime
            success(1:Ntime+15)=1.0           !D_OUT: success(I) -> /success_vec/
        elseif( ntime .lt. 1600) then
            ntime2=Ntime/2
            success(1:ntime2+15)=1.0
        endif
    endif
    
    ! Specify convergence limits dependent on if in the initial timeindependent regime or in the full time dependent regime. 
    IF(t .LE. 0) THEN
        MK=MK_stat                       !D_IN: /itt/ -> MK_stat
        if(SS .GT. 1) MK=MK*2
        ER_ss=ER_stat!/SS                !D_IN: /itt/ -> ER_stat
        ! A minimum amount of itterations are defined to ensure better stabillity. Set to 12 since the 9 first has a tighter limit on the pressure increment
        MMmin=12                                                       
    ELSE
        MK=MK_time
        ER_ss=ER_time!/SS
        ! Mimimum amount of itterations. This stablizes the code. ! When H00 is locked, H00 occilations cannot occure and MMmin can be lower.
        MMmin=12                                                       
    ENDIF
        

    NN=(NYs+1)/2        ! Node number at symmetry line
    M_conv=0            ! Start counter for whole itteration loop
    MM=0                ! Start counter for internal loop for convergence
    DWold=0.0           ! Old increment 
    ER=0.1              ! Start at C1=0.15
    iter_cvot=1         
    M_init=1
    sim_pass=1          ! initialize sim_pass
    
    ! The relaxation parameter C1 is scaled based on how close we are to the right solution.
10  IF(Tcvot .LT. 1 .and. Tcvot .NE. 0 .AND. ER .GT. 5E-4) THEN     
        C1=0.02
        iter_lim=2000
        WRITE(4,*)'Tcvot =', Tcvot, 'So C1 and C2 are =', C1
        
    ! To low pressures can cause problems with convergence, and also lub_param = 5 which is Non Newtonian has a harder time converging
    ELSEIF (    Mm .GT. 9 .AND. ER .LT. 0.00001 .and. Ph .GT. 2E9   .AND. Lub_param .NE. 5)THEN      !D_IN: /Ref/ -> lub_param; /Visc/ -> Ph
        C1=0.5
        iter_lim=1200
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(( Mm .GT. 6 .AND. ER .LT. 0.0001  .and. Ph .GT. 2E9   .AND. Lub_param .NE. 5 .and. lub_temp .NE. 1)  .or. (ER .LT. 0.00002 .and. Mm .GT. 2))THEN    !D_IN: /Temp_param/ -> lub_temp
        C1=0.4
        iter_lim=1000
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(( Mm .GT. 3 .AND. ER .LT. 0.005   .and. Ph .GT. 1.8E9 .AND. Lub_param .NE. 5 .and. lub_temp .NE. 1)  .or. (ER .LT. 0.00015 .and. Mm .GT. 2)      .or. (ER .LT. 0.0002 .and. Mm .GT. 40))THEN
        C1=0.3
        iter_lim=800
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(( Mm .GT. 3 .AND. ER .LT. 0.005   .and. Ph .GT. 1.8E9 .AND. Lub_param .NE. 5 .and. lub_temp .NE. 1)  .or. (ER .LT. 0.00025 .and. Mm .GT. 2)       .or. (ER .LT. 0.0004 .and. Mm .GT. 30))THEN
        C1=0.22
        iter_lim=600
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(( Mm .GT. 1 .AND. ER .LT. 0.1     .and. Ph .GT. 1.8E9 .AND. (Ta .LE. T40 .OR. tcvot .EQ. 0))         .or. (ER .LT. 0.0004 .and. Mm .GT. 2)       .or. (ER .LT. 0.0008 .and. Mm .GT. 20))THEN
        C1=0.15
        iter_lim=400
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(( Mm .GT. 1 .AND. ER .LT. 0.1     .and. Ph .GT. 1.8E9 .AND. (Ta .LE. T40 .OR. tcvot .EQ. 0))         .or. (ER .LT. 0.0006 .and. Mm .GT. 2)       .or. (ER .LT. 0.001 .and. Mm .GT. 10))THEN
        C1=0.10
        iter_lim=200
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(  (lub_param .NE. 5 .and. Ph .GT. 1.8E9)                                                             .OR. (ER .LT. 0.001 .and. Mm .GT. 1))THEN
        C1=0.07
        iter_lim=100
        WRITE(4,*)'C1 and C2 are ', C1
    ELSEIF(  (lub_param .NE. 5 .and. Ph .GT. 1.8E9)                                                             .OR. (ER .LT. 0.01 .and. Mm .GT. 1))THEN
        C1=0.04
        iter_lim=50
        WRITE(4,*)'C1 and C2 are ', C1
    ELSE
        C1=0.02
        iter_lim=50
        WRITE(4,*)'C1 and C2 are ', C1
    ENDIF
        
    temp_conv = 0
    
    ! Ensure that the pressure is not lower att he central line if analytic defect. This should be updated since analytical defects can be indents.     
    IF( (M_conv == 25 .or. M_conv == 50) .and. t .LT. 0.8*Ntime .and. ( asp_shape .lt. 200 .or. asp_shape .gt. 209) ) CALL P_centre_fix       !D_IN: /Ref/ -> asp_shape
    
    ! Call SingleGrid to itterate the new pressure profile. This is done KK times. 
    CALL SingleGrid(t,H00,DW,kvot,DT2, SS, NYs, C1, T40, Tcvot, lim, M_conv)                    

    ! Call ERP in order to evaluate the difference between the new and the old pressure. 
    CALL ERP(ER,Nys,max_diff, SS)                                                   ! CALL: max_diff is the output and the commonblock Error is  updated
    
    ! Write the convergence status
    WRITE(4,*)' t=', t,' M_conv= ',M_conv,' SS= ', SS
    WRITE(4,*)' ER= ', ER, 'DW', DW,' H00=', H00                                     
    WRITE(4,*)' The maximum pressure diff in a node is =', max_diff  
    IF( bad_h_nr_line     .GT. 0 .OR. bad_h_itt_line  .GT. 2) write(4,*) ' Bad h nr line= ', bad_h_nr_line,  'bad h itt line= ', bad_h_itt_line                                             !D_IN: /Numbers_of_bad/ -> bad_h_nr_line, bad_h_itt_line;
    
    !If we have contact in nodes. These values are only reset if we have contact
    IF( Max_contact .gt. 0.0) THEN                                                           
        IF( Max_contact .GT. 0.0 )                            write(4,*) ' Max_contact= ',   Max_contact                                                                                    !D_IN: /Numbers_of_bad/ -> Max_contact
        IF( bad_h_nr_node .GT. 0 .OR. bad_h_itt_node  .GT. 2) write(4,*) ' Bad h nr node= ', bad_h_nr_node, 'bad h itt node= ', bad_h_itt_node
    ENDIF
    
    IF( bad_p_node      .GT. 0)                               write(4,*) ' Bad p node= ',    bad_p_node                                                                                     !D_IN: /Numbers_of_bad/ -> bad_p_node
    IF( lim             .GT. 0)                               write(4,*) ' Number of pressure limits= ', lim                          
    IF( bad_h_end_node  .GT. 0)                               write(4,*) ' Bad lubricationg heights after contact routine = ', bad_h_end_node                                               !D_IN: /Numbers_of_bad/ -> bad_h_end_node
    IF( shear_iter      .GT. 0)                               write(4,*) ' Nonconverged shear thinning at ', shear_iter, ' times. The largest remaining increment was ', shear_incre        !D_IN: /shear_nbr/ -> shear_incre
    IF( temp_iter_t     .GT. 0)                               write(4,*) ' Nonconverged temperature at ', temp_iter_t, ' times'                                                             !D_IN: /shear_nbr/ -> temp_iter_t
    WRITE(4,*) 
    
    ! If the full data is requested, the writ H00 if it is not done before. 
    if ((t .EQ. 0) .and. collect_full_data .and. (.not. H00_collected)) then
        write (60,*) H00
        H00_collected = .true.
    ENDIF
    
    ! Reset contact and temperature counters
    bad_p_node      =0
    bad_h_nr_line   =0
    bad_h_itt_line  =0
    Max_contact     =0.0
    shear_incre     =0
    shear_iter      =0
    temp_iter_t     =0
    
    ! Adjust H00 based on the pressure levels. 
    IF( t .LE. -11 .or. (t .LT. -1 .and. Dw_meth == 1))Then
        DDW=abs((DW-DWold)/(DW))
        DWold=DW
        nr_DW=2
        if( shear_thin .GT. 0) nr_DW=5                                                                                                                      !D_IN: /Ref/ -> shear_thin
        
        ! If the Pressure is good enough or we done enough internal itterations, then uppdate the offset H00 
        IF(DDW .LT. 0.01 .AND. (ER .LT. 0.001 .OR. Mm .GT. 20) .and. (Mm .GT. nr_DW/C1 .or. DW .GT. 0.05)) THEN    
            WRITE(4,*)'Updating H00. At Mm =',Mm 
            
            ! If normal TEHL simulations
            IF( Geom .LE. 5) then                                                                                                                           !D_IN: /Ref/ -> Geom
                CALL H00upd(H00, t, DW, SS, NYs, H00ink)                    ! Uppdate the pressure based on the loadbalance                                 !CALL:	H00 is changed in this function. Also H00ink is an output as well.
                IF( H00ink .GT. 1E-6)THEN                                   ! The increment has to be big enough to count
                    CALL HREE(H00,t,SS, NYs, T40, Tcvot, 0, 0, M_conv)      ! Update lubrication parameters                                               
                    CALL pastupd(0, NX, NN, SS)                             ! Update the parameters of the past time to get rid of time dependent effects.    
                    Mm=0
                ENDIF
                
            ! Geom = 6 was set-up in a try to replicate the experiments by Höglund 2006 of a ball shot against a lubricated surface.    
            ELSE 
                DW = 0      
            ENDIF    
        ENDIF
        
        ! Reset the pressure and the other params to Hertzian pressure since not converged enough    
        IF(ER .GT. 0.01 .AND. MM .GT. 30*M_init) THEN                          
            WRITE(4,*)'Resetting the pressure. At Mm =',Mm
            CALL H00upd(H00, t, DW, SS, NYs, H00ink)                  
            
            ! Uppdate the offset with caution since not a converged solution. This ensures taht something has change from previous itterations 
            H00ink = H00ink * 0.1  
            
            ! Reset the pressure
            CALL INITI(H00)                                                                              
            
            ! Update the aprameter of the previous time step to get rid of any time dep effects
            CALL pastupd(0, NX, NN, SS)                                                                                                                                                 
            Mm=0                                                ! Reset internal counter
            M_init=M_init+1                                    
        ENDIF
    ENDIF
    
    ! Update the temperature fields.  
    ! sim_pass is an output parameter staing if the temperature simulations succeded or not. 
    if( lub_temp == 1 .and. SS .LE. 2) CALL Lubrication_temp(SS, NYs, C1, t, 0, sim_pass, iter_lim, iter_cvot, temp_conv, temp_convm)                       !CALL:  This one calls numeros subroutines
        
    Mm      = Mm+1
    M_conv  = M_conv+1
            
    if (ER .GT. 0.5 .AND. SS .EQ. 1 .AND. M_conv .GT. 10) THEN
        WRITE(4,*)'ER is greater than 0.5, at time',t
        WRITE(4,*)'Therfore we break'
        stop 'ER greater than 0.5'                                      
                
    ELSEIF(ER .GT. 1 .and. t .LT. 1 .AND. M_conv .GT. 10)THEN
        WRITE(4,*)'ER is greater than 0.1, at time',t
        WRITE(4,*)'Therfore we reset the pressure and restart on the next level'
                
        IF(SS .GT. 1) SS=SS/2
        Mm=-1
        
        ! Reset the problem since too fara off with the solution
        CALL INITI(H00)                                                 
        CALL HREE(H00,t,SS, NYs, T40, Tcvot, 0, 0, M_conv)          
        CALL pastupd(0, NX, NN, SS)                                   
        if( lub_temp == 1 .and. SS .LE. 2) CALL Lubrication_temp(SS, NYs, C1, t, 0, sim_pass, iter_lim, iter_cvot, temp_conv, temp_convm)  
        GOTO 444                                                        
                
    ENDIF
            
    if (isnan(ER)) THEN
        WRITE(4,*)'ER is a NaN, at time',t
        WRITE(4,*)'Therfore we break'
        stop 'ER is a NaN'                                               
    ENDIF
                
    ! Convergence check
    IF(M_conv .LT. MMmin )                                                                              GOTO 10         ! Make sure we do at least MMmin itterations at each level.	
    IF(M_conv .LT. MK .AND. ER .GT. ER_ss)                                                              GOTO 10         ! Check for a global convergence of the pressure
    IF(M_conv .LT. MK .AND. max_diff .GT. 20*ER_ss)                                                     GOTO 10         ! Check for a local convegence of the pressure
    IF(M_conv .LT. MK .AND. ABS(DW) .GT. 0.00005  .AND. t .LE. -11)                                     GOTO 10         ! Check for a converged load ballance if time independent.
    IF(M_conv .LT. MK .AND. ABS(DW) .GT. 0.00005  .AND. t .LT. -1 .and. Dw_meth == 1)                   GOTO 10         ! Check for a converged load ballance if time independent.
    IF(M_conv .LT. MK .AND. lim .GT. KK_time )                                                          GOTO 10         ! No pressure truncation allowed on all subiterations during the last iteration
    IF(M_conv .LT. MK .AND. Mm  .LT. 5 )                                                                GOTO 10         ! Do not proceed to close after a lubrication ofset updation
    IF(M_conv .LT. MK .AND. lub_temp == 1 .and. sim_pass == 0 )                                         GOTO 10         ! Calcluate new pressure profile based on updated temperature      
    IF(t .Lt. 0) WRITE(4,*)'temp_conv = ',temp_conv,  'temp_convm = ',temp_convm
    IF(M_conv .lT. MK .AND. SS .LE. 2 .AND. (max(temp_conv, temp_convm) .GT. 0.03 .AND. t .Lt. 0))      GOTO 10         ! Not converged temperaure                   
    
    ! If the temperature simulations has not had time to run all neded steps but alla ther parameters has converged
    ! Allowe then for 4 times more itterations of the temperature if at last global iteration step
    IF(lub_temp == 1 .and. sim_pass == 0)then                                                                           
        iter_cvot=4                                                                                                     
        Call Lubrication_temp(SS, NYs, C1, t, 1, sim_pass, iter_lim, iter_cvot, temp_conv, temp_convm)         
    endif
            
444 continue 
    
    !If we're not successfull in finding the new pressure. Then stor this failure in the success vector 
    IF(ER .GT. ER_ss .OR. max_diff .GT. 20*ER_ss .or. Mm .EQ. -1)THEN                                       
        if(Ntime .lt. 800) then
            success(iteration+1)=0
        elseif (Ntime .lt. 1600) then
            if(MOD(t,2)==0) success(iteration+1)=0
        endif
        WRITE(4,*)'Fail at time ',t
    ELSE
        WRITE(4,*)'SUCCESS :D at time ',t
    ENDIF
        
    IF( lub_temp == 1 ) then
        IF(sim_pass ==1)then
            WRITE(4,*)'SUCCESS :D for temperature '
        else
            WRITE(4,*)'Fail :( for temperature '
            if(Ntime .lt. 800) then
                ! Reducing the success vector with 0.5 so it can be identified of the temperature filed was solved or not. 
                success(iteration+1)=success(iteration+1)-0.5
            elseif (Ntime .lt. 1600) then
                if(MOD(t,2)==0) success(iteration+1)=success(iteration+1)-0.5
            endif
        endif
    endif
    
    ! The rough surface simulations need more time steps since the rough surfaces have a longer surface to be analised. 
    IF( asp_shape .LT. 200 .or. asp_shape .GT. 209) then
        if( t .EQ. Ntime)then !Finished
            WRITE(*,*)'Finished with the simulation :D'
            WRITE(4,*)'Finished with the simulation :D'
            WRITE(4,*)'Sucsess matrix is:'
            WRITE(4,*) success(1:t)
        endif
    else
        if( t .EQ. Ntime + (135-15)*Ntime/Nx)then !Finished
            WRITE(*,*)'Finished with the simulation :D'
            WRITE(4,*)'Finished with the simulation :D'
            WRITE(4,*)'Sucsess matrix is:'
            WRITE(4,*) success(1:t)
        endif
    endif
        
    ! Check that the output file is not to large. If it is, then something probably have gone wrong. 
    if(.not. single_step_only) Call Stop_to_large_out(t)                 
    return        
        
    END
    