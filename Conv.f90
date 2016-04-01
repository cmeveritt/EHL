! Generates a converged solution with help of multigrid methods
    Subroutine Conv(t, kvot, DT2, H00, iteration)
        COMMON      /Itt/MK_stat, MK_time, ER_stat, ER_time, Ntime, KK          ! Itteration paramters
        COMMON      /G0DT/G0,DT                                              ! G0 and DT
        real    ER_stat, ER_time
        real    GO, DT
        integer Ntime, MK_stat, MK_time, KK
        integer t, iteration
        real    kvot, DT2, H00
    real, dimension(:), allocatable ::      success
    allocate(       success(Ntime+15))
    
        !To count successfull itterations
        success(1:Ntime+15)=1.0
        
        

    
        M=0                                                                 ! Start counter for whole itteration loop
        MM=0                                                                ! Start counter for internal loop for convergence
        DWold=0.0
10          CALL ITER(t,H00,DW,kvot,DT2)                                    ! Itterate a new pressure KK times

            CALL ERP(ER)                                                    ! Call ERP to calculate how much the pressure changed
            WRITE(4,*)'ER=',ER,'M=',M,'t=',t,'H00=',H00, 'DW', DW           ! Write it uot both to file and comand window
            WRITE(*,*)'ER=',ER,'M=',M,'t=',t,'H00=',H00, 'DW', DW
            
            IF( t .LE. 0) THEN
                DDW=abs((DW-DWold)/(DW**2))
                DWold=DW
                IF(DDW .LT. 0.1 .AND. (ER .LT. 0.01 .OR. Mm .GT. 20)) THEN  ! If the Pressure is good enough or we done enough internal itterations, uppdate the offset H00 
                    CALL H00upd(H00, t, DW)                                 ! Uppdate the pressure based on the loadbalance
                    CALL HREE(H00,t)                                        ! Update lubrication parameters
                    CALL pastupd                                            ! Update the parameters of the past time. 
                ENDIF
            
                IF(ER .GT. 0.01 .AND. MM .GT. 30) THEN                      ! Reset the pressure and the other params to Hertzian pressure since not converged enough
                    CALL INITI
                    CALL HREE(H00,t)             
                    CALL pastupd                                            
                    Mm=0                                                    ! Reset internal counter
                ENDIF
            ENDIF
    
            Mm=Mm+1
            M=M+1
            
            IF(M .LT. MK_stat .AND. ER .GT. ER_stat)GOTO 10             ! Itterate again untill the pressure converged and the load is in ballance, or we run out of itterations
            IF(M .LT. MK_stat .AND. DW .GT. 0.0001  )GOTO 10
            IF(M .LT. 22 )GOTO 10                                       ! Make sure we do at least 21 itterations at each level. 
            
        IF(ER .GT. ER_stat)THEN                                         !Print if we're succesfull or not
            success(iteration+1)=0
            WRITE(*,*)'Fail at time ',t
            WRITE(4,*)'Fail at time ',t
        ENDIF
        
        if (isnan(ER)) stop 'ER is a NaN'                               ! Stop if ER is NaN
        CALL StressXZ                                                   ! Calculate the stresses in the xz plane
        CALL OUTPUT                                                     ! Converged or ran out of tries. Output whatever we have       
        
        
    END
    