! Subrutine for the main structure of the program
    SUBROUTINE Time(H00)
        implicit none          
        COMMON      /Grid/NX,NY,X0,XE,DX                                        ! Grid parameters
        COMMON      /Ref/ref, PAIAK, tmeth                                           ! Coise of equations based on referense
        COMMON      /Itt/MK_stat, MK_time, ER_stat, ER_time, Ntime, KK          ! Itteration paramters
        COMMON      /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0           ! Lubrication parameters
        COMMON      /HREEp/Emat, pg, Tg,YF                                       ! Parameters for HREE subrutine
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast  ! Past
        COMMON      /Setup/X,Y, P0                                       ! Parameters whihc do not change over time
        COMMON      /Residual/A5dxx,A5dx,A5dt                                   ! Parameters for the residual
        COMMON      /G0DT/G0,DT                                              ! G0 and DT
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        real         X(1:300),Y(1:300), P0(1:300,1:300)
        real         A5dxx(1:300,1:300),A5dx(1:300,1:300),A5dt(1:300,1:300)
        real         Emat(1:300,1:300) ,pg(1:300,1:300), Tg(1:300,1:300), YF(1:300,1:300)
        integer         t, Mm, M, iteration, tmeth
        real            isnan, H00, G0, ER, DW, DWold, DDw, DT, DTT
        integer         NX, NY
        real            X0,XE,DX                                                ! Grid parameters
        integer         ref
        real            PAIAK, abs                                          ! Coise of equations based on referense
        integer         Ntime, MK_stat, MK_time, KK
        real            ER_stat, ER_time                                        ! Itteration paramters
        real            kvot,DT2
        real            ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
        

        ! Ensure they are zero from the start
        A5dxx   =A5dxx*0.0
        A5dx    =A5dxx*0.0
        A5dt    =A5dxx*0.0
        xi      =A5dxx*0.0
        xipast  =A5dxx*0.0
       
        Iteration=1
        
        ! Starting calculations	
        ! Calculate the deformation matrix AK
        CALL SUBAK
        
        ! Setup grid and start guess. 
        CALL INITI
        CALL VIside
        CALL Pastupd
        
        ! ----------------------------------------------------------------------------------
        ! --------------Calculate static sulution------------------------------------------
        !-----------------------------------------------------------------------------------
  
        t=-1                                                            ! Controles asperity position
        CALL HREE(H00,t)                                         ! Updates Lubrication Parameters given a pressure and a lubrication offset
        CALL Pastupd 
        kvot=1.0                                                          ! How much weigth is put on the current time step. For static, kvot =1, for timedep acc C-N kvot=0.5      
        DT2=0.0                                                         ! =1/(2.*DT). For timeindependent DT=inf so DT2=0 
                                                                       ! Start counter for internal loop for H00 updation
        Call Conv(t, kvot, DT2, H00, iteration)
        
        ! ----------------------------------------------------------------------------------
        ! --------------Gently add the timedependence ------------------------------------------
        !-----------------------------------------------------------------------------------
        t=0                                                             ! Jump one timestep for numerical reasons
        kvot=1.0
        DO iteration = 1,10                                               ! Add the timedependance over a certain number of steps. If changing here, remember to switch length oc success
            DTT=1000.0/(iteration**3)*DT                                     ! Gradually increasing the timeaffects. 
            DT2=1./(2*DTT)
            
                CALL Pastupd                                            ! Update parameters of past timestep
                
                Call Conv(t, kvot, DT2, H00, iteration)
            
        ENDDO 
        
        
        ! ----------------------------------------------------------------------------------
        ! --------------Calculate time dependent  sulution---------------------------------
        !-----------------------------------------------------------------------------------
        IF( tmeth .GT. 1) kvot=0.5
        DT2=1.0/(2*DT)
        
        ! For all timesteps, calculate the timedepending solution
        DO t=0,Ntime
            
            iteration=iteration+1
            ! Stor solution from previus timestep
            CALL Pastupd
            
            ! Solve the current timestep
            Call Conv(t, kvot, DT2, H00, iteration)
        
        ENDDO
               
        RETURN
    END
! **********************************************************!************************************************************************************************************************************************