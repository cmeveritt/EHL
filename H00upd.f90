  ! To update the offset H00
    SUBROUTINE H00upd(H00, t, DW)
        implicit none
        COMMON      /G0DT/G0,DT                                              ! G0 and DT
        COMMON      /COMH/RAD(1:300,1:300)
        COMMON      /H00/ H00past, DWpast, Sloadpast
        COMMON      /Grid/NX,NY,X0,XE,DX                                                        ! Grid parameters
        COMMON      /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0                           ! Lubrication parameters
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                              ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast   ! Past
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        integer     I, J, loadupd, t,NX,NY
        real        H0, H00, H00past, HMIN, H00ink
        real        RAD,X0,XE,DX,DT
        real        Sload, Sloadpast
        real        load, loadpast, G0, DW, DWpast, loaddiff
        real        ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0
        real        SUM, minval
        save        /H00/

            
            !loadupd=0
            HMIN=minval(H(1:NX,1:NY))
    
            IF( HMIN .LT. 5E-6) HMIN = 5E-6
            IF( HMIN .GT. 0.4)  HMIN = 0.4
            
            load=sum(P)                     ! The load
           
            Sload=DX*DX*load/G0             ! Normalized load

            DW=Sload-1.0                    ! 1-scaled loadsum

            !loadDiff=Sload-Sloadpast
            
            ! If the pressure is resonable, use the following to update the filmthickness
            !IF( DW .LT.  0.01 .AND. loadDiff .LE. 0.01)  loadupd=1
            !IF( DW .GT. -0.01 .AND. loadDiff .GE. -0.01) loadupd=1

            
            !IF( loadupd .EQ. 1) THEN
                IF( H00 .EQ. H00past) THEN
                    H00ink=DW*HM0r*HMIN
                    IF( H00ink .GT.  0.9*HMIN) H00ink= 0.9*HMIN
                    IF( H00ink .LT. -0.9*HMIN) H00ink=-0.9*HMIN
                ELSE  
                    H00ink=-DW/(DW-DWpast)*(H00-H00past)      
                    IF( H00ink .GT.  0.9*HMIN) H00ink= 0.9*HMIN
                    IF( H00ink .LT. -0.9*HMIN) H00ink=-0.9*HMIN
                ENDIF
            !ENDIF
            
            ! Limiting the increment
            !IF(DW .LT.  -0.1) H00ink = 0.6*H00ink
            IF(ABS(H00ink) .GT.  0.5*sqrt(ABS(DW))) H00ink = 0.5*sign(sqrt(abs(DW)),DW)           ! Ensuring not to bigg steps with aspects to the loadballance !SIGN(A,B) returns the value of A with the sign of B
            
            !IF(t .GE. 0) H00ink=H00ink*0.1                                                         ! Lower the possibility for H00 updation if in the timedependent simulation. 
                                                                                                    ! Becouse then the timederivative of H00 will also affect the solution. 
            
            ! Ensuring the right derivative
            IF( DW .GT. 0.0 .and. H00ink .LT. 0.0) H00ink=-H00ink*0.5                               ! Since a decrease in H00 will cause an increase in the load
            IF( DW .LT. 0.0 .and. H00ink .GT. 0.0) H00ink=-H00ink*0.5
            
            ! Uppdation H00
            H00past=H00
            H00=H00+H00ink
            
            !Updata others
            DWpast=DW
            Sloadpast=Sload
            
            ! Divergence check
            if (H00 .GT. 5) stop 'Too High H00'                                    
            if (H00 .LT. -5) stop 'Too Low H00'
            
            
        RETURN
    END