  ! To update the film thickness offset H00
    ! Initially simpler equations were used boot these have prooven to be more stable
    SUBROUTINE H00upd(H00, t, DW, SS, NYs, H00ink)
        implicit none
        include     'inc_Com_H00.h'
        include     'inc_CurrentH.h'
        include     'inc_Grid.h'
        include     'inc_Ref.h'
        include     'inc_Visc.h'                                                  ! Grid parameters
        ! Input
        integer     NYs,SS
        real        DW
        ! Calculations
        integer     NN, I, J
        real        HMIN
        ! Other                    ! Current timestep
        integer     t
        !Output
        real        H00,  H00ink
        save        /Com_H00/

        NN=(NYs+1)/2
        
        ! The minimum film thickness Hmin is used later to ensure that the film thickness is not decreased too much.  
        ! Start guess and maximum value
        Hmin=0.4                                    
        
        ! Only look for Hmin at used gridpoint for current level. 
        DO J=1,NN,SS
            Do I= 1,NX,SS                                                     !D_IN: /Grid/ -> Nx
                IF( HMIN .GT. H(I,J)) then                                    !D_IN: /CurrentH/ -> H(I,J)
                    HMIN = H(I,J) 
                endif  
            ENDDO
        ENDDO
    
        ! Limit the value of the minimum film thickness.    
        IF( HMIN .LT. 5E-6) HMIN = 5E-6              

        IF( H00 .EQ. H00past) THEN                                            !D_IN: /Com_H00/ -> H00past
            H00ink=DW*HM0r*HMIN                                               !D_IN: /Visc/ -> HM0r
        ELSE  
            H00ink=-DW/(DW-DWpast)*(H00-H00past)                              !D_IN: /Com_H00/ -> DWpast
        ENDIF
        
        ! Rescaling the increment based on how much the load differs from the referential load. 
        IF( Geom .EQ. 5) H00ink=H00ink*((PH*1E-9/2.3893356) )                 !D_IN: /Ref/ -> Geom; /Visc/ -> PH           
            
        ! Different limits if we're close to or fara away from load equilibrium
        IF( abs(DW) .LT. 0.05) then
            IF( H00ink .GT.  10*HMIN)  H00ink=  10*HMIN
            IF( H00ink .LT. -0.9*Hmin) H00ink= -0.9*HMIN
        else
            IF( abs(H00ink) .GT.  0.1) H00ink= sign(0.1,H00ink)
        endif
            
        ! If the lubrication film should decrease, be more carfull.     
        IF( H00ink .LT. -0.9*HMIN) H00ink=-0.9*HMIN
        IF(ABS(H00ink) .GT.  0.5*sqrt(ABS(DW))) H00ink = 0.5*sign(sqrt(abs(DW)),DW)         ! Ensuring not to bigg steps with aspects to the loadballance !SIGN(A,B) returns the value of A with the sign of B
            
        IF(t .GT. -11) H00ink=H00ink*0.1                                                    ! Lower the possibility for H00 updation if in the timedependent simulation.  ! Becouse then the timederivative of H00 will also affect the solution. 
            
        ! Ensuring the right derivative
        IF( DW .GT. 0.0 .and. H00ink .LT. 0.0) H00ink=-H00ink*0.5                           ! Since a decrease in H00 will cause an increase in the load
        IF( DW .LT. 0.0 .and. H00ink .GT. 0.0) H00ink=-H00ink*0.5
            
        ! Uppdation H00
        H00past=H00                                                           !D_OUT: H00past -> /Com_H00/
        H00=H00+H00ink
            
        !Updata others
        DWpast=DW                                                             !D_OUT: DWpast -> /Com_H00/
            
        ! Divergence check
        if (H00 .GT. 20) stop 'Too High H00'                                    
        if (H00 .LT. -5) stop 'Too Low H00'
        
        ! Update the parameters of the past timestep. This is used to get rid of unnessesery time derivatives in the statinary part of the simulation
        If( t .LT. 0) call pastupd(0, NX, NN, SS)                        !CALL: No comment
        RETURN
    END