! The solution of a cylinder contact uses fewer lines in the y (TD) direction in the begining, in order to save time finding a converged solution. 
 !   This subroutine extrapolates to the rest of the Y nodes
    ! If too few lines are used the code have had a hard time converging. 
    Subroutine NYinterp(H00, t,NYs, T40, Tcvot)  
    implicit none
    	include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
	    include     'inc_CurrentT_met.h' 
        include     'inc_Grid.h'
        integer SS, I, J, SSN, NYs,k
        real        H00, T40, Tcvot
        integer     t, N_met
        save        /CurrentP/
        save        /CurrentT/
        save        /CurrentT_met/ 
        
        n_met=39
        
        DO J=NYs+1,NY                                                                                                                                  !D_IN: /Grid/ -> Ny
            DO I=1,NX                                                                                                                                  !D_IN: /Grid/ -> Nx
                P(I,J)=P(I,5+J-NYs)             ! Filling out the rest of the nodes with pressures in the Y-direction                                  !D_OUT: P -> /CurrentP/
                temp(I,J)=temp(I,5+J-NYs)                                                                                                              !D_OUT: temp -> /CurrentT/
            ENDDO
        ENDDO
        
        DO k=1,N_met
        DO J=NYs+1,NY
            DO I=1,NX
                temp_ma(I,J,k)=temp_ma(I,5+J-NYs,k)            ! Filling out the rest of the nodes with pressures in the Y-direction                  !D_OUT: temp_ma -> /CurrentT_met/
                temp_mb(I,J,k)=temp_mb(I,5+J-NYs,k)                                                                                                   !D_OUT: temp_mb -> /CurrentT_met/
            ENDDO
        ENDDO
        ENDDO
        
        SS=1
        CAll HREE(H00,t,SS, NY, T40, Tcvot, 0, 0, 0)                                                                                                  !CALL: No comment
        
        
    RETURN
    END
    