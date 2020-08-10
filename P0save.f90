 ! Subrutine for saving the time independent pressure
	SUBROUTINE P0save()
    	implicit none
        include     'inc_CurrentP.h'
        include     'inc_Grid.h'
        include     'inc_InitialP.h'
        integer     I, J
        SAVE        /InitialP/
        

        P0 = P                      !D_IN: /CurrentP/ -> P         !D_OUT: P0 -> /InitialP/

        
        
        RETURN
    END