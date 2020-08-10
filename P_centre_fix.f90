! A subroutine to increase the pressure of the center line. Since sometimes the pressure gets a too low value at the asperity centre. 
    ! This should only be used for asperities which are highest at the centre line. 
    Subroutine P_centre_fix
    implicit none
    include     'inc_CurrentP.h'
    include     'inc_Grid.h'
    integer     NN, I
    Save    /CurrentP/
    
    NN=(NY+1)/2       !D_IN: /Grid/ -> Ny
    
    DO I=1,NX         !D_IN: /Grid/ -> Nx
        IF( ( P(I,NN) .LT. P(I,NN-1)) .and. (P(I,NN) .GT. 0.8) ) P(I,NN)=P(I,NN-1)  !D_OUT: P -> /CurrentP/
    ENDDO
    
    
    Return
    end
    
    