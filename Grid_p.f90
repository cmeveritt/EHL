! Subroutine for updating gridparameters which are depending on the stepsize ss
Subroutine Grid_p(SS)
    implicit none
    include     'inc_Grid.h'
    include     'inc_Grid2.h'
    ! Input
    integer     SS
    ! Output
    SAVE       /Grid2/ 
    
    ! Grid parameters
    DX1=1./DX*1./SS    !D_IN: /grid/ -> Dx    !D_OUT: DX1 -> /Grid2/
    DX2=DX*DX*SS*SS                           !D_OUT: DX2 -> /Grid2/
    DX3=1./DX2                                !D_OUT: DX3 -> /Grid2/
    DX4=0.3*DX2             ! Only used for determening if to use Gauss-Sidel or Jacobi method  !D_OUT: DX4 -> /Grid2/
        
    END
    
        