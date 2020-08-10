! Calculate the dimensionless film thickness, H, based on the previusly read surface, the curvature, RAD, the elastic deformation, W, and the film thickness offset, H00
	SUBROUTINE Asp_surf_calc(NX,NYs, NN,SS, t, H00)
    implicit none 
    ! Input
    real            position, H00
    integer         NX, NYs, NN, SS, F, T
    ! Calculations
    real            H0
    integer         i, j, jj
    include     'inc_Asp.h'
    include     'inc_COMH.h'
    include     'inc_Current.h'
    include     'inc_CurrentH.h'
    include     'inc_Rough_surf.h'
    ! Output
    save            /currentH/

    do j = 1, NN
        H0 = 0
        JJ=NYs-J+1
        do i = 1, NX
            H0      = Asp_surf_tot(800-Nx-t+i,j )                             !D_IN: /Rough_surf/ -> Asp_surf_tot(800,301))
            H(I,J)  = RAD(I,J) + W(I,J)  + H00 + H0                           !D_IN: /Current/-> W(I,J); /COMH/ -> RAD(I,J);                                      !D_OUT: H(I,J) -> /CurrentH/
            H(I,JJ) = H(I,J)
        enddo
    enddo

    RETURN
    END
