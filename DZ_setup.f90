! A subroutine to define the vertical distance between the metal nodes 
    ! DZ_method = 1 defines the vertical distance as doubeling all the time. 
    ! This since in the end it will give the same total vertical distance as used på DZ_method=0
    ! The distance is 80 x DX which is about 400 um for the gear ivestigated investigated by Everitt and Alfredsson 2018 Surface initiation of rolling contact fatigue at asperities considering slip, shear limit and thermal elastohydrodynamic lubrication
    ! DZ_vec(i) is the distance from node i-1 to node i. 
    SUBROUTINE DZ_setup(DZ_method2, b, DX) !D_IN: as arguments Dz_method contact half width and distance between nodes in X)
    implicit none
    include     'inc_DZ_com.h'
    real        DX, b
    integer     DZ_method2, i 
    SAVE        /DZ_com/
    
    DZ_method = DZ_method2
    
    ! Dz is the vertical distance between nodes in the metal
    IF( DZ_method == 1) then
        ! Increase the distance gradually down into the material. An increasing derivative of the distance could also have been used since less and lesss happends down in the material
        DZ_vec(1) = 0.5                 ! 0.5 comes from an average lubrication height. 
        DZ_vec(2) = 0.5                 !D_OUT: DZ_vec -> /DZ_com/
        do i = 3,40
            DZ_vec(i)=DZ_vec(i-1)+0.5
        enddo
        DZ_vec=DZ_vec*1e-6              ! Rescaling to microns
        
    ELSE
        ! Use the same as the horizontal distance   
        DZ_vec=0*DZ_vec+2*DX*b
    ENDIF

    RETURN
    END
    
