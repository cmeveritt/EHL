! Subroutine for reading the data of a previously measured rough surface
! The format of the imported surface has to be 135x134 data points stating the height in um at each position. 
! This is only executed if asp_shape in input8.csv has a certain value
	SUBROUTINE Asp_read(NX,NY, Rx, b) !D_IN: arguments -> number of nodes in X and Y, Radius of the rolling object, and contact half width
        implicit none 
        real                        Asp_surf_r(135,134), Asp_surf(155,134)
        real                        factori, factorj
        real                        Rx, b
        integer                     Nx_asp, Ny_asp, i, j, NX, NY,NN
        integer                     jj, ii, j2
        include     'inc_Rough_surf.h'
        include     'inc_Ref.h'
        SAVE        /Rough_surf/
        
    open (unit = 101, file = "Fortran_surface.csv")
    open (unit = 102, file = "Used_Fortran_surface.dat")
    
    Nx_asp = 135  ! Number of colums in file
    NY_asp = 134  ! Number of rows in the file
    NN     = (NY+1)/2
    
    ! Read the data to Asp_surf
    do i = 1,Ny_asp               
        read(101,*) Asp_surf_r(:,i)
    enddo
    
    ! Turn the surface upside down, since bumps decrease the lubrication height, and scale it
    Asp_surf_r=-Asp_surf_r* surf_scale !D_IN: /Rough_surf/ -> surf_scale
    
    ! Extract the wanted data and smoothen the data at the borders
    do j = 1, NN
        jj=j-4
        if( j .le. 4) jj=1
        factorj = 1

        do i = 1, Nx_asp+20
            factori = 1
            ii      = i-10
            if( i .le. 10)       ii=1
            if( ii .gt. 135)    ii=135
            
            ! Slowly introduce the surface roughness from the side.
            if( i .lt. 20 ) then
                factori = i*1.0/20.0
            elseif( i .gt. Nx_asp ) Then
                factori = (Nx_asp + 20 - i) * 1.0/20.0
            endif
            
            ! Original surface rotaion
            if( asp_shape == 200 .or. asp_shape == 201) then
                Asp_surf(i,j) = asp_surf_r(ii,jj+shift_y) * factori !D_IN: /Ref/ -> shift_y
                
            ! Changing the rolling direction
            elseif( asp_shape == 202 .or. asp_shape == 203 ) then
                Asp_surf(i,j) = asp_surf_r(135+1-ii,jj+shift_y) * factori
                
            ! Rotate the surface 90 degrees    
            elseif( asp_shape == 204) then   
                if(ii .gt. 134) ii=134      ! To not end up outside the matrix since asp_surf_r is not equal in x and y. 
                Asp_surf(i,j) = asp_surf_r(jj+shift_y,134+1-ii) * factori
                
            ! Rotate the surface 90 degrees   but run backwards 
            elseif( asp_shape == 205) then   
                if(ii .gt. 134) ii=134      ! To not end up outside the matrix sinse asp_surf_r is not equal in x and y. 
                Asp_surf(i,j) = asp_surf_r(jj+shift_y,ii) * factori
                
            ! Rotate the surface 90 degrees in the other directions    
            elseif( asp_shape == 206) then   
                if(ii .gt. 134) ii=134
                Asp_surf(i,j) = asp_surf_r(NN+1-jj+shift_y,ii) * factori
            
            ! Rotate the surface 90 degrees in the other directions and run backwards 
            elseif( asp_shape == 207) then   
                if(ii .gt. 134) ii=134
                Asp_surf(i,j) = asp_surf_r(NN+1-jj+shift_y,134+1-ii) * factori
            endif
            
            ! Use a constant second derivative over the three central  nodes. 
            if( J==NN) Asp_surf(i,j) = ( 4* Asp_surf(i,j-1) - Asp_surf(i,j-2)) / 3  
            
        enddo
    enddo
    
    ! Rescaling surface to dimensionless form. The data should have the units um
    Asp_surf=Asp_surf*Rx/(b*b) * 1e-6
   
    do j=1, NN
        ! Slowly introduce the surface roughness from the side. At the centrale line it sould be the full roughness
        if( j .lt. 9 .and. asp_shape .ne. 201 .and. asp_shape .ne. 203) then
            factorj = (j-1.0)*1.0/8.0
        else
            factorj = 1.0
        endif
                
        ! Use constant values of the surface at borders. This is good for ground surfaces
        if( j .LT. 8 .and. (asp_shape == 201 .or. asp_shape == 203))then 
            j2=8
        else
            j2=j
        endif

        do i=1,Nx_asp+20
            JJ=NY-J+1
            Asp_surf_tot(800-NX-Nx_asp+i,j)  = Asp_surf(i,j2)*factorj !D_OUT:  Asp_rough_tot(I,J) -> /Rough_surf/
            Asp_surf_tot(800-NX-Nx_asp+i,jj) = Asp_surf(i,j2)*factorj
        enddo
    enddo
    
    ! Print out the actual used surface
    DO I=1,800
        WRITE(102,110) (Asp_surf_tot(I,J),J=1,NN)
    ENDDO
        
    110 FORMAT(2001(E13.6,1X))
    RETURN
    END
