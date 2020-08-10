! Subroutine for reading the data of a specific step. Used if restarting a simulation at a specific time. 
    SUBROUTINE Get_single_step_data(t, H00)
    implicit none
    include     'inc_Single_step.h'
    include     'inc_Current.h'
    include     'inc_CurrentT.h'
    include     'inc_CurrentH.h'
    include     'inc_CurrentP.h'
    include     'inc_CurrentRO.h'
    include     'inc_Init_H00.h'
    include     'inc_Itt.h'
    include     'inc_Grid.h'
    include     'inc_Ref.h'
    include     'inc_InitialP.h'
    include     'inc_CurrentT_met.h'
    include     'inc_Past.h'
    integer     dummy, p0_index, NX_m, NY_m, i, j, k
    integer     t 
    real        H00
    ! Output
    save        /Init_H00/
    save        /InitialP/
    save        /Past/
    save        /CurrentP/
    save        /CurrentH/
    save        /CurrentT/
    save        /Current/
    save        /CurrentRO/
    save        /CurrentT_met/
    
    !------------------------------------------------------------------------------
    !----------------Opening unopened but needed files-----------------------------
    !------------------------------------------------------------------------------
    OPEN(60, FILE = 'Film_initial.DAT', STATUS = 'UNKNOWN')    
    OPEN(54,FILE='MetalB_temp.DAT',STATUS='UNKNOWN')
    OPEN(55, FILE = 'EDAx.DAT', STATUS = 'UNKNOWN')
    OPEN(56, FILE = 'EDAy.DAT',STATUS = 'UNKNOWN')
    OPEN(57, FILE = 'EPSx.DAT', STATUS = 'UNKNOWN')
    OPEN(58, FILE = 'EPSy.DAT', STATUS = 'UNKNOWN')
    OPEN(59, FILE = 'RO.DAT', STATUS = 'UNKNOWN')
    OPEN(61, FILE = 'xi.DAT', STATUS = 'UNKNOWN')
    OPEN(62, FILE = 'w.DAT', STATUS = 'UNKNOWN')
    OPEN(63, FILE = 'drodP_mat.DAT', STATUS = 'UNKNOWN')
    
    !------------------------------------------------------------------------------
    !----------------Calculating required numbers----------------------------------
    !------------------------------------------------------------------------------
    NX_m = (NX + 1)/2
    NY_m = (NY + 1)/2
    
    if(multi_grid_param .NE. 1) then
        p0_index = 11 + 1 - 1
        t = selected_step - 11 - 1
        if(t > 10 .and. Ntime .GT. 100) p0_index = p0_index + 11 
    else if (NX .gt. 200 .and. NY .gt. 160) then
        p0_index = 11 + 4 - 1
        t = selected_step - 11 - 4
        if(t > 10 .and. Ntime .GT. 100) p0_index = p0_index + 11
    else 
        p0_index = 11 + 3  - 1
        t = selected_step - 11 - 3
        if(t > 10 .and. Ntime .GT. 100) p0_index = p0_index + 11
    endif
    
    !------------------------------------------------------------------------------
    !----------------Moving reading pointers to the selected data------------------
    !------------------------------------------------------------------------------
    
    do i = 1,NX * (p0_index - 1) + p0_index
        read(10,*) !This is Pressure.DAT
    ENDDO
    
    do i = 1,NX * (selected_step - 3) + selected_step - 2
        read(8,*) !This is Film.DAT
    ENDDO
    
    do i = 1,NX * (selected_step - 2) + selected_step - 1
        read(23,*) !This is Temp.DAT
    ENDDO
    
    do i = 1,NX * (selected_step - 2) 
        read(55,*) !EDAx.DAT
        read(56,*) !EDAy.DAT
        read(57,*) !EPSx.DAT
        read(58,*) !EPSy.DAT
        read(59,*) !RO.DAT
        read(61,*) !xi.DAT
        read(62,*) !w.DAT
        read(63,*) !dRodP_mat.DAT
    ENDDO
    
    DO K=1,39 * (selected_step - 2)
        DO I=1,NX_m
            read(51,*) !MetalA_temp.DAT
            read(54,*) !MetalB_temp.DAT
        ENDDO
    ENDDO
    
    !------------------------------------------------------------------------------
    !----------------Reading data--------------------------------------------------
    !------------------------------------------------------------------------------
    read(60,*) H00_t0
    
    H00 = H00_t0
    
    !read P0
    do i = 1,Nx
        read(10,*)Dummy,(P0(i,j),j=1,NY)
    ENDDO
    
    !adjust the reading index again on Pressure.dat
    do i = 1,NX * ((selected_step - p0_index - 1) - 2) + (selected_step - p0_index - 2)
        read(10,*)
    ENDDO
    
    !Reading 1 step in the past
    do i = 1,Nx
        read(10,*) Dummy,(Ppast(i,j),j=1,NY)
        read(8,*) Dummy,(Hpast(i,j),j=1,NY)
    ENDDO
    !Pold = Ppast
    
    !move index one line down
    read(10,*)
    read(8,*)
    
    !Read all 2D arrays for the current time step
    do i = 1,Nx
        read(10,*) Dummy,(P(i,j),j=1,NY)
        read(8,*) Dummy,(H(i,j),j=1,NY)
        read(23,*) Dummy, (temp(i,j),j=1,NY)
        read(55,*) (EDAx(I,J),J=1,NY)
        read(56,*) (EDAy(I,J),J=1,NY)
        read(57,*) (EPSx(I,J),J=1,NY)
        read(58,*) (EPSy(I,J),J=1,NY)
        read(59,*) (RO(I,J),J=1,NY)
        read(61,*) (xi(I,J),J=1,NY)
        read(62,*) (w(I,J),J=1,NY)
        read(63,*) (dRodP_mat(I,J),J=1,NY)
    ENDDO
    
    !Read the 3D arrays for the current time step
    DO K=1,39
            DO I=1,NX_m
                read(51,*)(temp_ma(I,J,k),J=1,NY_m)
                read(54,*)(temp_mb(I,J,k),J=1,NY_m)
            ENDDO
    ENDDO
    
    110 FORMAT(2001(E13.6,1X)) 
    RETURN
END