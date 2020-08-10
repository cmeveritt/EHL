! Subrutine for controling the updates of the film thicknes, the viscocity and the density based on the given pressure
	SUBROUTINE HREE(H00,t,SS, NYs, T40, Tcvot, k, k_use, M_conv)
        implicit none 
        include     'inc_Current.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_Grid.h'
        include     'inc_Itt.h'
        include     'inc_Limit_pressure.h' 
        include     'inc_Numbers_of_bad.h'
        include     'inc_Outp.h'
        include     'inc_Ref.h'
        include     'inc_Pres_ave_param.h'
        include     'inc_Cores_used.h'
        ! Input
        integer     NYs, SS, k, k_use, M_conv
        real        H00, T40, Tcvot
        Integer     t
        ! Calculations
        integer     NN, I, J, JJ
        integer     number_n
        real        av_p, av_T,  SRR
        integer     position
        ! Output
        integer     test, test2
        SAVE        /Current/                    ! Current timestep
        SAVE        /CurrentP/
        SAVE        /CurrentT/

        NN=(NYs+1)/2
        
        !To update the pressure outside of the contact.
        If(Geom .EQ. 1 .OR. Geom .EQ. 5) then       
            Call Pressure_outside(t,SS, NYs, k, k_use)  !CALL: If geom = 1 || 5; D_OUT: (t = -11, SS = 1, NYs = Ny, k = 1.0, k_use = 1) -> When called by INITI
        endif
  
        ! Check that the pressure is too high or below zero
        DO J=1,NN,SS
            DO I=1,NX,SS
                IF(P(I,J) .GT. Plim)THEN                                                        !D_IN: /CurrentP/ -> P; /Limit_pressure/ -> Plim
                    P(I,J)=Plim                    
                    Bad_p_node=Bad_p_node+2     ! Since occores on both sides of the mirror     !D_OUT: bad_p_node -> /Numbers_of_bad/ 
                elseIF(P(I,J).LT.0.0)THEN
                    P(I,J)=0.0
                    Bad_p_node=Bad_p_node+2
                ENDIF
            ENDDO
        ENDDO
        
        ! Calculate the elastic deformation --------------------------------------------------
        if (.not. use_fft) then
            CALL VI(SS, NYs, NN, t, k, 1)          ! k indicates if we should update the side-preasure deformation or not, the 1 says that we should update all deformations  !CALL: SS, NYs, NN = (NYs+1)/2, t, k, full = 1
        else
            CALL VIFFT(SS, NYs, NN, t)
        endif

    ! Calculate the film thicknes due to the asperity
    CALL Asp_calc(H00, t, SS, NYs, NN)                                                !CALL: H00, t, SS, NYs, NN = NYs+1/2  
    
    ! Ensure that the the minimum film thickness is keept al over
    CALL Contact_routine(H00, t, SS, NYs, NN, Plim)                                   !CALL: H00, t, SS, NYs, NN = NYs + 1/2, Plim
              
    ! Calculate the viscocity, the density and the dimentionles parameter EPS 
    CALL Newtonian(SS, NYs, T40, Tcvot, NN, t, k, M_conv)                             !CALL: SS, Nys, T40, Tcvot, NN, t, k, M_conv
    
    ! If including shear thinning, calculate the reduction here
    ! shear thinning has a different subroutine for every reference model
    ! Non of the shear thinning models are working correctly
    If(     shear_thin .EQ. 1   .and. SS .LE. 2 )then               
        Call Ehret_shear_thin(SS, NYs, T40, Tcvot, NN)     !CALL: SS, NYs, T40, Tcvot, NN
    ELSEIF( shear_thin .EQ. 2   .and. SS .LE. 2) Then
        Call Liu_shear_thin(SS, NYs, T40, Tcvot, NN)       !CALL: SS, NYs, T40, Tcvot, NN
    ELSEIF( shear_thin .EQ. 3   .and. SS .LE. 2) Then
        Call Temp_Liu_shear_thin(SS, NYs, T40, Tcvot, NN, t) !CALL: SS, NYs, T40, Tcvot, NN, t   
    ELSEIF( shear_thin .EQ. 4   .and. SS .LE. 2) Then
        ! Temperature dependent newtonian model so do nothing
    ENDIF

    ! Mirror all updated parameters around the X-axis------------------------------------------------------------------------------------------------------------------
67  DO J=1,NN-ss,SS
            JJ=NYs-J+1
            DO I=1,NX,SS
                P(I,JJ)=P(I,J)                                  !D_OUT: P(I,J) -> /CurrentP/
                RO(I,JJ)=RO(I,J)                                !D_OUT: R0(I,J) -> /Current/
                xi(I,JJ)=xi(I,J)                                !D_OUT: XI(I,J) -> /Current/
                EDAx(I,JJ)=EDAx(I,J)                            !D_OUT: EDAx(I,J) -> /Current/
                EDAy(I,JJ)=EDAy(I,J)                            !D_OUT: EDAy(I,J) -> /Current/
                EPSx(I,JJ)=EPSx(I,J)                            !D_OUT: EPSx(I,J) -> /Current/
    	        EPSy(I,JJ)=EPSy(I,J)                            !D_OUT: EPSy(I,J) -> /Current/
            ENDDO
        ENDDO

    RETURN
    END
