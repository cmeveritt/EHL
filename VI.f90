! Subroutine for calculating the elastic deformation based upon the solutions for point loads on elastic half-spaces
    ! This is a very timeconsuming way of obtaining the elastic deformations. The VI_FFT subroutine is much faster since it uses fast fourier transforms
	SUBROUTINE VI(SS, NYs, NN, t, iter, full)    
        implicit none
        include     'inc_COMAK2D.h'
        include     'inc_Contact_mat.h'
        include     'inc_Current.h'
	    include     'inc_CurrentP.h'
        include     'inc_Grid.h'
        include     'inc_Itt.h'
        include     'inc_Ref.h'
        include     'inc_Pres_ave_param.h'
        include     'inc_P_line_side.h'
        include     'inc_Cores_used.h'
        include     'inc_Output_control.h'
        ! Input
        integer     SS, NYs, NN, iter, t, full
        ! Calculation parameters
        integer     K, L, I, J, IK, JL, JJ
        integer     Js, Je, number_n
        integer     LL, JL1, JL2
        real        H0
        ! Output
        real ::               W_test(NX,NYs), err_array(NX,NYs), dummy(NX,NYs)
        
        save        /current/
        save        /outcontrol/
        
        if(.NOT. VI_collected) then
          open(88, file = 'contact.DAT', status = 'unknown')
          do i = 1,NX,SS
              write(88,109) (contact(i,j),J=1,NN)
          enddo
          write(88,109) full,iter,NYs
          write(88,110) PAIAK
                 
          write(88,110) (P_line(i),i=1,NX)
109       FORMAT(2001(i5,1X))          
110       FORMAT(2001(E13.6,1X))
          
          VI_collected = .true.  
        endif
        
        ! Calculate the elastic deformation from the sidepressures 
        ! if a ball ----------------------------------------------------------------------------------------------
        IF(Geom .EQ. 2 .OR. Geom .EQ. 3 .OR. Geom .EQ. 6) THEN                                      !D_IN: /Ref/ -> Geom
        !$OMP PARALLEL DO                   &
        !$OMP& IF(use_multiple_cores)       &    
        !$OMP& PRIVATE(J,I,L,JL,LL,K,IK,H0) &
        !$OMP& SHARED(NN,SS,NX,contact,full,NYs,W,DX,PAIAK,P,AK)
        DO J=1,NN,SS
            DO I=1,NX,SS                                                                            !D_IN: /Grid/ -> Nx
                H0=0.0
              
                IF( contact(i,j) .ge. 1 .OR. full == 1 ) then   !Update the elastic deformation     !D_IN: /Contact_mat/ -> contact(I,J)
                    !Deformation from pressure inside the area
                    DO L=1,NYs,SS
                        JL=IABS(J-L)/SS
                        LL=L
                        IF( LL .GT. NN) LL = NYs-L+1            ! Make sure to use updated pressures
                        DO K=1,NX,SS
                            IK=IABS(I-K)/SS
                    
                            H0=H0+AK(IK,JL)*P(K,LL)                                                 !D_IN: /CurrentP/ -> P(I,J); /COMAK2D/ -> AK(I,J)
                        ENDDO
                    ENDDO
                    
                W(I,J)=H0*DX*SS*PAIAK  ! See Eq 2.12 in the new book. The scaling factor between cylinder and ball is pi/4. V_cyl=pi/4*V_ball due to the scaling R/a^2   !D_IN: /Ref/ -> PAIAK    !D_OUT: W(I,J) -> /Current/     this is the elastic deformation matrix
                ENDIF

               
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
    ELSE 
        if( full == 1) then
        if( t .LT. -10 .or. &
        (asp_shape .EQ. 6 .and. contact_alg .eq. 1) .or. &                                      !D_IN: /Ref/ -> asp_shape, Contact_alg
        ( MOD(iter,2) == 0  .and. t .LE. 40 .and. p_ave_param ==1 ) .or. &                      !D_IN: /pres_ave_param/ -> p_ave_param
        ( asp_shape   == 1 .and.  MOD(iter,2) == 0 )  .or. & ! The folowing asp_shapes has line asperities affecting the pressure outside
        ( asp_shape   == 3 .and.  MOD(iter,2) == 0 )  .or. &
        ( asp_shape   == 6 .and.  MOD(iter,2) == 0 )  .or. &
        ( asp_shape   == 7 .and.  MOD(iter,2) == 0 )  .or. &
        ( asp_shape == 201 .and.  MOD(iter,2) == 0 ) )Then   ! Update side elastic deformation if in multigrid calculations or during updating process
            !$OMP PARALLEL DO               &
            !$OMP& IF(use_multiple_cores)   & 
			!$OMP& PRIVATE(J,I,L,LL,K,IK,H0)&
            !$OMP& SHARED(NN,SS,NX,NYs,P_line,Wside,AK)
			 DO J=1,NN,SS
                DO I=1,NX,SS
                    H0=0.0
                    DO L=-NX+J-1,1-SS,SS
                        LL=abs((J-L)/SS)                   ! Number of steps between pressure and current node
                        DO K=1,NX,SS
                            IK=IABS(I-K)/SS
                            H0=H0+AK(IK,LL)*(P_line(K))                                         !D_IN: /P_line_side/ -> P_line(K);
                        enddo
                    enddo
                    DO L=NYs+SS,NX+J,SS
                        LL=abs((L-J)/SS)                   ! Number of steps between pressure and current node
                        DO K=1,NX,SS
                            IK=IABS(I-K)/SS
                            H0=H0+AK(IK,LL)*(P_line(K))
                        enddo
                    enddo
                    Wside(i,j)=H0                                                                !D_OUT: Wside(I,J) -> /Current/ 
                ENDDO
            ENDDO
            !$OMP END PARALLEL DO
        ENDIF                                   ! Make use of the already calculated side deformations
        ENDIF
            !$OMP PARALLEL DO                   &
            !$OMP& IF(use_multiple_cores)       & 
		    !$OMP& PRIVATE(J,I,L,JL1,LL,K,IK,H0)&
            !$OMP& SHARED(NN,SS,NX,contact,full,NYs,W,W_test,DX,PAIAK,P,AK,Wside)  !&
            !!$OMP& SCHEDULE(GUIDED, NN/(2*cores) + 1)
            DO J=1,NN,SS
                DO I=1,NX,SS
                       IF( contact(i,j) .ge. 1 .OR. full == 1) then !Update the elastic deformation
                            H0=0.0
                            DO L=-(J-1),Nys-J,SS
                                LL=abs((L)/SS)                  ! Number of steps between pressure and current node
                                JL1=L+J
                                IF( JL1 .GT. NN) JL1 = NYs-JL1+1 ! To make sure only using updated pressure values
                                DO K=1,NX,SS
                                    IK=IABS(I-K)/SS
                                    H0=H0+AK(IK,LL)*(P(K,JL1))
                                enddo
                            enddo
                            H0=H0+Wside(I,J)
                            W(I,J)=H0*DX*SS*PAIAK  ! See Eq 2.12 in the new book of P. Huang. The scaling factor between cylinder and ball is pi/4. V_cyl=pi/4*V_ball due to the scaling R/a^2
                       ENDIF
                ENDDO
            ENDDO
            !$OMP END PARALLEL DO
    ENDIF

       !$OMP PARALLEL DO            &
       !$OMP& IF(use_multiple_cores)& 
       !$OMP& PRIVATE(J,JJ,I)       &
	   !$OMP& SHARED(NN,SS,NYs,NX,W,W_test)
        DO J=1,NN,SS
            JJ=NYs-J+1
            DO I=1,NX,SS
                W(I,JJ)=W(I,J)
            ENDDO
        ENDDO
        !$OMP END PARALLEL DO

        RETURN
    END