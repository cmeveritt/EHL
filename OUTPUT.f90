! Subroutine for printing the solution at each time step
	SUBROUTINE OUTPUT
        implicit none
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
	    include     'inc_CurrentT_met.h'
        include     'inc_CurrentRO.h'
        include     'inc_Grid.h'
        include     'inc_Residual.h'
        include     'inc_Setup.h'
        include     'inc_Output_control.h'
        include     'inc_Single_step.h'
        
        integer     NN, I, J, Nz, k, Depth_number
        integer     NN_m, NY_m, NX_m
        real        A

        NY_m=(ny+1)/2                    !D_IN: /Grid/ -> Ny 
        NX_m=(nx+1)/2                    !D_IN: /Grid/ -> Nx
        NN_m=(NY_m+1)/2
        NN=(NY+1)/2
        A=0.0
        
        if (.not. single_step_only) then 
        
            WRITE(8,110)A,(Y(I),I=1,NY)     !D_IN: /Setup/ -> Y(I)
            DO I=1,NX
                WRITE(8,110)X(I),(H(I,J),J=1,NY)      !D_IN: /CurrentH/ -> H(I,J)
            ENDDO
        
            WRITE(10,110)A,(Y(I),I=1,NY)
            DO I=1,NX
                WRITE(10,110)X(I),(P(I,J),J=1,NY)        !(10 is the file number, 110 is the format number)     !D_IN: /Setup/ -> X(I); /CurrentP/ -> P(I,J)
            ENDDO
        
            ! Print oil temperature
            WRITE(23,110)A,(Y(I),I=1,NY)
            DO I=1,NX
                WRITE(23,110)X(I),(temp(I,J),J=1,NY)        !(10 is the file number, 110 is the format number)  !D_IN: /CurrentT/ -> temp(i,j)
            ENDDO
        
            ! Print the temperature at the symmetry line
			Depth_number = 30
            WRITE(50,110)A,(Y(I),I=1,NX)
            DO K=1,Depth_number
                WRITE(50,110)X(Depth_number+1-K),(temp_ma(I,NN_m,Depth_number+1-k),I=1,NX)        !(10 is the file number, 110 is the format number)   !D_IN: /CurrentT_met/ -> Temp_ma   ! ??? What is I here?
            ENDDO
            WRITE(50,110)X(K),(temp(I,NN),I=1,NX)        !(10 is the file number, 110 is the format number)
            DO K=1,Depth_number
                WRITE(50,110)X(K),(temp_mb(I,NN_m,k),I=1,NX)        !(10 is the file number, 110 is the format number)         !D_IN: /CurrentT_met/ -> Temp_mb
            ENDDO
        
            ! Print the temperature in the body witht the asperity
            if(collect_full_data) then
                DO K=1,39
                    DO I=1,NX_m
                        WRITE(51,110)(temp_ma(I,J,k),J=1,NY_m)        !(10 is the file number, 110 is the format number)
                        WRITE(54,110)(temp_mb(I,J,k),J=1,NY_m)
                    ENDDO
                ENDDO
            
            
                DO I = 1,NX
                    WRITE(55,110) (EDAx(I,J),J=1,NY)
                    WRITE(56,110) (EDAy(I,J),J=1,NY)
                    WRITE(57,110) (EPSx(I,J),J=1,NY)
                    WRITE(58,110) (EPSy(I,J),J=1,NY)
                    WRITE(59,110) (RO(I,J),J=1,NY)
                    WRITE(61,110) (xi(I,J),J=1,NY)
                    WRITE(62,110) (w(I,J),J=1,NY)
                    WRITE(63,110) (dRodP_mat(I,J),J=1,NY)
                ENDDO 
        
            Else
            
                DO K=1,39
                    WRITE(51,110)A,(Y(I*2-1),I=1,NN_m)
                    DO I=1,NX_m
                        WRITE(51,110)X(I*2-1),(temp_ma(I,J,k),J=1,NN_m)        !(10 is the file number, 110 is the format number)
                    ENDDO
                ENDDO
            
            Endif
            
        Else
            !---Film thickness---
            WRITE(65,110)A,(Y(I),I=1,NY)    
            DO I=1,NX
                WRITE(65,110)X(I),(H(I,J),J=1,NY)    
            ENDDO
             !----Pressure ----
            WRITE(64,110)A,(Y(I),I=1,NY)
            DO I=1,NX
                WRITE(64,110)X(I),(P(I,J),J=1,NY)
            ENDDO
            !------film temperature------
            WRITE(66,110)A,(Y(I),I=1,NY)
            DO I=1,NX
                WRITE(66,110)X(I),(temp(I,J),J=1,NY)         !D_IN: /CurrentT/ -> temp(i,j)
            ENDDO
            !------reduced metal A temperatue matrix----
            DO K=1,39
                    WRITE(67,110)A,(Y(I*2-1),I=1,NN_m)
                    DO I=1,NX_m
                        WRITE(67,110)X(I*2-1),(temp_ma(I,J,k),J=1,NN_m)       
                    ENDDO
            ENDDO
        Endif
            
        
        
    110 FORMAT(2001(E13.6,1X))
        RETURN
    END