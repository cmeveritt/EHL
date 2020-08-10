 ! Subrutine for generating a good new pressure guess by moving the pressure peaks one step forward. 
    ! Could prabably be improved by also slightly increasing the pressures in the inlet and decreasing then in the outlet, especcially a bit after the pressure peak. 
	SUBROUTINE Initial_guess()
    	implicit none
        include     'inc_CurrentP.h'
        include     'inc_Error.h'
        include     'inc_Grid.h'
        include     'inc_Itt.h'
        include     'inc_InitialP.h'
        include     'inc_Ref.h'
        include     'inc_Limit_pressure.h' 
        SAVE        /Error/
        ! Input
        integer     NN
        ! Calculations 
        integer     I, J, JJ
        real        Pasp(NX,NY), Pasp2(NX,NY)                                           !D_IN: /Grid/ -> Nx, Ny
        ! Other                  
        real        asp_step                                      
        ! Output
        SAVE        /currentP/
        
        NN=(NY+1)/2
        
        ! Calculate the pressure from the asperity which is the current pressure P, minus P0 which is a converged pressure for the smooth surfaces    
        DO J=1,NN
            DO I=1,NX-5
                Pasp(I,J)=P(I,J)-P0(I,J)                                                !D_IN: /CurrentP/ -> P; /InitialP/ -> P0
            ENDDO
        ENDDO
        
        ! How many nodes does the surface move each time step
        asp_step = F                     
        
        ! Move the preassure peaks the same distance as the surface will move
        DO J=1,NN
            DO I=6,NX-5
                IF( asp_step .LT. 1)then
                    Pasp2(I,J)=(1-asp_step)*Pasp(I,J)+(asp_step)*Pasp(I-1,J)

                elseif( asp_step .LT. 2)then
                    Pasp2(I,J)=(2-asp_step)*Pasp(I-1,J)+(asp_step-1)*Pasp(I-2,J)
                    
                elseif( asp_step .LT. 3)then
                    Pasp2(I,J)=(3-asp_step)*Pasp(I-2,J)+(asp_step-2)*Pasp(I-3,J)   
                    
                elseif( asp_step .LT. 4)then
                    Pasp2(I,J)=(4-asp_step)*Pasp(I-3,J)+(asp_step-3)*Pasp(I-4,J)
                    
                elseif( asp_step .LT. 5)then
                    Pasp2(I,J)=(5-asp_step)*Pasp(I-4,J)+(asp_step-4)*Pasp(I-5,J)
                else 
                    Pasp2(I,J)=Pasp(I,J)                                                !If this big time step, the code will not do anything. 
                endif
                
                P(I,J)=P0(I,J)+Pasp2(I,J)                                               !D_OUT: P(I,J) -> /CurrentP/
                
                If( P(i,j) .LT. 0.0) P(i,j)  = 0.0
                if( P(i,j) .gt. Plim) P(i,j) = Plim
                
                JJ=NY-J+1
                P(I,JJ)=P(I,J)
            ENDDO
        ENDDO
        
        ! Update Pold, which is the reference pressure for the current itteration
        Pold=P                                                                          !D_OUT: Pold -> /Error/
        RETURN
    END