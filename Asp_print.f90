! Subroutine for printing out the simulated surface so it can be visulized
	SUBROUTINE Asp_print(DX2, Y0, NN)
    implicit none 
    ! Input
    real            DX2, Y0
    ! Calculation Parameters
    integer         NN
    integer          I, J, JJ
    real            H23, PAI, asph1, H00
    real            position, Ab, Wb, bump,  bump2,bumpL, H0, asp
    ! Other
    real            aa, sign                                               
    DATA            PAI/3.14159265/
    include     'inc_Grid.h'
    include     'inc_Itt.h'
    include     'inc_Asp.h'
    include     'inc_Ref.h'
    include     'inc_CurrentH.h'
    include     'inc_Current.h'
    include     'inc_COMH.h'
    include     'inc_Setup.h'
    include     'inc_Init_H00.h'

    open (unit = 102, file = "Used_Fortran_surface.dat")
    H00=0

    ! Defining an arbitraty position 
    i=floor(3.0*NX/4.0)
    position=X(i)

    ! Choice of asperity chape ------------------------------------------------------------------------------------------------ 
    ! The parameter asp_shape controuls which equations should be used 
    ! References to the shapes are added in the asp_calc subroutine
    IF( asp_shape  .EQ. 1) THEN
        DO  J=1,NN
            DO  I=1,NX 
                asp=asph*(10**(-10*((X(I)-position)/aspw)**2))   !D_IN: /Setup/ -> X(I); /asp/ -> asph, aspw
                bump=asp*cos(2*PAI*((X(I)-position)/aspw))
                
                H0=bump
                H(I,J)=H0+H00    
                
            ENDDO
        ENDDO
        
    ELSE IF(asp_shape .EQ. 3)THEN
        DO  J=1,NN
            DO  I=1,NX 
                Ab=asph
                Wb=aspw
                bump=Ab*10**(-10*((X(I)-position)/Wb)**2)*cos(2*PAI*(X(I)-position)/Wb)
                H0=-bump
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ELSE IF(asp_shape .EQ. 6) THEN   
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( abs(X(I)-position)/aspw .LT. 1.0)THEN                  
                    bump=asph/2*cos((X(I)-position)*PAI/(aspw))+asph/2;
                ENDIF
                H0=-bump
                H(I,J)=H0  +H00
            ENDDO
        ENDDO   
        
    ELSE IF(asp_shape .EQ. 7) THEN 
        DO  J=1,NN
            DO  I=1,NX 
                bumpL=0
                IF ( (X(I)-position)**2/aspw2 .LT. 1.0)THEN
                    bumpL=asph2/2*cos((X(I)-position)*PAI/(aspw2))+asph2/2;             !D_IN: /asp/ -> aspw2, asph2
                ENDIF
                bump=0
                IF ( ((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                H0=-bump-bumpL
                H(I,J)=H0  +H00
            ENDDO
        ENDDO   
        
    ELSE IF(asp_shape .EQ. 8) THEN 
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( 2*((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                bump2=0
                IF ( 2*((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw2 .LT. 1.0)THEN
                    bump2=asph2/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*2*PAI/(aspw2))+asph2/2;
                ENDIF
                H0=-bump-bump2
                H(I,J)=H0  +H00
            ENDDO
        ENDDO   
        
    ELSE IF(asp_shape .EQ. 13)THEN
        DO  J=1,NN
            DO  I=1,NX 
                Ab=asph
                Wb=aspw
                bump=Ab*10**(-10*((X(I)-position)**2+Y(J)**2)/Wb**2)*cos(2*PAI*sqrt((X(I)-position)**2+Y(J)**2)/Wb)                      !D_IN: /Setup/ -> Y(J)
                H0=-bump
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ELSE IF(asp_shape .EQ. 12 .OR. asp_shape .EQ. 122)THEN   
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( 2*((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN                                                !D_IN: /asp/ -> asp_ratio
                    bump=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ELSE IF(asp_shape .EQ. 112)THEN    
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( 2*((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN                      
                    asph1=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2
                    asph1=asph1*(2.0)**(0.5)                                                                    
                    aa=X(I)-position
                    IF ( 2*abs(aa)/aspw .LT. 0.5)THEN                                                          
                    bump=asph1*sin(aa*2*PAI/(aspw));
                    ELSE                                                                                        
                        IF (X(I) .LT. Position) then
                            bump=asph1/2*cos(2*aa*2*PAI/(aspw))-asph1/2                                        
                        else
                            bump=-asph1/2*cos(2*aa*2*PAI/(aspw))+asph1/2 
                        endif
                    endif
                ENDIF
                H0=-bump
                H(I,J)=H0   +H00
                
            ENDDO
        ENDDO 
        
    ELSE IF(asp_shape .EQ. 32)THEN  
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( 2*((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2
                ELSE IF ( 2*((X(I)-position+aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+1*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2
                ELSE IF ( 2*((X(I)-position+2*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+2*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2 
                ELSE IF ( 2*((X(I)-position+3*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+3*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2
                ELSE IF ( 2*((X(I)-position+4*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+4*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2
                ENDIF
                H0=-bump
                H(I,J)=H0   +H00
                
            ENDDO
        ENDDO 
        
    ELSE IF(asp_shape .EQ. 33)THEN    
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( 2*((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw)) 
                ELSE IF ( 2*((X(I)-position+aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+1*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))
                ELSE IF ( 2*((X(I)-position+2*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+2*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))
                ELSE IF ( 2*((X(I)-position+3*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+3*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))
                ELSE IF ( 2*((X(I)-position+4*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position+4*aspw)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))
                ENDIF
                H0=-bump
                H(I,J)=H0   +H00
                
            ENDDO
        ENDDO      
        
    ELSE   
        DO  J=1,NN
            DO  I=1,NX 
                bump=0
                IF ( 2*((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                H0=-bump
                H(I,J)=H0  +H00 
            ENDDO
        ENDDO   
    ENDIF
    
    DO J=1,NN
            JJ=NY-J+1
            DO I=1,NX
                H(I,JJ)=H(I,J)
            ENDDO
    ENDDO
    
    ! Print the surface
    !Addin extra lines so that the surface has the same size as the on printed if rough surface is used
    DO I=1,800
        if (i .lt. NX) then 
            WRITE(102,110) (H(I,J),J=1,NN)
        else
            WRITE(102,110) (H(NX,J),J=1,NN) 
        endif
    ENDDO
        
110     FORMAT(2001(E13.6,1X))
        close(102)
   
    RETURN
    END
