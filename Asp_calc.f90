! Subroutine for calculating the shape of the surface roughness. 
    !The shape is either a anayltical equatio och measurements of a surface
	SUBROUTINE Asp_calc(H00, t, SS, NYs, NN)
    implicit none 
    ! Input
    real            H00
    integer         t, SS, NYs, NN
    ! Calculation Parameters
    integer          I, J, JJ
    real            H23, PAI, asph1
    real            position, Ab, Wb, bump, bump2, bumpL, H0, asp
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
    ! Output
    save            /currentH/
    save            /Init_H00/
        
    ! Calculate the dimensionless position of the surface defect. 
    position=X0+(XE-X0)*t/Ntime ! = Ua/Um*t/Ntime*Tend                   !D_IN: /Grid/ -> XE, X0; /Itt/ -> Ntime

    ! This part is used for simulatios where the load increases over the time. This part is not fully developed yet
    H23=-2.0 ! This is a constant for how much the film thickness hsould change. Idealy whould be to have it as in input parameter
    If( asp_shape .EQ. 122)then                                          !D_IN: /Ref/ -> asp_shape
        if( t .EQ. 0) H00_t0=H00                                         !D_OUT: H00_t0 -> /Init_H00/
        if( t .GT. 0) then
            H00=H00_t0+(H23-H00_t0)*t/Ntime
        endif
    endif
    
    ! Choice of asperity chape ------------------------------------------------------------------------------------------------ 
    ! The parameter asp_shape controuls which equations should be used 
    
    ! No asperity 
    IF(asp_shape .EQ. 0 .OR. t .LT. -10)THEN              
        DO  J=1,NN,SS
            DO  I=1,NX,SS                                                !D_IN: /Grid/ -> NX
                H0=RAD(I,J)+W(I,J)                                       !D_IN: /COMH/ -> RAD; /Current/ -> W
                H(I,J)=H0   +H00                                         !D_OUT: H(I,J) -> /CurrentH/                
            ENDDO
        ENDDO
        
    ! Asperity shape based upone Venner CH and Lubrecht AA. Transient analysis of surfacefeatures in an EHL line contact in the case of sliding. Trans ASME J Tribol 1994; 116: 186–193.
    ELSE IF( asp_shape  .EQ. 1) THEN
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                asp=asph*(10**(-10*((X(I)-position)/aspw)**2))          !D_IN: /Setup/ -> X(I); /asp/ -> asph, aspw
                bump=asp*cos(2*PAI*((X(I)-position)/aspw))
                H0=RAD(I,J)+W(I,J)+bump
                H(I,J)=H0+H00    
            ENDDO
        ENDDO
        
    ! Input parameters acc Holmes et al. Transient EHL point contact analysis 2003. Also the same as Venner 1994 ball and cylinder  
    ELSE IF(asp_shape .EQ. 3)THEN
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                Ab=asph
                Wb=aspw
                bump=Ab*10**(-10*((X(I)-position)/Wb)**2)*cos(2*PAI*(X(I)-position)/Wb)
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ! The asperity shape used by D. Hannes in D. Hannes 2011 Rolling contact fatigue crack path prediction but formulates as a line defect
    ELSE IF(asp_shape .EQ. 6) THEN   
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                bump=0
                IF ( abs(X(I)-position)/aspw .LT. 1.0)THEN                  
                    bump=asph/2*cos((X(I)-position)*PAI/(aspw))+asph/2;
                ENDIF
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0  +H00
            ENDDO
        ENDDO   
    
    ! The asperity shape used by D. Hannes in Rolling contact fatigue crack path prediction 2011 but formulates as a line defect plus on as a point defect
    ELSE IF(asp_shape .EQ. 7) THEN   
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                ! line
                bumpL=0
                IF ( (X(I)-position)**2/aspw2 .LT. 1.0)THEN
                    bumpL=asph2/2*cos((X(I)-position)*PAI/(aspw2))+asph2/2;             !D_IN: /asp/ -> aspw2, asph2
                ENDIF
                ! point
                bump=0
                IF ( ((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                H0=RAD(I,J)+W(I,J)-bump-bumpL
                H(I,J)=H0  +H00
            ENDDO
        ENDDO   
    
    ! Two times the asperity shape used by D. Hannes in D. Hannes 2011 Rolling contact fatigue crack path prediction 
    ELSE IF(asp_shape .EQ. 8) THEN
        DO  J=1,NN
            DO  I=1,NX 
                ! Point
                bump=0
                IF ( 2*((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                ! Point
                bump2=0
                IF ( 2*((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw2 .LT. 1.0)THEN
                    bump2=asph2/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*2*PAI/(aspw2))+asph2/2;
                ENDIF
                H0=RAD(I,J)+W(I,J)-bump-bump2
                H(I,J)=H0  +H00
                
            ENDDO
        ENDDO
        
    ! Input parameters acc Holmes et al. Transient EHL point contact analysis 2003. Also the sam as Venner 1994 ball and cylinder But reshaped as a point defect
    ELSE IF(asp_shape .EQ. 13)THEN
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                Ab=asph
                Wb=aspw
                bump=Ab*10**(-10*((X(I)-position)**2+Y(J)**2)/Wb**2)*cos(2*PAI*sqrt((X(I)-position)**2+Y(J)**2)/Wb)                      !D_IN: /Setup/ -> Y(J)

                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0   +H00
                
            ENDDO
        ENDDO 
        
    ! The asperity shape used by D. Hannes in Rolling contact fatigue crack path prediction 2011 with variable length to width ratio   
    ELSE IF(asp_shape .EQ. 12 .OR. asp_shape .EQ. 122)THEN    
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                bump=0
                IF ( 2*((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN                                                !D_IN: /asp/ -> asp_ratio
                    bump=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0   +H00
                
            ENDDO
        ENDDO 
        
    ! The asperity shape used by D. Hannes in Rolling contact fatigue crack path prediction 2011 with variable length to width ratio   
    ELSE IF(asp_shape .EQ. 112)THEN   
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                bump=0
                IF ( 2*((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN                       ! We're at the asperity
                    asph1=asph/2*cos(((X(I)-position)**2+(asp_ratio*Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2
                    asph1=asph1*(2.0)**(0.5)                                                                    ! To make highest poing = asph and lowest = -asph
                    aa=X(I)-position
                   
                    IF ( 2*abs(aa)/aspw .LT. 0.5)THEN                                                           ! Central part
                    bump=asph1*sin(aa*2*PAI/(aspw));
                    
                    ELSE                                                                                        ! At the outer region, half the amplitude and the wavelength to get continiuos derivatives.   
                        IF (X(I) .LT. Position) then
                            bump=asph1/2*cos(2*aa*2*PAI/(aspw))-asph1/2                                         ! rise or lower the mean value depending on if in front or behind the center of the asperity
                        else
                            bump=-asph1/2*cos(2*aa*2*PAI/(aspw))+asph1/2 
                        endif
                    endif
                ENDIF
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ! The asperity shape used by D. Hannes in Rolling contact fatigue crack path prediction 2011 with changes length to width ratio for multiple defects   
    ELSE IF(asp_shape .EQ. 32)THEN  
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
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
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ! The asperity shape used by D. Hannes in Rolling contact fatigue crack path prediction 2011 with changes length to width ratio for multiple defects but with no offset. I.E half asperity half dent.       
    ELSE IF(asp_shape .EQ. 33)THEN  
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
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
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0   +H00
            ENDDO
        ENDDO 
        
    ! Using real surface measurements imported by the Asp_read routine   
    ELSEIF( asp_shape .ge. 200 .and. asp_shape .le. 209) then 
        if( t .LT. 1) then
            DO  J=1,NN,SS
                JJ=NYs-J+1
                DO  I=1,NX,SS 
                    H0=RAD(I,J)+W(I,J)
                    H(I,J)=H0   +H00
                    H(I,JJ)=H(I,J)
                ENDDO
            ENDDO
        else
            Call Asp_surf_calc(NX,NYs,NN,SS, t, H00)             !CALL: Only if asperity surface is imported Nx, Nys, NN = (Nys+1)/2, SS, t, H00
        endif
            
    ! The asperity shape used by D. Hannes in Rolling contact fatigue crack path prediction 2011 formulated as a indent    
    ELSE   
        DO  J=1,NN,SS
            DO  I=1,NX,SS 
                bump=0
                IF ( 2*((X(I)-position)**2+(Y(J))**2)**(0.5)/aspw .LT. 1.0)THEN
                    bump=asph/2*cos(((X(I)-position)**2+(Y(J))**2)**(0.5)*2*PAI/(aspw))+asph/2;
                ENDIF
                H0=RAD(I,J)+W(I,J)-bump
                H(I,J)=H0  +H00 
            ENDDO
        ENDDO   
    ENDIF
    
    !Mirror the surface around the symmetry line
    DO J=1,NN,SS
            JJ=NYs-J+1
            DO I=1,NX,SS
                H(I,JJ)=H(I,J)
            ENDDO
    ENDDO

    RETURN
    END
