! Grid interpolation used to go from courser to finer grid. This is only used in the begning where the step size SS is decreased.  
    Subroutine Fine(H00,SSN, kvot, DT2, t, NYs, T40, Tcvot)  
    implicit none
    include     'inc_Com_H00.h'
    include     'inc_CurrentP.h'
    include     'inc_CurrentT.h'
    include     'inc_CurrentT_met.h' 
    include     'inc_Grid.h'
    include     'inc_Temp_param.h'
    ! Input
    integer     SS, NYs
    ! Calculations
    integer     I, J, SSN,k, NN_m, NY_m, NX_m
    real        H00
    ! Other       
    integer     t
    real        P3, kvot, DT2, T40, Tcvot
    ! Output
    save        /Com_H00/
    save        /CurrentP/
    save        /CurrentT/
    save        /CurrentT_met/ 
        
    SS=SSN*2
    IF(SSN .LE. 1)SSN=1
        
    ! Linear interpolation of P --------------------------------------------------------
    ! All exept the two far away edges
    DO J=1,NYs-SS,SS
        DO I=1,NX-SS,SS                                                                                          !D_IN: /Grid/ -> Nx                                                                                                  
            P(I,J)          =       P(I,J)                                                                       !D_OUT: P(I,J) -> /CurrentP/
            P(I+SSN,J)      =   0.5*(P(I,J)+P(I+SS,J))                                                           

            P(I,J+SSN)      =   0.5*(P(I,J)+P(I,J+SS))
            P(I+SSN,J+SSN)  =   0.25*(P(I,J)+P(I+SS,J)+P(I,J+SS)+P(I+SS,J+SS))
        ENDDO
    ENDDO
        
    ! The edge of J=NYs
    J=NYs
    DO I=1,NX-SS,SS
        P(I,J)          =        P(I,J)
        P(I+SSN,J)      =   0.5*(P(I,J)+P(I+SS,J))
    ENDDO
        
    !The edge of I=NX
    I=NX
    DO J=1,NYs-SS,SS
        P(I,J)          =        P(I,J)
        P(I,J+SSN)      =   0.5*(P(I,J)+P(I,J+SSN))
    ENDDO
        
    ! The end node already exists. 
    !P(NX,NYs)=P(NX,NYs)
        
    if(lub_temp==1) then                                                                                         !D_OUT: lub_temp -> /temp_param/
        ! Linear interpolation of the Temperature in the oil--------------------------------------------------------
        ! All exept the two far away edges
        DO J=1,NYs-SS,SS
            DO I=1,NX-SS,SS
                temp(I,J)          =       temp(I,J)                                                                 !D_OUT: temp -> /CurrentT/
                temp(I+SSN,J)      =   0.5*(temp(I,J)+temp(I+SS,J))

                temp(I,J+SSN)      =   0.5*(temp(I,J)+temp(I,J+SS))
                temp(I+SSN,J+SSN)  =   0.25*(temp(I,J)+temp(I+SS,J)+temp(I,J+SS)+temp(I+SS,J+SS))
                
            ENDDO
        ENDDO
        
        ! The edge of J=NYs
        J=NYs
        DO I=1,NX-SS,SS
                temp(I,J)          =        temp(I,J)
                temp(I+SSN,J)      =   0.5*(temp(I,J)+temp(I+SS,J))
        ENDDO
        
        !The edge of I=NX
        I=NX
        DO J=1,NYs-SS,SS
                temp(I,J)          =        temp(I,J)
                temp(I,J+SSN)      =   0.5*(temp(I,J)+temp(I,J+SSN))
        ENDDO
        
        ! Linear interpolation of the Temperature in the metal--------------------------------------------------------
        ! All exept the two far away edges
        NY_m=(nys+1)/2
        NN_m=(NY_m+1)/2
        NX_m=(NX+1)/2
        DO k=1,10
        DO J=1,NY_m,SS
            DO I=1,NX_m-SS,SS

                temp_ma(I,J,k)          =       temp_ma(I,J,k)                                                       !D_OUT: temp_ma -> /CurrentT_met/
                temp_ma(I+SSN,J,k)      =   0.5*(temp_ma(I,J,k)+temp_ma(I+SS,J,k))

                temp_ma(I,J+SSN,k)      =   0.5*(temp_ma(I,J,k)+temp_ma(I,J+SS,k))
                temp_ma(I+SSN,J+SSN,k)  =   0.25*(temp_ma(I,J,k)+temp_ma(I+SS,J,k)+temp_ma(I,J+SS,k)+temp_ma(I+SS,J+SS,k))
                
                temp_mb(I,J,k)          =       temp_mb(I,J,k)                                                       !D_OUT: temp_mb -> /CurrentT_met/
                temp_mb(I+SSN,J,k)      =   0.5*(temp_mb(I,J,k)+temp_mb(I+SS,J,k))

                temp_mb(I,J+SSN,k)      =   0.5*(temp_mb(I,J,k)+temp_mb(I,J+SS,k))
                temp_mb(I+SSN,J+SSN,k)  =   0.25*(temp_mb(I,J,k)+temp_mb(I+SS,J,k)+temp_mb(I,J+SS,k)+temp_mb(I+SS,J+SS,k))
                
            ENDDO
        ENDDO
        
        ! THe edge of J=NYs
        J=NY_m
        DO I=1,NX_m-SS,SS
                temp_ma(I,J,k)          =        temp_ma(I,J,k)
                temp_ma(I+SSN,J,k)      =   0.5*(temp_ma(I,J,k)+temp_ma(I+SS,J,k))
                
                temp_mb(I,J,k)          =        temp_mb(I,J,k)
                temp_mb(I+SSN,J,k)      =   0.5*(temp_mb(I,J,k)+temp_mb(I+SS,J,k))
        ENDDO
        
        !The edge of I=NX
        I=NX_m
        DO J=1,NY_m-SS,SS
                temp_ma(I,J,k)          =        temp_ma(I,J,k)
                temp_ma(I,J+SSN,k)      =   0.5*(temp_ma(I,J,k)+temp_ma(I,J+SSN,k))
                
                temp_mb(I,J,k)          =        temp_mb(I,J,k)
                temp_mb(I,J+SSN,k)      =   0.5*(temp_mb(I,J,k)+temp_mb(I,J+SSN,k))
        ENDDO
        ENDDO
    endif
        
    ! Update other parameters 
     CALL HREE(H00, t, SSN, NYs, T40, Tcvot, 0, 1, 0)                                                               !CALL: No comment
    H00past=H00                                                                                                     !D_OUT: /Com_H00/ -> H00past
    Call pastupd(0, NX, (NYs+1)/2, SS)                                                                              !CALL: No comment
        
    RETURN
    END
    