    !Subroutine for calculating the efect a pressure at one node has on the deformation of one other. The data is stored in a matrix  where the indicies indicates the distance
	SUBROUTINE SUBAK
        implicit none
        include     'inc_Grid.h'
        include     'inc_COMAK2D.h'
	    include     'inc_COMAKline.h'
        integer     I, J, NX2, NX4, NX8
        real        XP,XM,YP,YM, AA1, AA2, AA3, AA4, X, Y
        real        S, ALOG, ABS
        SAVE        /COMAK2D/
        SAVE        /COMAKline/

        ! See Contact mechanics by KL Johnsson page 54 for more theory behind int. Its based on equations from Love (1929)
        
        S(X,Y)=X+SQRT(X**2+Y**2)        ! ??? What is happening here? Is this a transformation?
        DO I=0,NX                                                                                                                  !D_IN: /Grid/ -> NX
            XP=I+0.5
            XM=I-0.5
            
            DO J=0,I
                YP=J+0.5
                YM=J-0.5
                
                AA1=S(YP,XP)/S(YM,XP)
                AA2=S(XM,YM)/S(XP,YM)
                AA3=S(YM,XM)/S(YP,XM)
                AA4=S(XP,YP)/S(XM,YP)
                
                AK(I,J)=XP*ALOG(AA1)+YM*ALOG(AA2)+XM*ALOG(AA3)+YP*ALOG(AA4)                                                         !D_OUT: AK -> /COMAK2D/
                AK(J,I)=AK(I,J) !This might cause really bad cache utilization maybe can be re done as shown in figure 5.15 of Using OpenMP
            ENDDO
        ENDDO
        
        sum_Line_AK8 = AK(0,0)                                                                                                      !D_OUT: sum_Line_AK8 -> /COMAKline/
        
        NX8=(NX-1)/8+1
        NX4=(NX-1)/4+1
        NX2=(NX-1)/2+1
        
        DO J=1,NX8           !Since we calculate the sidepressure on an area which is NX in both directions
            sum_Line_AK8 = sum_Line_AK8 + 2*AK(0,J) !One for each side
        ENDDO
        
        sum_Line_AK4=sum_Line_AK8                                                                                                   !D_OUT: sum_line_AK4 -> /COMAKline/
        
        DO J=NX8,NX4         
            sum_Line_AK4 = sum_Line_AK4 + 2*AK(0,J) !One for each side
        ENDDO
        
        sum_Line_AK2=sum_Line_AK4
        DO J=NX4,NX2        
            sum_Line_AK2 = sum_Line_AK2 + 2*AK(0,J) !One for each side                                                              !D_OUT: sum_line_AK2 -> /COMAKline/
        ENDDO
        
        sum_Line_AK=sum_Line_AK2 
        DO J=NX2,NX         
            sum_Line_AK = sum_Line_AK + 2*AK(0,J) !One for each side                                                                !D_OUT: sum_line_AK -> /COMAKline/
        ENDDO
                
        RETURN
    END