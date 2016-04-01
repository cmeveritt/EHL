    !Subroutine for calculating how much efect pressure at one node has on the deformation of one other and storing it in a matrix  
	SUBROUTINE SUBAK
        implicit none
        COMMON      /Grid/NX,NY,X0,XE,DX     
        COMMON      /COMAK2D/AK(0:300,0:300)
        real        AK
        integer     NX,NY, I, J
        real        XP,XM,YP,YM, AA1, AA2, AA3, AA4, X, Y
        real        S, ALOG, ABS
        real        X0,XE,DX
        SAVE        /COMAK2D/

        
        
        S(X,Y)=X+SQRT(X**2+Y**2)
        DO 10 I=0,NX
            XP=I+0.5
            XM=I-0.5
            
            DO 10 J=0,I
                YP=J+0.5
                YM=J-0.5
                
                AA1=S(YP,XP)/S(YM,XP)
                AA2=S(XM,YM)/S(XP,YM)
                AA3=S(YM,XM)/S(YP,XM)
                AA4=S(XP,YP)/S(XM,YP)
                
                AK(I,J)=XP*ALOG(AA1)+YM*ALOG(AA2)+XM*ALOG(AA3)+YP*ALOG(AA4)
    10	AK(J,I)=AK(I,J)
        RETURN
    END