!Subroutine for calculating how much efect pressure at one node has on the deformation of one other and storing it in a matrix in 1D
	SUBROUTINE SUBAK1D(N)
        implicit none
        real        AK1D
        integer I,N
        real    ABS, ALOG
        COMMON /COMAK1D/AK1D(0:300)
       
        DO 10 I=0,N
10           AK1D(I)=(I+0.5)*(ALOG(ABS(I+0.5))-1.0)-(I-0.5)*(ALOG(ABS(I-0.5))-1.0)
        RETURN
    END