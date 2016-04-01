! Subroutine for calculating the elastic deformation
	SUBROUTINE VI1D(NX,NY)
        implicit none
        real            P,V,AK1D
        integer     NX,NY, I,J,IJ
        real        DX, PAI1, C
        DIMENSION   P(NX,NY),V(NX)
        COMMON      /COMAK1D/AK1D(0:300)
        
        PAI1=0.318309886
        C=ALOG(DX)
        
        DO 10 I=1,NX
            V(I)=0.0
            DO 10 J=1,NX
                IJ=IABS(I-J)
10              V(I)=V(I)+(AK1D(IJ)+C)*DX*P(J,1)
        
        DO I=1,NX
            V(I)=-PAI1*V(I)
        ENDDO
        
        RETURN
    END