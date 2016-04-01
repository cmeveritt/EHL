! Subroutine for defining the output data
	SUBROUTINE OUTPUT
        implicit none
        COMMON      /Grid/NX,NY,X0,XE,DX                                                    ! Grid parameters
        COMMON      /HREEp/Emat, pg, Tg,YF                                                   ! Parameters for HREE subrutine
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                          ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast  ! Past
        COMMON      /Setup/X,Y, P0                                                          ! Parameters whihc do not change over time
        COMMON      /Residual/A5dxx,A5dx,A5dt                                               ! Parameters for the residual
        COMMON      /Stress/sig_x, sig_y, sig_z                   
        real         sig_x(300), sig_y(300), sig_z(300)
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        real         X(1:300),Y(1:300), P0(1:300,1:300)
        real         A5dxx(1:300,1:300),A5dx(1:300,1:300),A5dt(1:300,1:300)
        real         Emat(1:300,1:300) ,pg(1:300,1:300), Tg(1:300,1:300), YF(1:300,1:300)
        integer     NX,NY, NN, I, J
        real        X0, XE, DX, A

        
        NN=(NY+1)/2
        A=0.0
        
        WRITE(8,110)A,(Y(I),I=1,NY)
        DO I=1,NX
            WRITE(8,110)X(I),(H(I,J),J=1,NY)
        ENDDO
        
        WRITE(10,110)A,(Y(I),I=1,NY)
        DO I=1,NX
            WRITE(10,110)X(I),(P(I,J),J=1,NY)        !(10 is the file number, 110 is the format number)
        ENDDO
        
        WRITE(2,110)A,(Y(I),I=1,NY)
        DO I=1,NX
        WRITE(2,110)X(I),(EDAx(I,J),J=1,NY)        
        ENDDO
        
        WRITE(3,110)A,(Y(I),I=1,NY)
        DO I=1,NX
        WRITE(3,110)X(I),(EPSx(I,J),J=1,NY)        
        ENDDO
        
        WRITE(11,110)A,(Y(I),I=1,NY)
        DO I=1,NX
        WRITE(11,110)X(I),(A5dxx(I,J),J=1,NY)       
        ENDDO
        
        WRITE(12,110)A,(Y(I),I=1,NY)
        DO I=1,NX
        WRITE(12,110)X(I),(A5dx(I,J),J=1,NY)       
        ENDDO
        
        WRITE(13,110)A,(Y(I),I=1,NY)
        DO I=1,NX
        WRITE(13,110)X(I),(A5dt(I,J),J=1,NY)       
        ENDDO
        
        WRITE(14,110)A,(sig_x(I),I=1,NX)
        
    110 FORMAT(2001(E12.6,1X))
        RETURN
    END