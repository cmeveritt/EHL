! Subroutine for calculating the elastic deformation
	SUBROUTINE VI
        implicit none 
        COMMON      /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
        COMMON      /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referense
        COMMON      /COMAK2D/AK(0:300,0:300)
        COMMON      /HREEp/Emat, pg, Tg,YF                                       ! Parameters for HREE subrutine
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast  ! Past
        COMMON      /Setup/X,Y, P0                                       ! Parameters whihc do not change over time
        COMMON      /Residual/A5dxx,A5dx,A5dt                                   ! Parameters for the residual
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        real         X(1:300),Y(1:300), P0(1:300,1:300)
        real         A5dxx(1:300,1:300),A5dx(1:300,1:300),A5dt(1:300,1:300)
        real         Emat(1:300,1:300) ,pg(1:300,1:300), Tg(1:300,1:300), YF(1:300,1:300)
        real        AK
        integer     NX, NY, ref, K,L,I,J,IK,JL, tmeth
        real        X0, XE, DX, PAIAK, H0

        
        DO 40 I=1,NX
            DO 40 J=1,NY
                H0=0.0
              
                !Deformation from pressure inside the area
                DO 30 K=1,NX
                    IK=IABS(I-K)
                    DO 30 L=1,NY
                        JL=IABS(J-L)
                     
30              H0=H0+AK(IK,JL)*P(K,L)
                

        40      W(I,J)=H0*DX*PAIAK  ! See Eq 2.12 in the new book. The scaling factor between cylinder änd ball is pi/4. V_cyl=pi/4*V_ball due to the scaling R/a^2
        
        
        RETURN
    END