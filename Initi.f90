!Subroutine for initiating important parameters
    SUBROUTINE INITI
        implicit none
        COMMON      /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag    ! Lubrication parameters acc Yasutomi
        COMMON      /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
        COMMON      /Grid2/DX1,DX2,DX3,DX4
        COMMON      /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referense
        COMMON      /COMT/Temp
        COMMON      /COMH/RAD
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Past/Ppast,Hpast,ROpast,EPSxpast,EPSypast,EDAxpast,EDAypast,xipast,Wpast  ! Past
        COMMON      /Setup/X,Y, P0                                       ! Parameters whihc do not change over time
        COMMON      /AKparam/ AK00, AK10, AK20, AK30, BK00, BK10, BK20, BK30
        COMMON      /COMAK2D/AK(0:300,0:300)
        COMMON      /Error/Pold
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        real         X(1:300),Y(1:300), P0(1:300,1:300)
        real        RAD(1:300,1:300), Temp(1:300,1:300)
        real        Pold(1:300,1:300)
        real        AK,AK00, AK10, AK20, AK30, BK00, BK10, BK20, BK30
        Integer     NX, NY, ref, I, J, NN, JJ, tmeth
        real        Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
        real        X0,XE,Y0,DX
        real        DX1,DX2,DX3,DX4
        real        C,D, SQRT
        real        PAIAK
        SAVE        /Current/                                       ! Current timestep
        SAVE        /Past/                                          ! Past
        SAVE        /Setup/                                       ! Parameters whihc do not change over time
        SAVE        /AKparam/
        SAVE        /Grid/
        SAVE        /Grid2/
        SAVE        /Error/

     
        
        NN=(NY+1)/2
        DX=(XE-X0)/(NX-1.)
        Y0=-0.5*DX*NY+0.5*DX 
        
        DO I=1,NX
            X(I)=X0+(I-1)*DX
        ENDDO
        
        DO J=1,NY
            Y(J)=Y0+(J-1)*DX
        ENDDO
        
        IF(ref .EQ. 3 .OR. ref .EQ. 4)Then ! Ball
            DO  I=1,NX
                D=1.-X(I)*X(I)
                DO J=1,NN
                    C=D-Y(J)*Y(J)           ! Include or not depending on if Cylinder or Ball
                    IF(C.LE.0.0)P(I,J)=0.0
                    IF(C.GT.0.0)P(I,J)=SQRT(C)
                    RAD(I,J)=0.5*(X(I)*X(I)+Y(J)*Y(J))
                ENDDO
            ENDDO
        ELSE                    ! Cylinder
            DO  I=1,NX
                D=1.-X(I)*X(I)
                DO J=1,NN
                    C=D
                    IF(C.LE.0.0)P(I,J)=0.0
                    IF(C.GT.0.0)P(I,J)=SQRT(C)
                    RAD(I,J)=0.5*(X(I)*X(I))
                ENDDO
            ENDDO
        ENDIF

        DO I=1,NX
            DO J=NN+1,NY
                JJ=NY-J+1
                P(I,J)=P(I,JJ)
                RAD(I,J)=RAD(I,JJ)
            ENDDO
        ENDDO

        Ppast=P         ! Pressure at past timestep
        Pold=P          !Pressure at past itteration
        
        DO I=1,NX
            DO J=1,NY
                Temp(I,J)=Ta
            ENDDO
        ENDDO
        
        ! Parameters for the res subrutine in the ITER subrutine
        AK00=AK(0,0)
        AK10=AK(1,0)
        AK20=AK(2,0)
        AK30=AK(3,0)
        BK00=AK00-AK10
        BK10=AK10-0.25*(AK00+2.*AK(1,1)+AK(2,0))
        BK20=AK20-0.25*(AK10+2.*AK(2,1)+AK(3,0))
        BK30=AK30-0.25*(AK20+2.*AK(3,1)+AK(4,0))
       
        ! Grid parameters
        DX1=1./DX
        DX2=DX*DX
        DX3=1./DX2
        DX4=0.3*DX2
        
        RETURN
    END