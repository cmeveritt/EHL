! Subroutine for estimating the error based on how much the pressure has been updated
    SUBROUTINE ERP(ER)
        implicit none
        COMMON      /Grid/NX,NY,X0,XE,DX     
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON      /Error/Pold
        integer     NX, NY, I, J, NN 
        real        er, ABS, sum
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real        Pold(1:300,1:300)
        real        X0,XE,DX
        SAVE        /Error/
        
        ER=0.0
        SUM=0.0
        NN=(NY+1)/2
        DO 10 I=1,NX
            DO 10 J=1,NN
                ER=ER+ABS(P(I,J)-Pold(I,J))
                SUM=SUM+P(I,J)
    10	CONTINUE
        ER=ER/SUM
        
        !Uppdate the past pressure. 
        Pold=P
        RETURN
    END