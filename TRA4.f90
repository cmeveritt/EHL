 ! Subrutine for the pressureupdates
	SUBROUTINE TRA4(N)
        implicit none 
        COMMON              /CTRA4/ D, A
        COMMON              /TRA/  B
        real            D(1:300), A(1:5,1:300), B(1:3,1:300)
        integer         in1, N, INN, I 
        real            C
       
        C=1./A(3,N)
        B(1,N)=-A(1,N)*C
        B(2,N)=-A(2,N)*C
        B(3,N)=A(5,N)*C
        DO 10 I=1,N-2
            INN=N-I
            IN1=INN+1
            C=1./(A(3,INN)+A(4,INN)*B(2,IN1))
            B(1,INN)=-A(1,INN)*C
            B(2,INN)=-(A(2,INN)+A(4,INN)*B(1,IN1))*C
    10	    B(3,INN)=(A(5,INN)-A(4,INN)*B(3,IN1))*C
        D(1)=0.0
        D(2)=B(3,2)
        DO 20 I=3,N
    20	    D(I)=B(1,I)*D(I-2)+B(2,I)*D(I-1)+B(3,I)
        RETURN
    END