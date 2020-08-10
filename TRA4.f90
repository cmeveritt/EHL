 ! Subrutine for the pressureupdates for a given line in RD. The updates are given in the vector D(I)
    ! Based on the residual of the Reynolds equation on the deriatives of reynolds equation based upone the pressure of different nodes. 
	SUBROUTINE TRA4(N,SS)
        implicit none 
        include     'inc_CTRA4.h'
        ! Calculation parameters
        real            B(1:3,1:601)
        integer         in1, N, INN, I, SS 
        real            C
        ! Output
        SAVE            /CTRA4/
   
        C=1./A(3,N)                                                                    !D_IN: /CTRA4/ -> A
        B(1,N)=-A(1,N)*C
        B(2,N)=-A(2,N)*C
        B(3,N)=A(5,N)*C
        
        DO I=1,N-2*SS,SS
            INN=N-I+1-1*SS
            IN1=INN+1*SS
            C=1./(A(3,INN)+A(4,INN)*B(2,IN1))
            B(1,INN)=-A(1,INN)*C
            B(2,INN)=-(A(2,INN)+A(4,INN)*B(1,IN1))*C
    	    B(3,INN)=(A(5,INN)-A(4,INN)*B(3,IN1))*C
        ENDDO
        
        D(1)=0.0                                                                       !D_OUT: D -> /CTRA4/
        D(1+SS)=B(3,1+1*SS)
        DO I=1+2*SS,N,SS
    	    D(I)=B(1,I)*D(I-2*SS)+B(2,I)*D(I-1*SS)+B(3,I)
                    if ( isnan(D(I) ) ) then 
                        WRITE(4,*)'check 22225 Bad D=',D(I),' at i=', i
                        WRITE(4,*)'B(1,i)=', B(1,i)
                        WRITE(4,*)'B(2,i)=', B(2,i)
                        WRITE(4,*)'B(3,i)=', B(3,i)
                        
                        WRITE(4,*)'A(1,i)=', A(1,i)
                        WRITE(4,*)'A(2,i)=', A(2,i)
                        WRITE(4,*)'A(3,i)=', A(3,i)
                        WRITE(4,*)'A(4,i)=', A(4,i)
                        WRITE(4,*)'A(5,i)=', A(5,i)
                        
                        if( I .LT. N) then
                            WRITE(4,*)'A(1,i+ss)=', A(1,i+ss)
                            WRITE(4,*)'A(2,i+ss)=', A(2,i+ss)
                            WRITE(4,*)'A(3,i+ss)=', A(3,i+ss)
                            WRITE(4,*)'A(4,i+ss)=', A(4,i+ss)
                            WRITE(4,*)'A(5,i+ss)=', A(5,i+ss)
                        endif
                        
                        if( I .gt. 1)then
                            WRITE(4,*)'A(1,i-ss)=', A(1,i-ss)
                            WRITE(4,*)'A(2,i-ss)=', A(2,i-ss)
                            WRITE(4,*)'A(3,i-ss)=', A(3,i-ss)
                            WRITE(4,*)'A(4,i-ss)=', A(4,i-ss)
                            WRITE(4,*)'A(5,i-ss)=', A(5,i-ss)
                        endif
                    endif
            ENDDO
        
        RETURN
    END