! Subroutine for estimating the error based on how much the pressure has been updated
    SUBROUTINE ERP(ER,NYs,max_diff,SS)
        implicit none
        include     'inc_Grid.h'
        include     'inc_CurrentP.h'
        include     'inc_Error.h'
        integer      I, J, NN, NYs,SS 
        real        er, ABS, sum
        real        diff, max_diff
        SAVE        /Error/
        
        ER=0.0
        SUM=0.0
        max_diff=0.0;
        NN=(NYs+1)/2
        
        DO 10 J=1,NN,SS
            DO 10 I=1,NX,SS                                                          !D_IN: /Grid/ -> Nx
                diff        = ABS(P(I,J)-Pold(I,J))                                  !D_IN: /CurrentP/ -> P(I,J); /Error/ -> Pold
                ER          = ER+diff
                max_diff    = max(max_diff, diff)
                SUM         = SUM+P(I,J)
    10	CONTINUE
        ER=ER/SUM

        !Uppdate the past pressure. 
        Pold=P                                                                       !D_OUT: Pold -> /Error/
        RETURN
    END