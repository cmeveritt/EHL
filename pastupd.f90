! Subrutine for storing the current parrameters to the past time step
    SUBROUTINE pastupd(a, NX, NN, SS)
        Implicit none
        include     'inc_Average_past.h'
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
        include     'inc_CurrentT.h'
	    include     'inc_CurrentT_met.h' 
	    include     'inc_CurrentRO.h'
        include     'inc_Past.h'
	    include     'inc_Past2.h'
        include     'inc_PastT.h'
	    include     'inc_PastT_met.h'
        include     'inc_PastRO.h'
        include     'inc_Therm_matr.h'
	    include     'inc_Therm_matr_past.h'
  
        integer     a                                     !a is 0 if t<0 and 1 otherwise
        integer     aa, i,j,i0, i1, j0, j1, NX, NN, SS
        SAVE        /past/
        SAVE        /PAST2/
        save        /pastT_met/
        SAVE        /PastT/ 
        SAVE        /PastRO/    
        SAVE        /Therm_matr_past/
        SAVE        /Average_past/
                
        CP_O_past=CP_O        !D_IN: /Term_matr/ -> CP_O !D_OUT:  CP_O_past -> /Therm_matr_past/
        
        if(a .eq. 0) then   !Time "independent" simulations
        Hpast2   = H                         !D_IN: /CurrentH/ -> H            !D_OUT: Hpast2 -> /Past2/
        Ppast2   = P                         !D_IN: /CurrentP/ -> P            !D_OUT: Ppast2 -> /Past2/
        Hpast    = H                                                           !D_OUT: Hpast -> /Past/
        ROpast   = RO                        !D_IN: /Current/ -> RO            !D_OUT: ROpast -> /Past/
        Ppast    = P                                                           !D_OUT: Ppast -> /Past/
        EPSxpast = EPSx                      !D_IN: /Current/ -> EPSx          !D_OUT: EPSxpast -> /Past/
        EPSypast = EPSy                      !D_IN: /Current/ -> EPSy          !D_OUT: EPSypast -> /Past/
        EDAxpast = EDAx                      !D_IN: /Current/ -> EDAx          !D_OUT: EDAxpast -> /Past/
        EDAypast = EDAy                      !D_IN: /Current/ -> EDAy          !D_OUT: EDAypast -> /Past/
        xipast   = xi                        !D_IN: /Current/ -> xi            !D_OUT: xipast -> /Past/
        Wpast    = W                         !D_IN: /Current/ -> W             !D_OUT: wpast -> /Past/
        dRodPpast = dRodP_mat                !D_IN: /CurrentRO/ -> dRodP_mat   !D_OUT: dRodPpast -> /PastRO/
        
        temp_map=temp_ma                     !D_IN: /CurrentT_met/ -> Temp_ma  !D_OUT: temp_mao -> /PastT_met/
        temp_mbp=temp_mb                     !D_IN: /CurrentT_met/ -> Temp_mb  !D_OUT: temp_mbp -> /PastT_met/
        tempp=temp                           !D_IN: /CurrentT/ -> Temp         !D_OUT: tempp -> /PastT/
        
        elseif(a .eq. 1) then !Truly time dependent simulations      !In this case use the values in the /Past/ block to update the /Past2/ block and the /Currentx/ to update the /Past/
        Hpast2   = Hpast                     
        Ppast2   = Ppast
        Hpast    = H
        ROpast   = RO
        Ppast    = P
        EPSxpast = EPSx
        EPSypast = EPSy
        EDAxpast = EDAx
        EDAypast = EDAy
        xipast   = xi
        Wpast    = W
        dRodPpast = dRodP_mat
        
        temp_map = temp_ma
        temp_mbp = temp_mb
        tempp = temp
        
        elseif(a .eq. 2)then                   ! a = 2 means update past temp only 
            
        temp_map = temp_ma
        temp_mbp = temp_mb
        tempp = temp
        
        else
            WRITE(4,*)'Bad input for pastupd. a =', a
            stop 'Bad input for pastupd.' 
        endif
        
        if( a .NE. 2) then
        do j=1,NN,ss
            J0=J-1*SS
            IF( J0 .LE. 0)  J0=1
            J1=J+1*SS
            IF( J1 .GE. NN) J1=NN-SS
            do i=1,NX,ss
                I0=I-1*SS
                IF( I0 .LE. 0)  I0=1
                I1=I+1*SS
                IF( I1 .GE. NX) I1=NX
                
                P_ave_past(i,j)=(P(i0,j)+P(i,j0)+4*P(i,j)+P(i1,j)+P(i,j1))/8                !D_OUT: P_ave_past -> /Average_past/
            enddo
        enddo
        endif
        

        RETURN
    END