! Subrutine for uppdating the parameters of the past timestep
    SUBROUTINE pastupd
        Implicit none   
        COMMON      /Current/P,     H,     RO,     EPSx,     EPSy,     EDAx,     EDAy,     xi,     W,    Wside                      ! Current timestep
        COMMON      /Past/   Ppast, Hpast, ROpast, EPSxpast, EPSypast, EDAxpast, EDAypast, xipast, Wpast  ! Past                           
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         Ppast(1:300,1:300),Hpast(1:300,1:300),ROpast(1:300,1:300),EPSxpast(1:300,1:300),EPSypast(1:300,1:300),EDAxpast(1:300,1:300),EDAypast(1:300,1:300),xipast(1:300,1:300),Wpast(1:300,1:300)  ! Past
        SAVE        /past/
        
        
        Hpast=H
        ROpast=RO
        Ppast=P
        EPSxpast=EPSx
        EPSypast=EPSy
        EDAxpast=EDAx
        EDAypast=EDAy
        xipast=xi
        Ppast=P
        Wpast=W
        
        RETURN
    END