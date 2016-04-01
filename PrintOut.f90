! ******************************************************************************************************************************************** 
    SUBROUTINE Printout(H00,Tauc_real,asph_real, aspw_real,HM0,DT)
        Integer NX, NY
        integer Ntime, MK_stat, MK_time, KK
        integer term1, term2, term3, term4, ref, tmeth
        Real X0, XE, PAI
        Real RX, W0, Ua, Ub
        Real Z,  EDA0, Pref, alpha
        Real Elast1, Elast2, EE
        Real asph_real, aspw_real, ER_stat, ER_time
        Real DX, width
        Real asph, aspw,US, G0
        real ENDA
        real B, PH, H00, U, SRR
        real A1, A2, A3, W, ALFA, G, AHM, HM0, UTL, HM0r
        real tauc_real, tauc, tauS, taua, taua2,xilim
        real Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag
        real RA1,RA2, Tsize
        real kH, xH, gH, lH, L, M
        real    Uhd,Whd,Ghd
        COMMON /Visc/ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
        COMMON /NonNew/ tauc, taua, taua2, xilim                               ! Non Newtonian viscosity parameters
        COMMON /Yasutomi/Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag    ! Lubrication parameters acc Yasutomi
        COMMON /Rho/RA1,RA2                                             ! Density parameters
        COMMON /Grid/NX,NY,X0,XE,DX                                     ! Grid parameters
        COMMON /Itt/MK_stat, MK_time, ER_stat, ER_time, Ntime, KK   ! Itteration paramters
        COMMON /asp/asph,aspw                                           ! Apserity parameters
        COMMON /outp/ W0,EE,RX,US, Ua, B                                ! Output parameters
        COMMON /Tsize/Tsize                                             ! Size of time dep area
        COMMON /Holmes/kH, xH, gH, lH                                   ! Viscosity and denisty param acc Holems et al
        COMMON /Dimless/L,M                                             ! Diminsionless params
        COMMON /Method/term1,term2,term3,term4                          ! Controling the numerical method
        COMMON /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referense
        COMMON /Higginson/Uhd,Ghd,Whd
    
        ! Initaial data output -----------------------------------------------------------------------------------------------------------------
        OPEN(4,FILE='OUT.DAT',STATUS='UNKNOWN')
	    OPEN(8,FILE='FILM.DAT',STATUS='UNKNOWN')
	    OPEN(10,FILE='PRESSURE.DAT',STATUS='UNKNOWN')
        OPEN(2,FILE='EDA.DAT',STATUS='UNKNOWN')
        OPEN(3,FILE='EPS.DAT',STATUS='UNKNOWN')
        OPEN(5,FILE='Pside.DAT',STATUS='UNKNOWN')
        OPEN(11,FILE='A5dxx.DAT',STATUS='UNKNOWN')
        OPEN(12,FILE='A5dx.DAT',STATUS='UNKNOWN')
        OPEN(13,FILE='A5dt.DAT',STATUS='UNKNOWN')
        OPEN(14,FILE='Sig_x.DAT',STATUS='UNKNOWN')
        
        ! To Output file
        WRITE(4,*)'Grid infromation'
        WRITE(4,*)'NX,NY,X0,XE,DX,DT,KK'
	    WRITE(4,*) NX,NY,X0,XE,DX,DT,KK
        Write(4,*)'Ntime, MK_stat, MK_time, ER_stat, ER_time'
        Write(4,*) Ntime, MK_stat, MK_time, ER_stat, ER_time
        WRITE(4,*)''
    
        WRITE(4,*)'Geometry information'
        WRITE(4,*)'W0,EE,RX,US, Ua'
        WRITE(4,*) W0,EE,RX,US,Ua
        WRITE(4,*)'B, PH, H00'
	    WRITE(4,*) B, PH, H00
        WRITE(4,*)'HM0, HM0r'
	    WRITE(4,*) HM0, HM0r
        WRITE(4,*)''
        
        WRITE(4,*)'Density'
        WRITE(4,*)'RA1, RA2'
        WRITE(4,*) RA1, RA2
        
        WRITE(4,*)'Viscosity'
        WRITE(4,*) 'Z, EDA0, Pref, alpha,ENDA'
        WRITE(4,*) Z, EDA0, Pref, alpha,ENDA
        WRITE(4,*)''
        
        WRITE(4,*)'Solidification'
        WRITE(4,*)'Tauc_real, Taua, Taua2, xilim '
        WRITE(4,*) Tauc_real, Taua, Taua2, xilim
        WRITE(4,*)''
        
        WRITE(4,*)'Non Nwetonian shear thinning'
        WRITE(4,*)'Ta, Tg0, YA1, YA2, YB1, YB2'
        WRITE(4,*)'YC1, YC2, Yedag, Yalfag'
        WRITE(4,*)Ta, Tg0, YA1, YA2, YB1, YB2
        WRITE(4,*)YC1, YC2, Yedag, Yalfag
        WRITE(4,*)''
        
        WRITE(4,*)'Asperity'
        WRITE(4,*)'asph_real, aspw_real, asph, aspw'
        WRITE(4,*) asph_real, aspw_real, asph, aspw
        WRITE(4,*)'Tsize'
        WRITE(4,*) Tsize
        WRITE(4,*)''
        
        WRITE(4,*)'Holmes paramete4rs'
        WRITE(4,*)'lH, gH, xH, kH'
        WRITE(4,*) lH, gH, xH, kH
        WRITE(4,*)''
    
        WRITE(4,*)'Dimensionless params'
        WRITE(4,*)'L and M'
        WRITE(4,*) L, M
        WRITE(4,*)'U, G, W for line loads'
        WRITE(4,*) Uhd,Ghd,Whd
        WRITE(4,*)''
        
        WRITE(4,*)'Choise of method'
        WRITE(4,*)'The four numbers with the witght of the nodes'
        WRITE(4,*) term1, term2, term3, term4
        WRITE(4,*)'Te reference choosen, ref , PAIAK and tmeth'
        WRITE(4,*) ref, PAIAK, tmeth
        WRITE(4,*)''

        WRITE(*,*)'               Wait please'
    
        RETURN
    END
!