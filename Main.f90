! Thermal elastohydrodynamic calculation program
    ! This program was developed to handle thermal elastohydrodynamic contacts of varrying geometries
    ! To incoporate surface roughness and temperature generations
    ! It has been used to produce the data for among else the two folloeing publications in Tribology international:
    ! Contact fatigue initiation and tensile surface stresses at a point asperity which passes an elastohydrodynamic contact
    ! Surface initiation of rolling contact fatigue at asperities considering slip, shear limit and thermal elastohydrodynamic lubrication
    !
    ! The goal has been to generate a good enough code for the needed simulations. Therefore could the structure and the speed of the code be improved. 
    ! The code was developed by C.-M. Everitt in 2015-2020 based on the code publiched by P. Huang in 2013 in " Numerical calculations of lubrication: methods and programs". 
    ! The code was speeded up by Mr. G. Al Rheis in 2019, see the master thesis report "Parallelization of thermal elastohydrodynamic lubricated contacts simulation using OpenMP"
    !
	PROGRAM POINTEHL
#ifdef _OPENMP        
    Use omp_lib
    Use mkl_service
#endif
    Use IFPORT    
    implicit none                                                   ! The implicite non statement is used in all subroutines to ensure that all used parameters are initialized  http://www.obliquity.com/computer/fortran/common.html
    Integer L1NX,   L1NY,       DZ_method,  I_NUMBER                
    integer nr_N,   H00_method
    Real    PAI,    Y0,         H00,        HM0F
    Real    Elast1, Elast2      
    Real    asph_real,          aspw_real,  asph_real2,  aspw_real2
    Real    width,  alpha_bar,  Tauc_real
    real    U,      SRR,        Fe2,        Re 
    real    W,      ALFA,       G,          AHM,        HM0,    UTL ,       B_ref 
    real    Old_error,          old_C_true
    real    executionTime
    logical setAffinityStatus
    Common      /Error_conv/ Old_error, old_C_true
    include     'inc_Grid.h'
    include     'inc_Itt.h'
    include     'inc_Asp.h'
    include     'inc_Outp.h'
    include     'inc_Dimless.h'
    include     'inc_Method.h'
    include     'inc_Ref.h'
    include     'inc_Higginson.h'
    include     'inc_Com_H00.h'
    include     'inc_G0DT.h'
    include     'inc_Visc.h'
    include     'inc_Geom5.h'
    include     'inc_Temp_param.h'
    include     'inc_Glass_const.h'
    include     'inc_Hoglund_ball.h'
    include     'inc_Pres_ave_param.h'
    include     'inc_C_method.h'
    include     'inc_Numbers_of_bad.h'
    include     'inc_Rough_surf.h'
    include     'inc_Output_control.h'
    include     'inc_Single_step.h'
    include     'inc_Temp_conduction.h'
    include     'inc_Cores_used.h'
    
    !======================================================
    !------------------Affinity options--------------------
    !----------(Not needed in most use cases)--------------
    !======================================================
    
	!setAffinityStatus = SETENVQQ('KMP_HW_SUBSET=1t')                   !limit the threading to 1 thread per core.
    !setAffinityStatus = SETENVQQ('OMP_NESTED=true')                    !Turn on nested parallelism. Currently not used anywhere so its not needed.
    !if(.not. setAffinityStatus) print *, 'Failed to set the subset'
    
    !setAffinityStatus = SETENVQQ('KMP_AFFINITY=verbose,balanced')      !Print affinity report at the beginning of the program. Good to keep an eye on how threads are distributed if one is interested.
    !if(.not. setAffinityStatus) print *, 'Failed to set verbose'

    !======================================================
    !------------------------------------------------------
    !======================================================
    
   ! Open the file used for the general output regarding which time step the cod is at and the convergence 
    OPEN(4,FILE='OUT.DAT',STATUS='UNKNOWN')
    
    ! Open the fil containing the input data. Everything controlable is controled from the input file. 
    ! The data starts on line number two and is written on every third line so that the line above can be used for defining the name of the input fields
    open (unit = 7, file = "Input8_1.csv")
    
    ! First the general calculation methods and the geometry is read. 
    read (7,*)
    read (7,*) meth, tmeth, Geom, asp_shape, contact_alg, p_ave_param, shift_y, Multi_grid_param !D_IN: input8.csv
                                                                                                 !D_OUT: /Ref/
    PAI=3.14159265                      ! Only first 7 digits that counts since single precision
    
    ! The meth parameter defines which types of numerical method should be used for the spatial derivatives. 
    ! The backspace method defined for meth = 2 is recomended
    IF( meth .EQ. 1) THEN               ! 1st dx
        term1   =  0                    !D_OUT: /Method/
        term2   = -2
        term3   =  2
        term4   =  0
    ELSE IF( meth .EQ. 2)THEN           ! Backspace
        term1   =  1
        term2   = -4
        term3   =  3
        term4   =  0
    ELSE                                ! Central space Not working correctly
        term1   =  0
        term2   = -1
        term3   =  0
        term4   =  1
    ENDIF
    
    ! This line contains information for the grid size and timestepping
    read (7,*)
    read (7,*)
    read (7,*) NX, NY, X0, XE, F, DZ_method, Dw_meth !D_IN: input8.csv
                                                     !NX, Ny, X0, XE         !D_OUT:/Grid/
                                                     !F                      !D_OUT:/Ref/
                                                     !DZ_Method              !D_OUT:/Dz_com/
                                                     !Dw_meth                !D_OUT:/Itt/
    !Example values are 
    !   NX=150 which is the number of nodes in x = the rolling diredction RD
    !   NY=70  which is the number of nodes in y = the transverse diredction TD
    !   X0=-2.0 which is the dimensionless position of the inlet
    !   XE=1.5 which is the dimensionless position of the ooutlet
    !   F=2.0 which tells how many nodes the roughness hould move each timestep.
    !   DZ_method = 1 which increases the vertical distance between the nodes inside the metals as we move away from the lubricant
    !   DW_method = 1 which controles the method for finding load balance with help of the lubrication offset H00
    
    ! Grid refinment is uesd in the begining to speed up the process of achiving converged solution. 
    ! The number of nodes are therefore adjusted so that it will match with the courser grids. 
    ! A maximum of three refinements are used. The number of refinements are based upone the grid size. 
    L1NX = NX/8                         
    L1NY = NY/8
    L1NY = L1NY/2                       ! Ensuring odd number of nodes in Y-direction to ensure a single middel line. 
    L1NY = L1NY*2+1
    NX = (L1NX-1)*8+1             
    NY = (L1NY-1)*8+1
    DX = (XE-X0)/(NX-1.)                ! Space increment                       !DX             !D_OUT:/Grid/
    Y0 = -0.5*DX*NY+0.5*DX              ! Width in Y-direction
    width = DX*NY                       ! Width of model
    
    ! EHL parameters -------------------------------------------------------------
    !RX = 0.01 [m] is the radius of curvature in the rolling direction RD                        
    !W0 = 10.  [N/m] or [N] is the load parameter and the units depends thus upone the geometry. Its scaled with 1e5 later if modelling a cylinder          
    !Ua = 9.0  [m/s]    Velocity of the surface with the asperity  
    !Ub = 8.0  [m/s]    Velocity of the flat surface
    read (7,*)
    read (7,*)
    IF(geom .LE. 5)THEN
        read (7,*) RX, W0, Ua, Ub, lub_param, shear_thin, lub_temp, solid_material !D_IN: input8.csv
                                                                                   !Rx, W0, Ua, Ub             !D_OUT: /outp/
                                                                                   !Lub_param, shear_thin      !D_OUT: /ref/
                                                                                   !Lub_temp                   !D_OUT: /temp_param/
                                                                                   !solid_material             !D_OUT: /Glass_const/
        read (7,*)
        read (7,*)
        IF( geom .eq. 3) then
            read (7,*) b, by, Ph, RY    ! For geom 3 which is elliptical
                                        !D_IN: input8,csv
                                        !b, by, RY              !D_OUT: /outp/
                                        !Ph                     !D_OUT: /Visc/
            ph          = ph*1E9

        ELSE IF( Geom .EQ. 4)THEN ! Cylinder with given Ph 
            read (7,*) B_ref,  by, PH 
            Ph          =Ph*1e9       !Ph                       !D_OUT: /Visc/            
                                      !by                       !D_OUT: /outp/
            B_ref       =B_ref*1e-4
            
        ELSE IF( Geom .EQ. 5)THEN ! Cylinder with given Ph - only load change
            read (7,*) PH_new, by, PH       !by                 !D_OUT: /outp/
                                            !PH                 !D_OUT: /Visc/
                                            !PH_new             !D_OUT: /Geom5/
            PH_new      =PH_new*1e9
            Ph          =Ph*1e9
        else
            read(7,*)                   ! No extra geometry input
        endif
        
    ElSEIF(geom .EQ. 6)THEN             ! A try to model the experiments of a ball shot on a lubricated surface acc Höglund 1999 Influence of lubricant properties on elastohydrodynamic lubrication
        WRITE(4,*)'Warning! Not fully developed geoemtry module'
        read (7,*) RX, W0, V_ball, theta_ball, lub_param, shear_thin, lub_temp, solid_material  !D_IN: input8.csv
                                                                                                !RX, W0                     !D_OUT: /outp/
                                                                                                !Lub_param, shear_thin      !D_OUT: /ref/
                                                                                                !Lub_temp                   !D_OUT: /temp_param/
                                                                                                !solid_material             !D_OUT: /Glass_const/
                                                                                                !V_ball, theta_ball         !D_OUT: /Hoglund_ball/
        read (7,*)
        read (7,*)
        read (7,*) !b, by, Ph, RY
        
        theta_ball = theta_ball * PAI/180
        Vv_ball = V_ball * sin(theta_ball)
        Vh_ball = V_ball * cos(theta_ball)
        omega_ball = 0                  ! Initial rotaion speed
        Ub = Vh_ball
        Ua = 0
        
        M_ball = PAI * 4/3* RX**3 * 7850 ! Mass of ball
        I_ball = 2.0/5.0 * M_ball * Rx**2
        !Ub, Ua                                                         !D_OUT: /outp/
        !Vv_ball, Vh_ball, M_ball, I_ball, omega_ball                   !D_OUT: /Hoglund_ball/
    ELSE !Geom not 1-6
        WRITE(4,*)'Bad Geometric entry. Geom ='
	    WRITE(4,*) geom
        stop 'Bad Geometric entry'
    endif
    
    Um=(Ua+Ub)/2                        ! Mean velocity of fluid                            [m/s]
    !D_OUT: /outp/

    ! Asperity Data -------------------------------------------------------------
    !asph_real  = 1.5 [um]                 ! Real asperity height                             
    !aspw_real  = 50  [um]                 ! Real asperity radius 
    !asph_real2 = 1.5 [um]                 ! Real asperity height of second defect                           
    !aspw_real2 = 50  [um]                 ! Real asperity radius of second defect 
    !asph_ratio = 1                        ! Width ratio of first defect. Ratio between widt in RD and TD                           
    !Hminimim   = 80                       ! Dimensionless limit on lowes alloweble film thickness. If lower than this metal to metal contact is assumed. This is scaled with 1e-6 later
    !surf_scale = 1                        ! Sacles the surface if using measurements from a real surface
    read (7,*)
    read (7,*)
    read (7,*)  asph_real, aspw_real,asph_real2, aspw_real2, asp_ratio, Hminimum, surf_scale !D_IN: index8.csv
                                                                                             !asp_ratio, Hminimum           !D_OUT: /asp/
                                                                                             !surf_scale                    !D_OUT: /Rough_surf/
    ! Converting
    asph_real  = asph_real*1E-6
    aspw_real  = aspw_real*1E-6
    asph_real2 = asph_real2*1E-6
    aspw_real2 = aspw_real2*1E-6
    Hminimum   = Hminimum *1E-6
    
    
    ! Elastic modulus -------------------------------------------------------------
    ! Elast1=206E9                      ! Elastic modulus of first body                     [Pa]
    ! Elast2=206E9                      ! Elastic modulus of second body                    [Pa]
    ! EE=2.26E11                        ! Equvivalent E =E'                                 [N/m^2]
    read (7,*)
    read (7,*)
    read (7,*) Elast1, Elast2, EE                !D_IN: input8.csv
                                                 !EE                 !D_OUT:/outp/ 
    Elast1=Elast1*1E9
    Elast2=Elast2*1E9
    IF(EE .EQ. 0) THEN
        EE= 2/((1-0.3**2)/Elast1+(1-0.3**2)/Elast2)             ! If asperity case, ref2, Equvivalent Elastic modulus [Pa] Se Contact mech, KL Johnsson page 92. The book writes without a 2
    ElSE
        EE = EE*1E9                                             ! Equvivalent Elastic modulus [Pa]
    ENDIF
    
    
    ! Parameter calculations based upone the selescted geoemtry
    IF(Geom .EQ. 1)THEN         ! Cylinder
        W0    = W0*1E5
        nr_N  = 2*NX*NY/(XE-X0)                                 ! Number of nodes in the contact
        PAIAK = 0.15915489024                                   !D_OUT: /Ref/
        G0    = width*PAI/2                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
                                                                !G0                 !D_OUT: /G0DT/
        B     = SQRT(8*W0*RX/(PAI*EE))                          ! Contact halfwidth                     [m] (This comes from /outp/)
	    PH    = 2*W0/(PAI*B)                                    ! Hertzian pressure                     [N/m^2]
        
    ELSE IF(Geom .EQ. 2) THEN   ! Ball
        PAIAK = 0.2026423
        nr_N  = 2*NX*NY/(XE-X0)*1.0/(-Y0)                       ! Number of nodes in the contact
        
        G0    = 2.0/3.0*PAI                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan or Contact mech by KL Johnsson page 92.
        PH    = (3*W0*EE**2/(2*PAI**3*RX**2))**0.3333333        ! Hertzian pressure                     [N/m^2] se Contact mech by KL Johnsson page 93.
        B     = PAI*PH*RX/EE                                    ! Contact halfwidth                     [m] se Contact mech by KL Johnsson page 92.
        
    ELSE IF(Geom .EQ. 3) THEN   ! Elliptical Ball.              ! So far just based on EQs for a ball but with different radious
        PAIAK = 0.2026423
        nr_N  = 2*NX*NY/(XE-X0)*1.0/(-Y0)                       ! Number of nodes in the contact
        
        !G0=2.0/3.0*PAI                                         ! For elliptical contacts, this is defined later in Initi
        !Read as input                                          ! Hertzian pressure                     [N/m^2]
        !Read as input                                          ! Contact halfwidth                     [m]

        Re          = (Rx*Ry)**0.5
        if (b==0) then !No value for b given
            W0          = ph**3 * pai**3 * Re**2 / ( 6.0 * EE**2) !2.0/3.0*Ph*pai*b*by                             ! Applied load, see equation (4.27) page 92 in KL Johnsson Contact mech
            b           = sqrt(3*w0/(ph*2*pai*(Rx/Ry)**(-2.0/3.0)))
            by          = b*(Rx/Ry)**(-2.0/3.0)
        else
            W0=2.0/3.0*Ph*pai*b*by 
        endif
        
    Else IF(Geom .EQ. 4) THEN         ! Cylinder with given Ph
          
          ! PH    = 2*W0/(PAI*B)                                  ! Hertzian pressure   !Given in input.  
          PAIAK = 0.15915489024
          W0    = 2* PH**2 * PAI * RX / EE                        ! Applied load acc Johansson page 101 given that EE = 2*E'
          B     = SQRT(8*W0*RX/(PAI*EE))                          ! Contact halfwidth                     [m]
           
          XE    = XE * B_ref/B                                    ! Rescaling so that the real dx has the same length, not the dimensionless DX
          X0    = X0 * B_ref/B
          DX    = (XE-X0)/(NX-1.)                                 ! Space increment
          Y0    = -0.5*DX*NY+0.5*DX                               ! Width in Y-direction
          width = DX*NY                                           ! Width of model
          G0    = width*PAI/2                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
          nr_N  = 2*NX*NY/(XE-X0)                                 ! Number of nodes in the contact
          Geom  = 1                                               ! To continiue the calculations with a cylinder    
    Else IF(Geom .EQ. 5) THEN           ! Cylinder with given Ph - only load change
            
            W0    = 2* PH**2 * PAI * RX / EE
            nr_N  = 2*NX*NY/(XE-X0)                                 ! Number of nodes in the contact
            PAIAK = 0.15915489024
            G0    = width*PAI/2*(PH_new/PH)**2                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
            B     = SQRT(8*W0*RX/(PAI*EE))                          ! Contact halfwidth of Ph=2.3893356                    [m]
	        !PH    = 2*W0/(PAI*B)                                    ! Hertzian pressure                     [N/m^2]
   
    ELSE IF(Geom .EQ. 6) THEN   ! Höglund Ball Experimen
        W0 = 2 * Vv_ball * M_Ball / (60*1e-6)                   ! The nominal load is here calculated by asuming constant load over 60 us and the vertical spead compleatly changes direction. 
        
        PAIAK = 0.2026423
        nr_N  = 2*NX*NY/(XE-X0)*1.0/(-Y0)                       ! Number of nodes in the contact
        
        G0    = 2.0/3.0*PAI                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan or Contact mech by KL Johnsson page 92.
        PH    = (3*W0*EE**2/(2*PAI**3*RX**2))**0.3333333        ! Hertzian pressure                     [N/m^2] se Contact mech by KL Johnsson page 93.
        B     = PAI*PH*RX/EE  
   
    ELSE
        WRITE(4,*)'Bad Geometric entry. Geom ='
	    WRITE(4,*) geom
        stop 'Bad Geometric entry'
    ENDIF
    
    ! Accuracy paramters ------------------------------------------------------------
    ! MK_stat = 200      ! Maximum numbers of itterations for static solution
    ! MK_time = 200      ! Maximum numbers of itterations for timedep solution
    ! ER_stat = 1E-4     ! Maximum numbers of itterations for static solution
    ! ER_time = 2E-5     ! Maximum numbers of itterations for timedep solution
    ! KK=15              ! Number of internal itterations in each presure itteration
    
    read (7,*)
    read (7,*)
    read (7,*)     MK_stat, MK_time, ER_stat, ER_time, KK_stat, KK_time !D_IN: Input8.csv
                                                                        !D_OUT:/Itt/ 
    ER_stat = ER_stat*1e-3               
    ER_time = ER_time*1e-3
    
    old_C_true = 0.1        ! Setting initial reference value used in P_update
    bad_h_nr_node  =  10    ! Initiating for the contact routine
    bad_h_itt_node =  1     ! Initiating for the contact routine
    !D_OUT: /Number_of_bad/
    
    ! Parameters for restricting the relaxation parameter in P_update. 
    read (7,*)
    read (7,*)
    read (7,*)     C_meth, C_loc_, C_glob, C_min, C_max, umax_p  !D_IN: Input8.csv
                                                                !D_OUT: /C_method/
    
    ! Load ballance ------------------------------------------------------------------
    ! H00=-0.1                  ! Initial offset for lubrication height 
    ! HM0r                      ! Reduction factor for uppd of H00
    read (7,*)
    read (7,*)
    read (7,*)  H00, HM0f, H00_method, temp_r, Shear_band, Shear_xi !D_IN: Input8.csv
                                              !temp_r            !D_OUT: /temp_param/
    
    If ( geom .EQ. 6) Then ! Rescale h00 fore easier input
        H00=H00*Rx/b**2
    endif
    
    
    ! Choice of lubrication equation ------------------------------------------------
    Call Lubrication_def(HM0r,  Z, EDA0, alpha) !CALL no data in arguments                         ! ???Why pass arguments if they are not initialized here and if they are in a common block already?
    
    ! Read the rough surface if using that as asperities ------------------------------------------------
    if( asp_shape .ge. 200 .and. asp_shape .le. 209) then
        Call asp_read(NX,NY, Rx, b) !CALL D_OUT: Nx, Ny, Rx, b
    end if 
    
    
    ! Setup vertical space for metal nodes
    CALL DZ_setup(DZ_method, b, DX) !CALL D_OUT: DZ_method, b, Dx
    
    ! Pre calculations to generate dimensionless paramters ------------------------------------------------------------------------
    
    ! Geometry setup
    IF( Geom .LE. 5) then
        DT          = (Um*(XE-X0))/(Ntime*Ua)                       ! The time increment. T=t*Um/a. t_end=(XE-X0)*a/Ua => T_end=(XE-X0)Um/Ua   !Ntime is initialized in Lubrication_def
        !D_IN: /itt/ -> !Ntime
        !DT -> !D_OUT: /G0DT/
    ELSEIF( Geom .eq. 6) then
        Ntime=F                                                     !D_OUT:/itt/ 
        DT=60e-6/Ntime                                              !DT -> D_OUT: /G0DT/
    ENDIF
   
    
    DT0 = DT                                                        !D_OUT: /Hoglund_ball/
    
    SRR         = (Ua-Ub)/Um                                        ! Slide to Roll Ratio
    U2_over_Us  = Ua/(2.0*Um)                                       ! Alternative SRR number used by Venner
    !D_OUT: /outp/        
	
	AHM         = 1.0-EXP(-0.68*1.03)
	HM0         = 1.0!                                              ! Parameter for uppdating H0, the fim thicknes ofset value
    ! Alternative, 
    ! ALFA=Z*5.1E-9*A1
	! G=ALFA*EE
    ! HM0=3.63*(RX/B)**2*G**0.49*U**0.68*W**(-0.073)*AHM      
    HM0r        = HM0*HM0f                                          ! Reduced lubrication uppdation. ! ??? HM0r is updated without ever using it
    ENDA        = 12 *EDA0 *Um * RX**2 / ( B**3*PH)                 ! Epsilon in Reynolds equation. = 12.*(2*U)*(EE/PH)*(RX/B)**3  !EDA0 is defined in a common block and read from input8.csv
    !D_IN: /Visc/ -> EDA0            !ENDA -> !D_OUT: /Visc/ 
    if(lub_param .GE. 4 ) then    !It would make more sense if this came before the last two lines
        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx later since EDAO is dependent on the temperature
        !ENDA -> D_OUT: /Visc/
    endif
    
	UTL         = EDA0*Um*RX/(B*B*2.E7)
    !D_IN: /Visc/ -> EDA0
    H00past     = H00       !D_OUT: /Com_H00/  !I think this will be checked later when updating H00
    
    ! Moes and Bosma non-dimentional paramters
    ! L = alpha * EE * (2*EDA0*Um/(EE*RX))**0.25
    ! M = (W0/(EE*RX**2))*(EE*RX/(2*EDA0*Um))**0.75
    
    ! Higginson and Dowsson non-dim param from Spikes review from 2006 for line load
    Uhd         = Um*EDA0/(EE*RX)                                   ! Dimensionless speed parameter. Venner wrote U=EDA0*US/(2E*Rx) and US=2*Um.
    !D_OUT: /Higginson/
    Ghd         = alpha*EE                                          ! Dimensionless material parameter  !
    !alpha -> !D_IN:/Visc/       Ghd ->   !D_OUT: /Higginson/
    alpha_bar   = alpha*PH
    L           = Ghd*(2.0*Uhd)**0.25                               ! same as L above
    !D_OUT: /Dimless/

   IF(Geom .EQ. 1 .OR. Geom .EQ. 5)THEN
        !Acc C. H. Venner 1994 Transient Analysis of Surface line contact          
       Whd      = W0/(EE*RX) 
       M        = Whd/(2.0*Uhd)**(0.5) 
   ELSE
        !Acc C. H. Venner 1994 Nummerical simulations of a transverse ridge in a circular .......        
       Whd      = W0/(EE*Rx**2) 
       !D_OUT: /Higginson/
       M        = Whd/(2.0*Uhd)**(0.75)
       !D_OUT: /Dimless/
   ENDIF
   
    if( H00_method == 1) then           !This comes from Input8.csv
        ! H00 = -b^2/(2*Rx)*Rx/b^2  + Lubrication height from Dowson - Elastic deformation from KL Johnson    
        IF (geom == 1) then
            Fe2=0.3
            Re=Rx
            H00 = -0.5 + 2.65*Ghd**0.54*Uhd**0.7*Whd**(-0.13)*RX**2/(b**2) - 1.0 
                
        elseif (geom == 2) then
 
            Re=Rx
            H00 = -0.5 + 2.65*Ghd**0.54*Uhd**0.7*Whd**(-0.13)*RX**2/(b**2) - 1.0 !The elastic deflection is b^2/R
            
        elseif (geom == 3) then
            Fe2= 1.0-log10(sqrt((RX/RY))**0.5)/4 ! Error eltimate based on that Fe2(by/b=30)=0.65
            if( Fe2 .LT. 0.65) Fe2 = 0.65
            if( Fe2 .GT. 1) then
                Fe2= 1.0-log10(sqrt((RY/RX))**0.5)/4.0
                if( Fe2 .LT. 0.65) Fe2 = 0.65
                if( Fe2 .GT. 1) Fe2=1
            endif
            H00 = -0.5 + 2.65*Ghd**0.54*Uhd**0.7*Whd**(-0.13)*RX**2/(b**2) - ( 9 * W0**2 / ( 16 * EE**2*Re))**(1.0/3.0) *Fe2*RX/(b**2)+0.6

        endif
        !H00=  2.65*Ghd**0.54*Uhd**0.7*Whd**(-0.13)*RX**2/(b**2)
        !H00=- ( 9 * W0**2 / ( 16 * EE**2*Re))**(1.0/3.0) *Fe2*RX/(b**2) 
        
    endif
	
    ! Asperity  rescaling data
    asph        = asph_real*RX/B**2                                 ! Normalized height of asperity            [-]   
    aspw        = aspw_real/B                                       ! Normalized radius or wavelength of asperity            [-]
    asph2       = asph_real2*RX/B**2   
    aspw2       = aspw_real2/B
    !D_OUT: /asp/
    
    ! Generally, not all data is preinted ut. But if one wants to be able to restart the simulations from a specified timestep, all data is needed. 
    ! Thats controled by the following paramters
    read (7,*)
    read (7,*)
    read (7,*) collect_full_data, single_step_only, selected_step, use_fft
    
    p_update_initialized = .true.
    VI_collected = .true.
    H00_collected = .true.
    if(single_step_only) then
        collect_full_data = .false.
        p_update_initialized = .false.
        VI_collected = .false.
    endif
    if (collect_full_data) H00_collected = .false.

    ! Define if going to execute on one ore more processors. The code scales well up to at least 8 cores. 
#ifdef _OPENMP 
    read (7,*)
    read (7,*)
    read (7,*) use_multiple_cores, cores
    
    if(use_multiple_cores) then
        if (cores .le. 1) then       
            cores = OMP_GET_MAX_THREADS() * cores
            !print *, 'Cores : ', OMP_GET_MAX_THREADS()
            CALL OMP_SET_NUM_THREADS(cores)  ! To dynamically set the number of threads this should be called once with omp_get_max_threads
            CALL MKL_SET_NUM_THREADS(cores)
        else
            CALL OMP_SET_NUM_THREADS(cores)
            CALL MKL_SET_NUM_THREADS(cores)
        endif
    endif    
#endif   
    !--------------------------------
    
    ! Printout the data
    CALL Printout(H00,Tauc_real,asph_real, aspw_real,HM0, alpha_bar) !CALL D_OUT: H00, asph_real, aspw_real,HM0, alpha_bar  ! !!! Tauc_real wasn't initialized in main for that reason its not in the data out list

    ! START simulation ---------------------------------------------------------------------------------------------------------------------------
    ! Only passing dimentionless paramters and setup parameters
    CALL TimeCalc(H00, aspw, executionTime) !CALL D_OUT: H00, aspw                                   !CALL: this handles getting the solutions for all time steps in the simulation

    OPEN(20,FILE='1_Finished_Wihoooo.DAT',STATUS='UNKNOWN')                           


#ifdef _OPENMP
    executionTime = OMP_GET_WTIME() - executionTime    !End measuring the execution time. Measuring starts from within the timecalc subroutine.
#endif

    ! To Output file
    WRITE(20,*)'This is awesome'
    WRITE(20,*)'Execution Time: ', executionTime
    Print *, 'Execution Time: ', executionTime

	STOP

    END   
