! Subrutine for updating the film  the viscocity and the density based on the given pressure 
    ! Based upon the Yasutomi et al. (1984) ewquations as formulated by Ehret et al. 1997 On lubricant transport conditions in elastohydrodynamic conjunctions
    ! This subroutine does not work properly
    ! Should add that if contact has occured the viscosity gets really high. Othervise the code will have a hard time converging
	SUBROUTINE Ehret_shear_thin(SS, NYs, T40, Tcvot, NN)
        implicit none 
        include     'inc_Contact_mat.h'
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
        include     'inc_Grid.h'
        include     'inc_NonNew.h'
        include     'inc_NonNew2.h'
        include     'inc_Outp.h'
        include     'inc_Ref.h'
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'

        ! Input
        integer     SS, NN
        real        dpdy, dpdx
        ! Calculations
        real        PAI
        integer     I,J, itt
        real        EDA2, EDA_cont
        real        taux, tauxp, tauy, Lub_h
        ! Other
        Integer     ink, NYs
        real        r, r_prim, taup
        real        temp,  Tau, EDA1, tau1, tau2
        real        av_p, T40, Tcvot
        real        Us, Ul, U_all(2)
        ! Output
        SAVE        /Current/                                                               
        DATA PAI/3.14159265/

        WRITE(4,*) 'Warning, the Ehret shear thinning subroutine does not work properly' 
        
        EDA_cont=Yedag                                                      !D_IN: /Yasutomi/ -> Yedag
        
        ! Yasutomi lubrication
        DO J=1,NN,SS
            DO I=1,NX,SS
                EDA1=EDAx(I,J)                                              !D_IN: /Current/ -> EDAx(I,J)

                tau1=EDA1/H(I,J)/tauc                                       !D_IN: /CurrentH/ -> H(I,J); /NonNew/ -> tauc
                if( tau1 .GT. 10000.0) THEN
                    tau=LOG(tau1+abs(tau1) )
                ELSE
                    tau2=tau1*tau1+1.0
                    tau=LOG(tau1+SQRT(tau2) )                               ! Absolute value since does not matter which direction the sliding is. tauc includes abs(a**2/(R*(Ua-Ub)))
                endif
                
        
                if( lub_param .EQ. 56 .and. EDA1 .LT. Yedag)then            !D_IN: /Ref/ -> Lub_param
                    if(i .eq. 1 .or. i .eq. NX) then
                        dpdx=0.0
                    else
                        dpdx    = (P(i+1*SS,j)-P(i-1*SS,j))/(2*DX*SS)       !D_IN: /Grid/ -> Dx; /CurrentP/ -> P(I,J)
                    endif
                    
                    if( j .eq. 1 .or. j .eq. nn) then
                        dpdy=0.0
                    else
                        dpdy    = (P(i,j+1*SS)-P(i,j-1*SS))/(2*DX*SS)
                    endif
                    
                    Lub_h   = H(i,j)*B**2/Rx/2                              !half the lubrication height since this will give me the highes shear value           !D_IN: /outp/ -> B, Rx
                    taux    = EDA1*Uslip/(Lub_h*2)+Lub_h*dpdx               ! taux                !D_IN: /NonNew2/ -> Uslip
                    tauy    = Lub_h*dpdy
                    tau     = sqrt(taux**2+tauy**2)                         ! equvivalent shear 
                    r       = 1
                    itt     = 1
                    taup    = tau
                    tauxp   = taux
                    !tau=tau*tauc_real      !Tau is scaled with tuc_limit to simplify convergence
                    do while (abs(r) .GT. 0.0001 .and. itt .LT. 20)
                        EDA2    = EDA1*tau/tauc_real/sinh(tau/tauc_real)                      !D_IN: /NonNew2/ -> tauc_real
                        if (tau .LE. 0.0) EDA2=EDA1
                        
                        taux    = EDA2*Uslip/(Lub_h*2)+Lub_h*dpdx           ! taux
                        tau     = 0.5*sqrt(taux**2+tauy**2)+0.5*taup                 ! equvivalent shear 
                        
                        tau2    = sqrt((1.01*taux)**2+tauy**2)              !Need a better estimate of the derivative. a numerical derivative shoud work fine. 
                        EDA2    = EDA1*tau2/tauc_real/sinh(tau2/tauc_real)
                        r       = (tau-taup)/tau                            !and a better residual function for whcih the derivativa can be defined

                        tauxp   = taux
                        taup    = tau
                        itt     = itt+1
                        if(itt .eq. 18)then
                            EDA2=H(i,j)
                            EDA1=P(i,j)
                        endif
                    enddo
                    tau=tau/tauc_real                                       !Rescaling to effective tau
                elseif(lub_param .EQ. 56 .and. EDA1 .GE. Yedag)then
                    tau=0.0                                                 !No reduction if solidifed
                elseif(lub_param .EQ. 57 .or. lub_param .EQ. 58) then       ! Implementation based on eq 2.6 tau=2eda*gamma7f(tau*) and that df/dt means df/dt
                    if ( uslip .eq. 0.0) then
                        tau=0.0
                    else
                        tau=asinh(EDA1/tauc_real * Uslip/(H(i,j)*B**2/Rx))  ! Scaled with tauc_real
                    endif
                elseif( lub_param .eq. 52) then                             ! Implementation based on that df/dt means df/dt* derived from the equations of EDA
                    if ( uslip .eq. 0.0) then
                        tau=0.0
                    else
                        r=1
                        itt=1
                        taup=0
                        do while (abs(r) .GT. 0.001 .and. itt .LT. 20 .and. taup .NE. tau)
                            r       = tau * tauc_real * cosh(tau) - EDA1 * Uslip/(H(i,j)*B**2/Rx)
                            r_prim  = tau * sinh(tau) + cosh(tau)
                            taup    = tau
                            tau     = tau -r / r_prim / tauc_real
                            itt     = itt +1
                        enddo
                    endif
                else
                    tau1=EDA1/H(I,J)/tauc
                    if( tau1 .GT. 10000.0) THEN
                        tau=LOG(tau1+abs(tau1) )
                    ELSE
                        tau2=tau1*tau1+1.0
                        tau=LOG(tau1+SQRT(tau2) )                           ! Absolute value since does not matter which direction the sliding is. tauc includes abs(a**2/(R*(Ua-Ub)))
                    endif
                endif
                
                IF (tau .EQ. 0.0) Then
                    EDAx(I,J)=EDA1                                          ! No Non-newtonian reduction.           !D_OUT: EDAx -> /Current/
                    EDAy(I,J)=EDA1                                                                                  !D_OUT: EDAy -> /Current/
                ELSEif (lub_param .EQ. 57 .or. lub_param .eq. 59)then
                    ! Now with trial equations where added an extra tau*tauc_real
                    EDAx(I,J)=EDA1/(sinh(tau)/tau*(1-tau)+tau*cosh(tau))    ! Rewritten to get better numbers. !EDA1/(f_tau+tau*df_tau/dtau)
                    EDAy(I,J)=EDA1*tau/(SINH(tau))                          ! f_tau
                ELSE
                    EDAx(I,J)=EDA1/(COSH(tau))                              ! Rewritten to get better numbers. !EDA1/(f_tau+tau*df_tau/dtau*)
                    EDAy(I,J)=EDA1*tau/(SINH(tau))                          ! f_tau
                ENDIF
                
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity       !D_IN: /Contact_mat/ -> contact
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))                                                                        !D_IN: /Current/ -> RO      !D_OUT: EPSx -> /Current/
                
                IF(contact(I,J) .EQ. 2 .and. EDAy(I,J) .LT. EDA_cont) EDAy(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSy(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))                                                                                                    !D_OUT: EPSy -> /Current/
                

                ! Solidification calculations
                ! Two different models were added but did not yet yeild satisfactory results. 
                if( xi_param .EQ. 1)Then                                !CM.Everitts arctan P model                                 !D_IN: /NonNew/ -> xi_param 
                    xi(I,J)     = atan(taua*P(I,J))*Uslip/(UM*PAI)                                                                  !D_IN: /NonNew/ -> taua; /outp/ -> Um        !D_OUT: xi(I,J) -> /Current/
                    IF( xi(I,J) .GT. 0.99 )then                         ! Double check since added after
                         xi(I,J) = 0.99
                    elseif( xi(I,j) .LT. 0)Then
                        xi(i,j) =0                                      
                    endif
                    
                Elseif( xi_param .EQ. 2)Then                            !CM.Everitts arctan dp/dx model
                    if(i .EQ. 1 .or. i .EQ. Nx) then
                        dpdx    = 0.0
                    else
                        dpdx    =(P(i+1*SS,j)-P(i-1*SS,j))/(2*DX*SS)
                    endif              
                    xi(I,J)     =atan(taua*dpdx)*Uslip/(UM*PAI) 
                    
                else                                                      !Accotding to P Ehret
                    xi(I,J)=taua*P(I,J)+taua2*P(I,J)**2 
                    IF(xi(I,J) .GT. xilim) xi(I,J)=xilim                  !D_IN: /NonNew/ -> xilim 
                    IF(xi(I,J) .LT. 0.0) xi(I,J)=0.0
                endif

        ENDDO
    ENDDO
    
    RETURN
    END
