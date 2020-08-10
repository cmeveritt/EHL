! A subroutine to calculate the increments of Epsilon based on apressure, p, change
    ! EDA = epsilon = density * film thickness^3 / (viscosity * lambda) term in the Reynolds equation
    ! Where lambda = 12 * average spped * radius^2 / contact halfwidth^3 * Hertzian pressure)
    ! d(eps(i,j))/d(P(k,l)) = incre
    Subroutine Epsilon_derivative(i,j,k,l, K00, K10, K20, K30, Tcvot, ss, incre, dRodP, d2ROdP2)
    implicit none
    include     'inc_COMAK2D.h'
    include     'inc_Contact_mat.h'
	include     'inc_Current.h'
	include     'inc_CurrentH.h'
	include     'inc_CurrentP.h'
    include     'inc_CurrentT.h'
    include     'inc_CurrentRO.h'
    include     'inc_Holmes.h'
    include     'inc_Ref.h'
    include     'inc_Rho.h'
    include     'inc_RLarsson.h'
    include     'inc_Visc.h'
    include     'inc_Yasutomi.h'
    include     'inc_Y_liu.h'
    !Input
    integer     i,j,k,l, ss
    real        K00, K10, K20, K30
    real        Tcvot
    ! Calculations
    real        eda01,eda02,eda03, EDA22, eda33
    real        EDA
    real        dEDAdP, dEDAdtemp, dRodtemp, dhdp, dTemp_dp, dYfdp
    real        h_cube, drodp1, dRodP0
    real        temp_y, YF, TG, T40, Pg
    real        dZdtemp, dEDA0dtemp, dEDAdtemp1, dEDAdtemp2, dEDAdtemp3, EpsT
    real        EDA04, Delta_s, dEDA03dtemp, dEDA04dtemp
    !Output     
    real        incre, dRodP, d2ROdP2
    save        /CurrentRO/
     
    ! Initialize parameters to zero, which is true if not updated later. 
    dTemp_dp    = 0
    dEDAdP      = 0
    dEDAdtemp   = 0
    dRodP       = 0
    dRodtemp    = 0
    incre       = 0
    h_cube      = H(i,j)**3                                                                                             !D_IN: /CurrentH/ -> H(I,J)                                                                                                                                  
    
    ! If along the x-line
    if( (     i==k+1*ss .or. i==k-1*ss) .and. j==l) then
        dhdp=K10
    else if( (i==k+2*ss .or. i==k-2*ss) .and. j==l) then
        dhdp=K20
    else if( (i==k+3*ss .or. i==k-3*ss) .and. j==l) then
        dhdp=K30  
    !If along the y-line
    elseif( (j==l+1*ss .or. j==l-1*ss) .and. i==k) then
        dhdp=K10
    else if( (j==l+2*ss .or. i==k-2*ss) .and. i==k) then
        dhdp=K20
    else if( (j==l+3*ss .or. i==k-3*ss) .and. i==k) then
        dhdp=K30  
    ! If at the node we have density and viscosity increments
    else if( i==k .and. j==l) then  
        dhdp=K00
        RL_Ta=temp(I,J)                                                                                                 !D_IN: /CurrentT/ -> temp(I,J)                       !D_OUT: RL_Ta -> /RLarsson/  (not saved, only used internaly)
        
        if(lub_param .EQ. 1)then                                                                                        !D_IN: /Ref/ -> lub_param
            RO(I,J) =1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J))                                                                  !D_IN: /Rho/ -> RA1, RA2; /Visc/ -> PH			     !D_OUT: RO(I,J) -> /Current/ (not saved, could maybe use the stored values instead? Or save these?)
            EDA     =EXP(alpha*Pref/Z*(-1+(1+P(I,J)*PH/Pref)**Z))                                                       !D_IN: /Visc/ -> alpha, Z, Pref
            
            ! Density derivative with respect to P
            dRodP   = RA1*PH/((1.0+RA2*PH*P(I,J))**2)
            incre   = (dRodP)*h_cube/(ENDA*EDA)                                                                         !D_IN: /Visc/ -> ENDA
            
            ! Second derivative of density with respect to P
            d2ROdP2      = -2*RA1*RA2*PH*PH / (1+RA2*PH*P(I,J))**3
            
        else if(lub_param == 3) then
            ! Roelands equation acc Holmes et al. 2003 Transient elastohydrodynamic point contact analysis using a new coupled differential de¯ection method Part 1: theory and validation
            EDA=EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1))                                                               !D_IN: /Visc/ -> EDA0; /Holmes/ -> kH, xH
            
            ! Density Formulation acc Homes et al. 2003
            RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
            
            ! Density derivative with respect to P
            dRodP   = RA1*PH/((1.0+RA2*PH*P(I,J))**2)
            incre   = (dRodP)*h_cube/(ENDA*EDA)
            
            ! Second derivative of density with respect to P
            d2ROdP2      = -2*RA1*RA2*PH*PH / (1+RA2*PH*P(I,J))**3
                
        else if(lub_param .EQ. 4)then
            ! Roelands equation acc R.Larsson 2000
            EDA0    = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)                                                            !D_IN: /RLarsson/ -> RL_G0, S0;               !D_OUT: EDA0 -> /Visc/  (not saved however)
            Z       = Dz+Cz*log10(1.0+RL_Ta/135.0)                                                                      !D_IN: /RLarsson/ -> Dz, Cz;                  !D_OUT: Z -> /Visc/	  (not saved however)
            EDA     = EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
            
            ! Derivative of EDA = Epsilon with respect to P
            dEDAdp  = EDA * (LOG(EDA0)+9.67)*(1.0+5.1E-9*P(I,J)*PH)**(Z-1.0)*Z*5.1E-9*PH
            
            ! Derivative of EDA = Epsilon with respect to temperature
            dEDA0dtemp  = EDA0*log(10.0)* RL_G0*S0*(RL_Ta/135.0+1.0)**(S0-1.0)*(1/135)
            dZdtemp     = Cz/((135.0+RL_Ta)*log(10.0))
            dEDAdtemp1  = dEDA0dtemp * EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
            dEDAdtemp2  = EDA*(1/EDA0)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)*dEDA0dtemp
            dEDAdtemp3  = EDA*(LOG(EDA0)+9.67)*((1.0+5.1E-9*P(I,J)*PH)**Z)*log(1+5.1E-9*P(I,J)*PH)* dZdtemp
            dEDAdtemp   = dEDAdtemp1 + dEDAdtemp2 + dEDAdtemp3

                
            ! Density Formulation acc R.Larsson 2000
            EpsT       = EpsT0*exp(-RL_c*P(I,J)*PH)                                                                     !D_IN: /RLarsson/ -> EpsT0
            RO(I,J)    = (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0))                               !D_IN: /RLarsson/ -> RL_T0
            
            ! Density derivatives with respect to P and temperature
            dRodP0      = RA1*PH/((1.0+RA2*PH*P(I,J))**2) 
            dRodP1      = dRodP0 * (1.0-EpsT*(RL_Ta-RL_T0))  
            dRodP1      = dRodP1 + (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(-EpsT*(-RL_c*PH)*(RL_Ta-RL_T0)) 
            dRodtemp    = (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J))) *(-EpsT)
            drodp       = (drodp1+drodtemp*dtemp_dp)
            
            incre = (dRodP)*h_cube/(ENDA*EDA)-(dEDAdP+dEDAdtemp*dtemp_dP)*Ro(i,j)*h_cube/(ENDA*EDA**2)

            ! Second derivative of density with respect to P
            dRodp1      = -2*RA1*RA2*PH*PH / (1+RA2*PH*P(I,J))**3  
            dRodp1      = dRodp1 + 2*dRodP0 * (-EpsT*(-RL_c*PH)*(RL_Ta-RL_T0)) 
            dRodp1      = dRodp1 + (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(-EpsT*(-RL_c*PH)**2*(RL_Ta-RL_T0)) 
            d2ROdP2      = (drodp1+0.0)
            
        else if(lub_param .EQ. 5) then !material description by Ehret 1997
            T40     = 40
            Tg      = Tg0+YA1*LOG(1.0+YA2*PH*P(I,J))                                                                    !D_IN: /Yasutomi/ -> YA1, Tg0, YA@ 
            temp_y  = ((T40*(1-Tcvot)+Tcvot*RL_Ta)-Tg)
            pg      = 1.0/YA2*(EXP((1.0/YA1)*temp_y)-1.0)
            YF      = 1.0-YB1*LOG(1.0+YB2*PH*P(I,J))                                                                    !D_IN: /Yasutomi/ -> YB1, YB2
            
            IF (P(I,J)*PH .GT. pg) THEN
                EDA     = Yedag*EXP(Yalfag*(P(I,J)*PH-pg))                                                              !D_IN: /Yasutomi/ -> Yedag, Yelfag
                dEDAdP  = Yalfag*PH*EDA
                
                if( P(i,j)*ph .LT. 1.3*pg)then                                                                          !Better to overestimate the derivative
                    dYfdp       = -YB1*YB2*PH/(YB2*PH*P(i,j)+1.0)
                    dEDAdP      = EDA*log(10.0) *( -YC1*temp_y/(YC2+temp_y*YF) + temp_y* YC1* temp_y*YF / ((YC2+temp_y*YF)**2) )* dYfdp                                !D_IN: /Yasutomi/ -> YC1, YC2
                    dEDAdtemp   = EDA*log(10.0) *( -YC1*YF/(YC2+temp_y*YF)     + YF* YC1* temp_y*YF / ((YC2+temp_y*YF)**2) )
                endif
                
                if ( ( EDA .GT. 1e15 )  ) then
                    EDA=1e15
                endif 
                
            ELSE
                EDA33   = Yedag*10**(-(YC1*temp_y*YF) / (YC2+temp_y*YF)) 
                EDA22   = Yedag*EXP(Yalfag*1)
                
                if (EDA22 .lt. EDA33) then
                    EDA         = EDA22
                    dEDAdP      = 0
                    dEDAdtemp   = 0
                else
                    EDA         = EDA33  
                    dYfdp       = -YB1*YB2*PH/(YB2*PH*P(i,j)+1.0)
                    dEDAdP      = EDA*log(10.0) *( -YC1*temp_y/(YC2+temp_y*YF) + temp_y* YC1* temp_y*YF / ((YC2+temp_y*YF)**2) )* dYfdp
                    dEDAdtemp   = EDA*log(10.0) *( -YC1*YF/(YC2+temp_y*YF)     + YF* YC1* temp_y*YF / ((YC2+temp_y*YF)**2) )
                endif
                if ( ( EDA .GT. 1e15 )  ) then
                    EDA=1e15
                endif 
                
            ENDIF
        
            dRodP1       = RA1*PH/((1.0+RA2*PH*P(I,J))**2)
            dRodtemp    = 0 
            drodp       = (drodp1+drodtemp*dtemp_dp)
            
            if( EDA .LT. 1e15) then
                incre = (dRodP)*h_cube/(ENDA*EDA)-(dEDAdP+dEDAdtemp*dtemp_dP)*Ro(i,j)*h_cube/(ENDA*EDA**2)
            else
                incre=(dRodP)*h_cube/(ENDA*EDA)
            endif

            ! Second derivative
            d2ROdP2      = -2*RA1*RA2*PH*PH / (1+RA2*PH*P(I,J))**3
            
            if ( isnan( incre )  ) then  
                WRITE(4,*)'check EPS_1112 incre, at I=', I
                call Stop_to_large_out(0) ! The zero is used since the time is not defined here
                incre=0
            endif 
            
        else if(lub_param .EQ. 7 .or. Lub_param .EQ. 8) then ! Material description by Liu and Wang

            eda01       = log(EDA0)+9.67
            eda02       = (1.0+5.1E-9*P(I,J)*PH)**Z
            eda03       = ((RL_Ta-138.0)/(RL_T0-138.0))**(-S0)
            EDA         = EDA0*exp(eda01 * (-1.0 + eda02 * eda03))
            
            dEDAdP      = EDA * eda01 * (1.0+5.1E-9*P(I,J)*PH)**(Z-1.0) * 5.1E-9*PH*Z *eda03
           
            dEDAdtemp   = EDA * eda01 * eda02 * (-S0)*((RL_Ta-138.0)/(RL_T0-138.0))**(-S0-1.0) * 1.0/(RL_T0-138.0)
            
            dRodP1       = RA1*PH/((1.0+RA2*PH*P(I,J))**2)
            dRodtemp    = - EpsT0
            drodp       = (drodp1+drodtemp*dtemp_dp)
            
            incre = (dRodP)*h_cube/(ENDA*EDA)-(dEDAdP+dEDAdtemp*dtemp_dP)*Ro(i,j)*h_cube/(ENDA*EDA**2)
            
            ! Second derivative
            d2ROdP2      = -2*RA1*RA2*PH*PH / (1+RA2*PH*P(I,J))**3
            
        else if(lub_param .EQ. 9) then ! Material description according to Bruyere

            eda01   = log(EDA0)+9.67
            eda02   = (1.0+P(I,J)*PH/(1.98*1E8))**Z
            eda03   = ((RL_Ta-138.0)/(RL_T0-138.0))
                
            delta_s = eda01*eda02*S0/(RL_T0-138.0)
            eda04   = delta_s*(RL_Ta-RL_T0)
                
            EDA     = EDA0*exp(eda01 * ((-1.0 + eda02 )* eda03)-eda04)
            dEDAdp  = EDA * eda01 * eda03* (1.0+P(I,J)*PH/(1.98*1E8))**(Z-1.0)*Z*PH/(1.98*1E8)*(1-S0)
            dEDAdtemp   = EDA*eda01*(-1.0+ eda02)* ( 1.0/(RL_T0-138.0) - delta_s )
                
            ! Density Formulation
            RO(I,J)     = (5.9*1E8 + RA1*PH*P(I,J))/(5.9*1E8 + RA2*PH*P(I,J)) - EpsT0*(RL_Ta-RL_T0)
            dRodP1       = 5.9*1E8*(RA1-RA2)*PH/((5.9*1E8+RA2*PH*P(I,J))**2)
            dRodtemp    = - EpsT0
            drodp       = (drodp1+drodtemp*dtemp_dp)
            
            incre = (dRodP)*h_cube/(ENDA*EDA)-(dEDAdP+dEDAdtemp*dtemp_dP)*Ro(i,j)*h_cube/(ENDA*EDA**2)
            
            ! Second derivative
            d2ROdP2      = -2.0*5.9*1E8*(RA1-RA2)*PH*RA2*PH/((5.9*1E8+RA2*PH*P(I,J))**3)
            
        else if(lub_param .EQ. 12 .or. Lub_param .EQ. 13) then ! Material description by Kim 2001
            eda01       = log(EDA0)+9.67
            eda02       = (1.0+5.1E-9*P(I,J)*PH)**Z
            eda03       = S0*(Rl_Ta-RL_T0)
            EDA         = EDA0*exp(eda01 * (-1.0 + eda02) - eda03)
            dEDAdP      = EDA * eda01 * (1.0+5.1E-9*P(I,J)*PH)**(Z-1.0) * 5.1E-9*PH*Z
            dEDAdtemp   = EDA * (-S0)
            
            dRodP1       = RA1*PH/((1.0+RA2*PH*P(I,J))**2)*(1-EpsT0*(RL_Ta-RL_T0))
            dRodtemp    = - EpsT0*(1+RA1*PH*P(i,j)/(1+RA2*PH*P(I,J)))
            drodp       = (drodp1+drodtemp*dtemp_dp)
            
            incre = (dRodP)*h_cube/(ENDA*EDA)-(dEDAdP+dEDAdtemp*dtemp_dP)*Ro(i,j)*h_cube/(ENDA*EDA**2)
            
            ! Second derivative
            d2ROdP2      = -2*RA1*RA2*PH*PH / (1+RA2*PH*P(I,J))**3*(1-EpsT0*(RL_Ta-RL_T0))
        else
            WRITE(4,*)'Bad lubrication number in subroutine Epsilon_derivative. lub_param ='
            WRITE(4,*) lub_param
            stop 'Bad Lubrication'
        endif
    
    ! Else we're too far awa for the derivatives to count
    else 
        dhdp=0.0
    endif
    
    if ( isnan( d2ROdP2 )  ) then
        WRITE(4,*)'check eps_2234 d2ROdP2 )  , at I=', I
    endif 
    
    IF( k==i .and. l==j) then
        dRodp_mat(i,j)=dRoDp                                         !D_OUT: dRodP_mat -> /CurrentRO/
    endif
    
    if (EDAx(i,j) .GT. 100)then                                      !D_IN: /Current/ -> EDAx
        incre=incre + 3*H(i,j)**2*dhdP*Ro(i,j)/(ENDA*EDAx(i,j))      
    endif
                
    return
    end
    
    
    
    
            
        
    