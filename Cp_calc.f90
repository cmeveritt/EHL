! A subroutine to calculate the material parameters Cp and Kappa in the oil
    ! So far only implemented for lub_param=4. This condition is not enforced.
    ! The equations are based upone R. Larsson et al. 2000 Lubricant properties for input to hydrodynamic and elastohydrodynamic lubricatio analyses. 
    subroutine cp_calc(NX, NN,SS)
    implicit none
    include     'inc_Current.h'
	include     'inc_CurrentP.h'
    include     'inc_CurrentT.h'
    include     'inc_Therm_cond.h'
    include     'inc_Therm_matr.h'
    include     'inc_Visc.h'
    include     'inc_Cores_used.h'
   
    ! Input
    integer     NX, NN, SS
    ! Calculations
    Integer     i, i1, j, j1
    real        Pi, Pi1, Pj1,  Roi1, Roj1
    real*8      tempi, tempi1, tempj1
    real        beta
    ! Output    
    save        /Therm_matr/
    
    ! Loop over the the whole half surface
    !$OMP PARALLEL DO SHARED(NN,SS,NX,P,Ph,temp,RO,k_O,Cp_O,DCp_dtemp,ka0,ka1,ka2,be0,be1,be2,CP0,CP1,CP2)  &
    !$OMP&            PRIVATE(J,I,J1,I1,Pi,Pi1,Pj1,tempi,tempi1,tempj1,Roi1,Roj1,beta) &
    !$OMP&  IF(use_multiple_cores)
    Do j=1,NN,SS
        J1=J+1*SS
        IF( J1 .GT. NN) J1=NN-SS
        Do i=1,NX,ss
            I1=I+1*SS
            IF( I1 .GT. NX) I1=NX
            
            ! Extract the pressure and the average presure in the RD and TD directions
            Pi      = P(i,j)                                            !D_IN: /CurrentP/ -> P(I,J) 
            Pi1     = 0.5*(P(i,j)+P(i1,j)) 
            Pj1     = 0.5*(P(i,j)+P(i,j1))
            
            ! Highest load the lubrication was tested for. See R. Larsson et al. 2000 Lubricant properties for input to hydrodynamic and elastohydrodynamic lubricatio analyses. 
            ! At higer load beta decreases for lub_param=4 PAOB
            IF(Pi*Ph .GT. 1.1e9)  Pi= 1.1e9/Ph                          !D_IN: /Visc/ -> Ph  
            IF(Pi1*Ph .GT. 1.1e9)  Pi1= 1.1e9/Ph
            IF(Pj1*Ph .GT. 1.1e9)  Pj1= 1.1e9/Ph
            
            ! Extract the temperature and the average temperature in the RD and TD directions
            tempi   = temp(i,j)                                         !D_IN: /CurrentT/ -> Temp
            tempi1  = 0.5*(temp(i,j)+temp(i1,j))
            tempj1  = 0.5*(temp(i,j)+temp(i,j1))
            
            ! Highest temperature the lubrication was tested for, see R. Larsson et al. 2000 Lubricant properties for input to hydrodynamic and elastohydrodynamic lubricatio analyses.
            IF(tempi .GT. 107)   tempi= 107     
            IF(tempi1 .GT. 107)  tempi1= 107
            IF(tempj1 .GT. 107)  tempj1= 107
            
            ! The density
            Roi1  = 0.5*(Ro(i,j)+Ro(i1,j))                              !D_IN: /Current/ -> Ro
            Roj1  = 0.5*(Ro(i,j)+Ro(i,j1))
        
            ! Calculate the thermal conductivity
            ! The material parameters are already normalized with Ph in this IF statement
            k_O(i*2-1,j*2-1)      = ka0 * (1+ ka1*Pi/(1+ka2*Pi))        !D_IN: /Therm_cond/ -> ka0, ka1, ka2; !D_OUT: K_o -> /Therm_matr/
            k_O(i*2,j*2-1)        = ka0 * (1+ ka1*Pi1/(1+ka2*Pi1))
            k_O(i*2-1,j*2)        = ka0 * (1+ ka1*Pj1/(1+ka2*Pj1))

            ! Calculate the specific heat capacity Cp
            beta                  = be0 * (1+ be1*Pi +be2*Pi*Pi)        !D_IN: /Therm_cond/ -> be0, be1, be2
            Cp_O(i*2-1,j*2-1)     = 1/Ro(i,j)*CP0* (1+ beta* (temp(i,j)-22)*  (1+ Cp1*Pi/(1 + Cp2*Pi)))    !D_IN: /Therm_cond/ -> Cp0, Cp1, Cp2; !D_OUT: Cp_O -> /Therm_matr/
            DCp_dtemp(i,j)        = 1/Ro(i,j)*CP0* (beta*  (1+ Cp1*Pi/(1 + Cp2*Pi)))                       !D_OUT: DCp_dtemp -> /Therm_matr/
            beta                  = be0 * (1+ be1*Pi1 +be2*Pi1*Pi1)
            Cp_O(i*2,j*2-1)       = 1/Roi1*CP0* (1+ beta* (tempi1-22)      *  (1+ Cp1*Pi1/(1 + Cp2*Pi1)))
            beta                  = be0 * (1+ be1*Pj1 +be2*Pj1*Pj1)
            Cp_O(i*2-1,j*2)       = 1/Roj1*CP0* (1+ beta* (tempj1-22)      *  (1+ Cp1*Pj1/(1 + Cp2*Pj1)))  
            
            enddo
        enddo
    
    !$OMP END PARALLEL DO
    
    return
    end
    
