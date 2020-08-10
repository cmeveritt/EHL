!Subroutine for initiating important parameters of the problem
    ! The initial guesses are based upon Hertz contact theory for dry contacts. 
    SUBROUTINE INITI(H00)
        implicit none
        include     'inc_AKparam.h'
        include     'inc_COMH.h'
        include     'inc_COMAK2D.h'
        
        include     'inc_Current.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_CurrentT_met.h'
        
        include     'inc_DZ_com.h'
        include     'inc_Error.h'
        
        include     'inc_Grid.h'
        include     'inc_G0DT.h'
        include     'inc_Geom5.h'
        
        include     'inc_Limit_pressure.h' 
        include     'inc_Outp.h'
        
        include     'inc_PastT.h'
        include     'inc_PastT_met.h'
        include     'inc_P_line_side.h'
        include     'inc_P_new_com.h'
        
        include     'inc_Ref.h'
        include     'inc_Residual.h'
        include     'inc_Setup.h'
        include     'inc_Temp_param.h'
        
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'
        include     'inc_Y_liu.h'
        ! Input
        real        H00
        real        R_rat, C_rat
        ! Calculation parameters
        Integer     I, J, NN, Nys, JJ, k, im, jm
        real        C,D
        real        time, umax, dxssb, SRR, temp_loc
        ! Other
        real        Y0, temp_convm
        ! Output
        SAVE        /Current/                                     
        SAVE        /CurrentP/ 
        SAVE        /Setup/                                         
        SAVE        /AKparam/
        SAVE        /Grid/
        SAVE        /Error/
        SAVE        /G0DT/
        SAVE        /Ref/
        SAVE        /COMH/ 
        SAVE        /CurrentT_met/
        SAVE        /PastT_met/
        SAVE        /CurrentT/
        SAVE        /PastT/
        save        /P_line_side/
        save        /Limit_pressure/ 
        save        /P_new_com/
        
        Plim=15                                                     ! Maximum normalized pressure. Used in HREE, P_update and Initial guess                 !D_OUT: Plim -> /Limit_pressure/
        Ry = 1
        By = 1
        R_rat=RX/Ry                                                 ! Elliptical ratio                                                                      !D_IN: /Outp/ -> Rx, Ry
        C_rat=b/By                                                  ! Elliptical ratio                                                                      !D_IN: /Outp/ -> B, By
        !R_rat=C_rat**2                                             ! something is wrong with this formulation
        NYs=Ny                                                                                                                                              !D_IN: /Grid/ -> Ny
        NN=(NY+1)/2
        Y0=-0.5*DX*NY+0.5*DX                                        
        
        umax=max(ua,ub)                                                                                                                                     !D_IN: /Outp/ -> Ua, Ub
        DXSSB=dx*1*b
        SRR=2*(Ua-Ub)/(Ua+ub) 
        
        DO I=1,NX
            X(I)=X0+(I-1)*DX                                                                                                                                !D_OUT: X -> /Setup/
        ENDDO
        
        DO J=1,NY
            Y(J)=Y0+(J-1)*DX                                                                                                                                !D_OUT: Y -> /Setup/
        ENDDO
        
        IF(Geom .EQ. 2 .or. Geom .EQ. 6)Then    ! Ball                                                                                                         !D_IN: /Ref/ -> Geom
            DO J=1,NN
                D=1.-Y(J)*Y(J)
                DO  I=1,NX
                    C=D-X(I)*X(I)               ! Include or not depending on if Cylinder or Ball
                    IF(C.LE.0.0)P(I,J)=0.0                                                                                                                  !D_OUT: P(I,J) -> /CurrentP/
                    IF(C.GT.0.0)P(I,J)=SQRT(C)
                    RAD(I,J)=0.5*(X(I)*X(I)+Y(J)*Y(J))                                                                                                      !D_OUT: RAD(I,J) -> /COMH/
                ENDDO
            ENDDO
            
        ELSE IF(Geom .EQ. 3)Then                ! Elliptical Ball
            DO J=1,NN
                D=1.-Y(J)*Y(J)*C_rat*C_rat
                DO  I=1,NX
                    C=D-X(I)*X(I)               ! Include or not depending on if Cylinder or Ball
                    IF(C.LE.0.0)P(I,J)=0.0
                    IF(C.GT.0.0)THEN
                        P(I,J)=SQRT(C)
                       
                    ENDIF
                    RAD(I,J)=0.5*(X(I)*X(I)+Y(J)*Y(J)*R_rat)!*R_rat)
                ENDDO
            ENDDO
                                                ! The integral of the dimensionless load. if b=by this should be =2/3 pai
        ELSE                                    ! Cylinder
            DO  I=1,NX
                D=1.-X(I)*X(I)
                DO J=1,NN
                    C=D
                    IF(C.LE.0.0)P(I,J)=0.0
                    IF(C.GT.0.0)P(I,J)=SQRT(C)
                    RAD(I,J)=0.5*(X(I)*X(I))
                ENDDO
                P_line(I)=P(I,1)
            ENDDO
            IF( Geom .EQ. 5) then
                P=P*(PH_new/PH)**2              ! A cours rescaling of the pressure
            endif
            
        ENDIF

        DO J=NN+1,NY
            JJ=NY-J+1
            DO I=1,NX
                P(I,J)=P(I,JJ)
                RAD(I,J)=RAD(I,JJ)
            ENDDO
        ENDDO
        
        if ( Geom .EQ. 6) then                  ! Rescale the initial guess since should start really low. 
            P=0.001*P
        endif
    
     
        Pold  = P                               ! Pressure at past itteration                                                                                 !D_OUT: Pold -> /Error/
        
        ! Calculate the applied load
        if(Geom .EQ. 3) G0=sum(P)*(DX*DX)                                                                                                                     !D_OUT: G0 -> /G0DT/

        ! Oil temperature
        if (lub_temp .NE. 0) then                                                                                                                             !D_IN: /temp_param/ -> lub_temp 
            ! Initializing the temperature fields to be the outside temperature
            Temp_ma=Ta                                                                                                                                        !D_OUT: Temp_ma -> /CurrentT_met/             !D_IN: /Yasutomi/ -> Ta
            Temp_mb=Ta                                                                                                                                        !D_OUT: Temp_ma -> /CurrentT_met/ 
            Temp_map=Ta                                                                                                                                       !D_OUT: Temp_map -> /PastT_met/
            Temp_mbp=Ta                                                                                                                                       !D_OUT: Temp_mbp -> /PastT_met/
            
            ! Temperature increase is quite proportionall to the pressure. So assume a initial temperature distribution based on pressure
            if( SRR == 0) then
                temp_r=5                                                                                                                                      !D_OUT: temp_r -> /temp_param/
            else if(lub_param .LT. 9)then                                                                                                                     !D_IN: /ref/ -> lub_param
                temp_r=temp_r*sqrt(abs(SRR))
            else 
                temp_r=30
            endif
            
            !Linear temperature increase
            IF( asp_shape .eq. 150)then                                                                                                                       !D_IN: /ref/ -> asp_shape
                asp_shape=12                                                                                                                                  !D_OUT: asp_shape -> /ref/
                DO J=1,NY
                    DO I=1,NX
                        
                        if(Geom == 2)then ! Ball
                            IF(X(I) .GT. -1.0 .and. X(I) .LE. 1.2 .and. abs(Y(J)) .LT. 1) then                                                                
                                Temp(I,J)=max(Temp(i,j), Ta+(X(I)+1.0)*120/2*(1-Y(J)**2))                                                                     !D_OUT: Temp(I,J) -> /CurrentT/
                            else
                                temp(i,j)=Ta                                                                                                                  !D_OUT: Tempp(I,J) -> /PastT/
                                tempp(i,j)=Ta
                            endif
                        else
                            IF(X(I) .GT. -1.0 .and. X(I) .LE. 1.2 ) then 
                                Temp(I,J)=max(Temp(i,j), Ta+(X(I)+1.0)*120/2)
                            else
                                temp(i,j)=Ta
                                tempp(i,j)=Ta
                            endif
                        endif
                    
                    ENDDO
                ENDDO
                
            !Quadratic temperature spike with shift backwards
            else IF( asp_shape .eq. 151)then 
                asp_shape=12
                DO J=1,NY
                    DO I=1,10
                        Temp(I,J)=Ta
                        Tempp(I,J)=Ta
                    ENDDO
                    
                    DO I=10,NX
                        Temp(I,J)=Ta+temp_r*P(i-5,j)
                        Tempp(I,J)=Ta+temp_r*P(i-5,j)
                    ENDDO
                ENDDO
                
            !Quadratic temperature spike with shift forwards
            else IF( asp_shape .eq. 152)then 
                asp_shape=12
                DO J=1,NY
                    DO I=NX-10,NX
                        Temp(I,J)=Ta
                        Tempp(I,J)=Ta
                    ENDDO
                    
                    DO I=1,NX-10
                        Temp(I,J)=Ta+temp_r*P(i+5,j)
                        Tempp(I,J)=Ta+temp_r*P(i+5,j)
                    ENDDO
                ENDDO
                
            ! Combination of linear and quadreatic
            else IF( asp_shape .eq. 153)then 
                asp_shape=12
                DO J=1,NY
                    DO I=1,NX
                        Temp(I,J)=Ta+temp_r*P(i,j)
                        
                        if(Geom == 2 .or. Geom == 6)then ! Ball
                            IF(X(I) .GT. -1.0 .and. X(I) .LE. 1.2 .and. abs(Y(J)) .LT. 1) then 
                                Temp(I,J)=max(Temp(i,j), Ta+(X(I)+1.0)*120/2*(1-Y(J)**2))
                            endif
                        else
                            IF(X(I) .GT. -1.0 .and. X(I) .LE. 1.2 ) then 
                                Temp(I,J)=max(Temp(i,j), Ta+(X(I)+1.0)*120/2)
                            endif
                        endif

                        Tempp(I,J)=temp(i,j)
                    ENDDO
                ENDDO
                
            ! Temperaure proportional Pressure
            else
                DO J=1,NY
                    DO I=1,(NX+1)/2
                        Temp(I,J)=Ta+temp_r*P(i,j)
                        Tempp(I,J)=Ta+temp_r*P(i,j)
                        
                    ENDDO
                    
                    ! Linear decrease of temperature from the central part of the contact
                    DO I=(NX+1)/2+1,NX
                        Temp(I,J)=Ta+temp_r*max(P(i,j),0.5+1.01*(1.0-(1.01*I)/NX))    
                        Tempp(I,J)=Temp(I,J)
                        
                    ENDDO
                ENDDO
            endif
            
            ! Update the temperature fields down in the metals based upone the temperature in the lubricants
            ! Since metal nodes closer to the oil they have to be warmer than in the other case      
            IF(DZ_method==1) THEN                      
                do k=1,30
                    do jm=1,(NY+1)/2
                        j=(jm-1)*2+1
                        do im=1,(NX+1)/2
                            i=(im-1)*2+1
                            temp_loc = Ta+(Temp(i,j)-Ta)*(1-k/30)
                            
                            Temp_ma(Im,Jm,k)=temp_loc
                            Temp_map(Im,Jm,k)=temp_loc
                        
                            Temp_mb(Im,Jm,k)=temp_loc
                            Temp_mbp(Im,Jm,k)=temp_loc
                        enddo
                    enddo
                enddo
            ELSE
                do k=1,5
                    do jm=1,(NY+1)/2
                        j=(jm-1)*2+1
                        do im=1,(NX+1)/2
                            i=(im-1)*2+1
                            temp_loc = Ta+(Temp(i,j)-Ta)/k
                            
                            Temp_ma(Im,Jm,k)=temp_loc
                            Temp_map(Im,Jm,k)=temp_loc
                        
                            Temp_mb(Im,Jm,k)=temp_loc
                            Temp_mbp(Im,Jm,k)=temp_loc
                        enddo
                    enddo
                enddo
            ENDIF
            
            time=10.0
        else
            ! Uniform tempreature in oil
            Temp=Ta
            Tempp=Ta
        endif 
        
        ! Uppdate the other parameters of the problen, such as the  viscosity and density
        CAll HREE(H00,-11, 1, NY, 40.0, 1.0, 0, 1, 1)                !Call D_OUT: H00, t=-11, SS = 1, Nys = Ny (even though its in /grid/), T40 = 40.0, k = 1.0, k_use = 1, M_conv = 1         
     
        ! Update the specific heat capacity
        CALL cp_calc(NX, NN,1)                                       !CALL: Nx, NN, SS = 1
        
        ! Metal temperature
        ! The time is scaled up to faster reach the time independent solution. The time independent solution could be obtained in a better manner. 
        if(lub_temp .NE. 0) then
            do I=1,10
                time=time*5
                CALL pastupd(2, NX, NN, 1)                           !CALL: a = 2, Nx, NN, SS = 1
                if( time .gt. 5*DXSSB/umax) time=5*DXSSB/umax        ! The maximum timestep is set so that the metal will not move more than 5 node step 
                CALL temp_calc_metal(time, 1, -11, NYs, temp_convm)  !CALL: dt_lim = time, SS = 1, t = -11, NYs, temp_convm 
                Call Met_temp_upd(time, NY, 1, TA, 0.0)              !CALL: dt_lim = time, NY, SS = 1, t0 = TA, time = 0.0 
                Call Met_temp_upd(time, NY, 1, TA, time)             !CALL: dt_lim = time, NY, SS = 1, t0 = Ta, time = time
            enddo
            I=1
            do while( I .LT. 20 .and. temp_convm .GT. 100)    ! !!! temp_convm is returned by temp_calc_metal
                time=time*5
                CALL pastupd(2, NX, NN, 1)
                if( time .gt. DXSSB/2*1/umax) time=DXSSB/2*1/umax  ! The maximum timestep is set so that the metal will not move more than half a node step 
                CALL temp_calc_metal(time, 1, -11, NYs, temp_convm) 
                Call Met_temp_upd(time, NY, 1, TA, 0.0)
                Call Met_temp_upd(time, NY, 1, TA, time)
                I=i+1
            enddo
        endif
        
        ! Parameters for the res subrutine in the ITER subrutine
        AK00=AK(0,0)       !D_IN: /COMAK2D/ -> AK(I,J) !D_OUT: AK00 -> /AKparam/
        AK10=AK(1,0)                                   !D_OUT: AK10 -> /AKparam/
        AK20=AK(2,0)                                   !D_OUT: AK20 -> /AKparam/
        AK30=AK(3,0)                                   !D_OUT: AK30 -> /AKparam/
        BK00=AK00-AK10                                 !D_OUT: BK00 -> /AKparam/
        BK10=AK10-0.25*(AK00+2.*AK(1,1)+AK(2,0))       !D_OUT: BK10 -> /AKparam/
        BK20=AK20-0.25*(AK10+2.*AK(2,1)+AK(3,0))       !D_OUT: BK20 -> /AKparam/
        BK30=AK30-0.25*(AK20+2.*AK(3,1)+AK(4,0))       !D_OUT: BK30 -> /AKparam/
        
        ! Ensure they are zero from the start
        A5dxx   = A5dxx*0.0                            !D_OUT: A5dxx -> /Residual/
        A5dx    = A5dxx*0.0                            !D_OUT: A5dx -> /Residual/
        A5dt    = A5dxx*0.0                            !D_OUT: A5dt -> /Residual/
        xi      = A5dxx*0.0                            !D_OUT: Xi -> /Current/
       
        CALL pastupd(0, NX, NN, 1)                     !CALL: a = 0, Nx, NN, SS = 1
        P_upd=P                                        !D_OUT: P_upd -> /P_new_com/
        
    if( asp_shape .ge. 200 .and. asp_shape .le. 209) then
        ! Then rough surface is used and that surface is already printed by the subroutine Asp_read
    else
        CALL Asp_print(DX, Y0, NN) ! Print the calculated asperity surface
    end if 

    RETURN
    END