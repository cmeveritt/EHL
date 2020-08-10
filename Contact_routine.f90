! Subrutine for solving the metal contacts. If meatal contact is reached this subroutine detects it and increases the pressure in the nodes in contact. 
    ! The pressure is raised untill the defined minimum film thickness is reached. 
	SUBROUTINE Contact_routine(H00, t, SS, NYs, NN, Plim)    
        implicit none 
        include     'inc_Asp.h'
        include     'inc_COMAK2D.h'
        include     'inc_COMAKline.h'
        include     'inc_Contact_mat.h'
	    include     'inc_CurrentH.h'
        include     'inc_CurrentP.h'
	    include     'inc_Grid.h'
	    include     'inc_Itt.h'
        include     'inc_Numbers_of_bad.h'
        include     'inc_P_line_side.h'
	    include     'inc_Ref.h'
        include     'inc_Cores_used.h'
        
        ! Input
        real        H00
        integer     SS
        integer     NN, NYs
        ! Calculation param
        real        sum_Line_AKSS, AK_sum
        Integer     I,J, JJ, iter, n_cont,n_cont_line, l_cont, I2, J2, IF_contact
        real        FC, inkre, minimal, H_diff
        real        AK_cvot, AK_cvot2, AK_cvot3, AK_cvot4, AK_cvot5, AK_cvot6, Time_scale
        integer     i3, j3, i4, j4
        ! Other                 
        Integer     t
        Real        Plim, test, test2(1:601), test3(1:601)
        ! Output
        SAVE        /Numbers_of_bad/
        SAVE        /CurrentP/                    
        SAVE        /CurrentH/
        SAVE        /Contact_mat/
        SAVE        /P_line_side/

        l_cont   = 0 
        contact  = 0               ! Resetting the contact matrix               !D_OUT: contact -> /Contact_mat/
        
        ! The matrix AK correlates how much a node will deform based on a pressure at the same or another node. 
        ! The distance is the indexis. So AK(3,4) gives the deformation of the node that is 3 steps away in RD and 4 steps away in TD. 
        ! AK_cvotX is thus the cvot between the deformation in the current node and the node next to it
        AK_cvot  = AK(1,0)/AK(0,0)                                              !D_IN: /COMAK2D/ -> AK
        AK_cvot2 = AK(2,0)/AK(0,0) 
        AK_cvot3 = AK(3,0)/AK(0,0)
        AK_cvot4 = AK(4,0)/AK(0,0)
        AK_cvot5 = AK(5,0)/AK(0,0)
        AK_cvot6 = AK(6,0)/AK(0,0)
        
        ! The deformation is sometimes calculated for the whole line of nodes in the TD. 
        ! In order to restrict the derivatives in the TD for simulations of cylinder contacts.
        IF(SS .EQ. 1) sum_Line_AKSS = sum_Line_AK*PAIAK*DX*SS                   !D_IN: /Ref/ -> PAIAK; /GRID/ -> DX; /COMAKline/ -> sum_line_ak, sum_line_ak2, sum_line_ak4, sum_line_ak8;
        IF(SS .EQ. 2) sum_Line_AKSS = sum_Line_AK2*PAIAK*DX*SS
        IF(SS .EQ. 4) sum_Line_AKSS = sum_Line_AK4*PAIAK*DX*SS
        IF(SS .EQ. 8) sum_Line_AKSS = sum_Line_AK8*PAIAK*DX*SS

    ! If cylinder, check j=NN and take that as valid for the whole line of nodes in TD. If timedependent, the pressure may be y-dependent if not line asperity as asp_shape=6      
    If((Geom .EQ. 1 .OR. Geom .EQ. 5) .and. (t .LE. 0 .OR. ((asp_shape .EQ. 1 .or. asp_shape .EQ. 3 .or. asp_shape .EQ. 6 ).and. contact_alg .eq. 1) ))then         !D_IN: /Ref/ -> Geom, asp_shape, contact_alg;
        ! Reset the itteration counter
        iter=0
        ! Check at J=NN and asume that if contact here, it's valid for the entire line. So increase the entire line load.
        ! The pressure of J=9 is used as side pressure.
134     J=NN                     
        
        ! FC is a relaxation factor to restrict the pressure increments so that the contact is not released straight away
        FC=0.7
        IF ( iter .GE. 1)  FC = 0.85
        IF ( iter .GE. 3)  FC = 1.0
        IF ( iter .GE. 5)  FC = 1.5
        
        IF( t .LE. 0)then
            ! At the start the minimim film thickness is scaled up to increase the convergence rates. Since no metal contact should occure for the smooth contact of the time independent case. 
            Time_scale = 5.0
        else
            Time_scale = 1.0
            FC = 0.3
            IF ( iter .GT. 0)  FC = 0.4
            IF ( iter .GT. 1)  FC = 0.6
            IF ( iter .GT. 4)  FC = 0.7
            IF ( iter .GT. 7)  FC = 1
            IF ( iter .GT. 10)  FC = 1.2
        endif
        
        DO I=1+SS,NX-SS,SS                                                          !D_IN: /Grid/ -> Nx
            ! If Contat add pressure to the whole line. This should only be used in the begining and not for the time dep part. Therefore a margine is used. Time dep simulation should not be this close to the surface   
            IF (H(I,J) .LT. Hminimum*Time_scale) THEN                               !D_IN: /CurrentH/ -> H(I,J); /asp/ -> Hminimum
                    l_cont=l_cont+1 
                    inkre  =  FC*(-(H(I,J)-Hminimum*Time_scale)/sum_Line_AKSS)
                    
                    ! Add a limit on the pressure update to sabilize the solution. Same as for point evaluation below
                    If( t .gt. -11) then
                        if( iter .LT. 3) then
                            if (inkre .GT. 0.03) inkre=0.03                         
                        elseif( iter .LT. 10) then
                            if (inkre .GT. 0.05) inkre=0.05
                        else
                            if (inkre .GT. 0.1) inkre=0.1
                        endif
                    else
                        if (inkre .GT. 0.1) inkre=0.1
                    endif
                    
                    ! Check to ensure positive pressure updates. Othervise something is wrong
                    if (inkre .LT. 0.0) then
                        WRITE(4,*)'Negative pressure increment at line', I ,J,  '. inkre was ', inkre,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                        inkre=0.0  
                        Call Stop_to_large_out(t)                  
                    endif
                    
                    ! Evaluate the change in film thickness and stor the maximum value
                    H_diff= H(I,J)-Hminimum
                    Max_contact=max(max_contact, -(H_diff))                !D_OUT: Max_contact -> /Numbers_of_bad/
                       
                    DO JJ=NYs,1,- SS
                        ! AK only counts mesh steps. So independent of whitch mesh zise we're at, one step away is always AK(1,0)
                        
                        ! If at time dependent simulation, update each node specific    
                        IF( t .GT. 0) THEN
                            P(I,JJ) = P(I,JJ) + inkre    !D_OUT:  p(I,J) -> /CurrentP/
                        
                        ! If not in the time dependent regime, set the whole line to the same value. If to low, the next line will hopefully take care of it.     
                        ELSE
                            P(I,JJ) = P_line(I) + inkre     !D_IN: /P_line_side/ -> P_line
                        ENDIF
                        
                        ! Ensure that the presure is not thigher than 80% of the maximum pressure allowed. 
                        IF(P(I,JJ) .GT. Plim*0.8) P(I,JJ)=Plim*0.8
                    enddo
                    
                    ! Ad increment to the pressure outside the contact. Only calide for cylinder contacts
                    P_line(I)   = P_line(I) + inkre                                 !D_OUT: P_Line -> /P_line_side/
                    
                    ! Update the elastic deformations
                    ! VIFFT uses the fast fourier transforms to faster calculate the elastic defromations. 
                    if(.not. use_fft) then
                        CALL VI(SS, NYs, NN, t, 0, 1)                               !CALL: SS, NYs, NN = NYs+1/2, t, iter = 0, full = 1
                    else
                       CALL VIFFT(SS, NYs, NN, t)
                    endif

                    ! Update the the surface shape including the film thickness
                    CALL Asp_calc(H00, t,SS, NYs, NN)                               !CALL: H00, t, SS, NYs, NN = Nys+1/2
                ENDIF
        ENDDO
        
        ! Count the number of itterations and evaluate if more ittrations are needed
        iter=iter+1
        minimal=Hminimum
        DO I=1,NX,SS
            DO J=1,NN,SS
                IF( H(I,J) .LT. minimal) minimal=H(I,J)
            ENDDO
        ENDDO
        IF(minimal .LT. Hminimum .AND. iter .LT. 20) Goto 134                       !Branch 134 - can be replaced with do-while
        
        ! Stor the number of contacts. ! Since iter is 1 the first time. Then this will only include extra itterations 
        bad_h_itt_line = bad_h_itt_line + iter-1                                    !D_OUT: bad_h_itt_line -> /Numbers_of_bad/
        bad_h_nr_line  = bad_h_nr_line  + l_cont                                    !D_OUT: bad_h_nr_line -> /Numbers_of_bad/
    ENDIF

    IF(iter .GE. 20) then
        WRITE(4,*)'WARNING. Not hight enough H after 20 contactsitterations of the line sort!'  
    endif
    
    ! Reset the counter and check the lubricationfilm at each node
    iter=0
    n_cont=0
    IF( t .LE. 0 .OR. contact_alg .eq. 1) THEN
555     IF (t .lt. Ntime/2) then                           
            IF(SS .GT.1 ) Then
                FC=1
            else
                FC=0.9
            endif

            IF ( iter .GT. 0)  FC = 1 
            IF ( iter .GT. 5)  FC = 1.5
            IF ( iter .GT. 10)  FC = 2.0
            IF ( iter .GT. 15)  FC = 2.5
            IF ( iter .GT. 20)  FC = 3.1
                        
        else ! Take it more carfull in the outlet since then the simulations are more unstable
            FC=0.5
            IF ( iter .GT. 0)  FC = 0.75
            IF ( iter .GT. 1)  FC = 0.9
            IF ( iter .GT. 2)  FC = 1
            IF ( iter .GT. 5)  FC = 1.5
            IF ( iter .GT. 10)  FC = 2.0
            IF ( iter .GT. 15)  FC = 2.5
            IF ( iter .GT. 20)  FC = 3.1            
        endif
    
        ! The contact routine performes badly if a lot of contacts. Therefore the pressure increments are reduced based on the previous number of contacts. 
        IF(bad_h_nr_node/bad_h_itt_node .gt. 100 .and. iter .lt. 5) then
            FC=FC*0.5
        elseif(bad_h_nr_node/bad_h_itt_node .gt. 50 .and. iter .lt. 2) then
            FC=FC*0.75
        endif
        
    ! Initially assume that there is no contact   
    ! Check the contact in different directions depending on the time step. Probably not needed but just to ensure no directional dependency     
    DO I=1+SS,NX-SS,SS
        IF_Contact = 0         
        IF (MOD(iter,4)==0 .or. MOD(iter+1,4)==0 ) THEN
            IF (MOD(iter,2)==0) THEN
                I2=NX+1-I
            else
                I2=I
            endif
        else
            IF (MOD(iter,2)==0) THEN
                I2=I
            else
                I2=NX+1-I
            endif
        endif
        
        n_cont_line = 0          
        DO J=1,NN,SS 
                If ( t .LT. 0) then  ! Start looking for contact from the centrum line. Makes a different if it´s a ball
                    J2=NN+1-J 
                ELSEIF (MOD(iter,4)==0 .or. MOD(iter+1,4)==0 ) THEN
                     IF (MOD(iter,2)==0) THEN
                         J2=J
                    else
                        J2=NN+1-J
                    endif
                 else
                    IF (MOD(iter,2)==0) THEN
                        J2=NN+1-J
                    else
                        J2=J
                    endif
                 endif
                 JJ=NYs+1-J2
                 
                ! Check if Contat but use a bit of margin here to not need to update the deformations all the time. 
                IF (H(I2,J2) .LT. 0.8*Hminimum) THEN 
                    AK_sum=1.0
                    IF_Contact=1
                    
                    ! Define contact in the contact matrix. This is later used when calculating the viscosities and ensures that there is a high viscosity in the contact. 
                    ! Real contact = 2, near contact = 1
                    contact(I2,J2)=2                                          
                    
                    ! Update deformation in neibouring nodes
                    do j3 = -4, 4
                        j4=j2+j3*SS
                        if (j4 .lt. 1) j4=1
                        if( j4 .gt. Nys) j4=NYs
                        do i3 = -4, 4
                            i4=i2+i3*SS
                            if(i4 .lt. 1) i4=1
                            if(i4 .GT. Nx) i4=Nx
                            
                            ! Since now near the contact, contact =1
                            if( contact(i4,j4)==0) contact(i4,j4) =1 
                        enddo
                    enddo
                    
                    IF ( J2 == NN-ss    .or. J2 == 1 +  ss)  AK_sum=AK_sum+AK_cvot2 !!FC=FC/(1+AK_cvot2)                            ! Due to the symmetry line
                    IF ( J2 == NN-2*ss  .or. J2 == 1 +2*ss)  AK_sum=AK_sum+AK_cvot3 !!FC=FC/(1+AK_cvot3)                            ! Due to the symmetry line
                    IF ( J2 == NN-3*ss  .or. J2 == 1 +3*ss)  AK_sum=AK_sum+AK_cvot4 !!FC=FC/(1+AK_cvot4)                            ! Due to the symmetry line
                    IF ( J2 == NN-4*ss  .or. J2 == 1 +4*ss)  AK_sum=AK_sum+AK_cvot5 !!FC=FC/(1+AK_cvot5)                            ! Due to the symmetry line
                    IF ( J2 == NN-5*ss  .or. J2 == 1 +5*ss)  AK_sum=AK_sum+AK_cvot6 !!FC=FC/(1+AK_cvot6)                            ! Due to the symmetry line
                    
                    ! If nearby nodes also in contact, the pressures will also be increased on these nodes.
                    IF (H(I2,J2+SS) .LT. Hminimum) AK_sum=AK_sum+AK_cvot !!FC=FC/(1+AK_cvot)                   
                    IF (H(I2+SS,J2) .LT. Hminimum) AK_sum=AK_sum+AK_cvot !!FC=FC/(1+AK_cvot)
                    IF (H(I2-SS,J2) .LT. Hminimum) AK_sum=AK_sum+AK_cvot !!FC=FC/(1+AK_cvot)
                    
                    ! Since the normal evaluation cant be performed for J=1 the values of j=9 are used for both sides. The pressure of j=9 is also used for the side pressures oc cylinder contacts. 
                    IF(J2 .GT. 1) THEN 
                        IF (H(I2,J2-SS) .LT. Hminimum) AK_sum=AK_sum+AK_cvot !!FC=FC/(1+AK_cvot)             
                    Else
                        IF (H(I2,J2+SS) .LT. Hminimum) AK_sum=AK_sum+AK_cvot !!FC=FC/(1+AK_cvot)             
                    ENDIF
                    
                    ! Evaluate the current and maximum penetration depth
                    H_diff= H(I2,J2)-Hminimum
                    Max_contact=max(max_contact, -(H_diff))
                    
                    ! Calculate the pressure increment including the relaxation factor FC
                    inkre=-FC/AK_sum*(H_diff)/(AK(0,0)*PAIAK*DX*SS) 
                    
                    ! Ad a limit on the pressure updation to sabilize the solution. 
                    If( t .gt. -11) then
                        if( iter .LT. 3) then
                            if (inkre .GT. 0.03) inkre=0.03                         
                        elseif( iter .LT. 10) then
                            if (inkre .GT. 0.05) inkre=0.05
                        else
                            if (inkre .GT. 0.1) inkre=0.1
                        endif
                    else
                        if (inkre .GT. 0.1) inkre=0.1
                    endif
                    
                    ! Check to ensure positive pressure updation
                    if (inkre .LT. 0.0) then
                        WRITE(4,*)'Negativ pressure inkrement at', I ,J,  '. inkre was ', inkre,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                        inkre=0.0                                           
                        Call Stop_to_large_out(t)
                    endif
                    
                    ! Increasing the pressure at al nodes in contact
                    P(I2,J2)    = P(I2,J2)+inkre 
                    
                    ! Count the number of contacts
                    n_cont      = n_cont+1                                         
                    n_cont_line = n_cont_line +1
                    
                ENDIF
                
                ! Evaluate the pressures compared to the global pressure limit. Here 30% higher pressures are allowed since the metal contact potentially could generate quite severe conditions. 
                ! Also ensure only positive pressures
                IF(P(I2,J2) .GT. Plim*1.3)then
                    P(I2,J2)=Plim*1.3 
                    Bad_p_node=Bad_p_node+1                               !D_OUT: Bad_p_node -> /Numbers_of_bad/
                elseIF(P(I2,J2).LT.0.0)THEN
                    P(I2,J2)=0.0
                    Bad_p_node=Bad_p_node+1
                ENDIF
                
                !Mirror the pressure
                P(I2,JJ)=P(I2,J2)
        ENDDO
                
        ! If contact somewhere along the line of constant x value, Update the deformation.
        IF( IF_Contact .EQ. 1) THEN             
            
            ! If analytical asperities in contact, the contact presure should be higher at the central line than the nearby line. 
            if( asp_shape .lt. 200 .or. asp_shape .gt. 209) CALL P_centre_fix    
            
            !Update the elastic deformation
            if( n_cont_line .LT. 5) Then
                if(.not. use_fft) then
                    ! Do not update the side pressures and do only update the nodes with contact(i,j) > 0
                    CALL VI(SS, NYs, NN, t, 1, 0)                                               !CALL: SS, Nys, NN = Nys+1/2, t, iter = 1, full = 0
                else
                    ! Since VIFFT uses fast fourier transforms the deformations of the whole surface has to be calculated
                    CALL VIFFT(SS, NYs, NN, t)
                endif
                
            else
                ! Update the pressure outside of the contact. Could be optimized
                If(Geom .EQ. 1 .OR. Geom .EQ. 5) CALL Pressure_outside(t,SS, NYs, 1, 1)         !CALL: t, ss, Nys, k = 1, k_use = 1
                
                ! Update the deformation of all nodes 
                if (.not. use_fft) then
                    CALL VI(SS, NYs, NN, t, 1, 1)                                               !CALL: SS, NYs, NN = Nys+1/2, t, iter = 1, full = 1
                else
                    CALL VIFFT(SS, NYs, NN, t)
                endif
            endif
            
            ! Calculate the new fiilmthickness based on the shape of the surface irregularities
            CALL Asp_calc(H00, t,SS, NYs, NN)                                                   !CALL: H00, t, SS, NYs, NN = NYs + 1 /2 
            IF_Contact=0
            n_cont_line=0
        ENDIF
    ENDDO

    ! Count the number of iterations and find the minimum film thickness
    iter=iter+1
    minimal=Hminimum
    DO I=1,NX,SS
        DO J=1,NN,SS
            IF( H(I,J) .LT. minimal) minimal=H(I,J)
        ENDDO
    ENDDO
            
    IF(iter .GE. 40) then
        WRITE(4,*)'WARNING. Not hight enough H after 40 contact itterations'  
    endif
    
    ! Stor the numer of contacts for the last itteration with contacts. These tells the program how delecate the problem is in are used for estimating the relaxation parameters
    IF( Max_contact .gt. 0) THEN  
        bad_h_nr_node  =  n_cont  
        bad_h_itt_node =  iter 
    ENDIF
    
    ! If to low H then do it again.
    IF( minimal .LT. Hminimum*0.7 .and. iter .lt. 100) Goto 555          
    
    if (iter .GE. 100 ) then
        WRITE(4,*)'Not converged contact at', I ,J,  '. inkre was ', inkre,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        Call Stop_to_large_out(t)
    endif
        
    ! Remove time dependent efects if in the time indep regime 
    IF( (n_cont .GE. 1 .or. l_cont .GE. 1) .and. t .LE. 0) CALL pastupd(0, NX, NN, SS)   !CALL: a = 0, Nx, NN, SS
    ENDIF
    
    bad_h_end_node=0                !D_OUT: bad_h_end_node -> /Numbers_of_bad/
    
    ! last resort fix if used 100 contact itterations. Even if the contact is not solven now it might be solved for futher itterations. 
    DO J=1,NN,SS   
        DO I=1+SS,NX-SS,SS
            IF (H(I,J) .LT. 0.6*Hminimum .and. t .LE. 0) THEN         
                !WRITE(4,*) 'WARNING! Now at end of contact routine and still bad H !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
                !WRITE(4,*)  'at I,J= ', I, J, 'H was ', H(I,J)
                !Call Stop_to_large_out(t)
                H(I,J)=0.6*Hminimum
                contact(I,J)=2
                bad_h_end_node=bad_h_end_node+1
            ENDIF
        enddo
    enddo
    
    RETURN
    END
