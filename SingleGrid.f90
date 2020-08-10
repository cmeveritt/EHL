! Subrutine for controling the pressure itterations at the current grid level of the pressure P based upon the Reynolds equation. 
    
	SUBROUTINE SingleGrid(t,H00,DW,kvot,DT2,SS, NYs, C1, T40, Tcvot, lim, M_conv)
    	implicit none
        include     'inc_Asp.h'
        include     'inc_CTRA4.h'
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
        include     'inc_Contact_mat.h'
        include     'inc_COMAK2D.h'
        include     'inc_G0DT.h'
        include     'inc_Geom5.h'
        include     'inc_Grid.h'
        include     'inc_Iterp.h'
        include     'inc_Itt.h'
        include     'inc_NonNew.h'
	    include     'inc_Method.h'
        include     'inc_Past.h'
        include     'inc_Rho.h'
        include     'inc_Ref.h'
        include     'inc_Residual.h'
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'
        include     'inc_P_new_com.h'
        real         Pcenter(NX,2), Clim(NX)                    !D_IN: /Grid/ -> Nx
        real         T40, Tcvot, NNSS, H00
        Integer      NN, t, SS, NYs, M_conv, C_meth
        integer     MM,K, J,JJ,J1,J0, I, I0,I1, I2, IA, lim, i4, K_use
        Real        PAI, eps1, eps2,  width, DW
        real        load, Sload, limit, limit2
        real        kvot, DT2, D2, D2p
        real        dHdx, dHdxp, Di_max
        real        KG1, C1, C2, DD, lim_check, C11             ! Apserity parameters
        real        AK_sum, AK_cvot, FC, incre
        real        test, test2(1:601), sum_d                                  
        SAVE      /Residual/                                    ! Parameters for the residual
        SAVE      /ITERp/
        SAVE      /CTRA4/
        SAVE      /Contact_mat/
        SAVE      /P_new_com/
        DATA      KG1/0/
        PAI     =3.14159265

        C11=C1
        C1=1.0                                                  ! Since adjusted in P_update
        C2=C1                                                   ! The lower factor, the more stable but the longer time
        lim=0
        C_meth = 1

        ! Initiating parameters
        KG1=1
    2   NN       = (NYs+1)/2              
        MM       = NX-1
        NNSS     = NN-SS
        AK_cvot  = AK(1,0)/AK(0,0)          !D_IN: /COMAK2D/ -> AK(I,J)
        
        
        If( t .LT. 0) then
            K_use=KK_stat                   !D_IN: /itt/ -> kk_stat
        else
            K_use=KK_time                   !D_IN: /itt/ -> kk_time
        endif
        K=1
        
        ! Itterate KK times untill checking if converged or not. 
        DO WHILE( K .LE. K_use)
            Di_max   = 0
            IF(C_meth .ge. 1) Then
                P_upd=P                     ! All updates of the pressure are stroed in P_upd. Then the subroutine P_update evaluates how much the pressure should be updated and decides on a suting relaxation factor.  !D_IN: /CurrentP/ -> P; !D_OUT: P_upd -> /p_new_com/
            endif
            
            DO I=1,NX,SS                    !Storing the old values of the pressure in the centerline and the line just before (NNSS). 
                Pcenter(I,1)=P(I,NNSS)
                Pcenter(I,2)=P(I,NN)
            ENDDO

            ! Start from NYs=2 ang to to the senter node in Y direction and uppdate the pressure
            DO 70 J=1+ss,NN+ss,SS
                ! Define nodenumbers for past and next nodes
                J0=J-1*SS
                IF( J0 .LT. 1) J0 = 1
                J1=J+1*SS
                JJ=NYs+1-J
                
                ! Check if its posible to save time by excluding nodes near the exit 
                IA=1*SS
    8           MM=NX-IA                           ! !!! Branch target 8
                
                IF(P(MM,J0).GT.1.E-6)GOTO 20       ! !!! Branch point to 20
                IF(P(MM,J).GT.1.E-6)GOTO 20
                IF(P(MM,J1).GT.1.E-6)GOTO 20
                IA=IA+1*SS
                IF(IA.LT.NX)GOTO 8                ! !!! Branch point to 8
                GOTO 70                           ! !!! Branch point to 70
    20          IF(MM.LT.NX-1*SS)MM=MM+1*SS 
        
                D2=0.5*(EPSx(1,J)+EPSx(1+SS,J))   !D_IN: /Current/ -> EPSx(I,J)
                D2p=0.5*(EPSxpast(1,J)+EPSxpast(1+SS,J)) !D_IN: /Past/ -> EPSxpast(I,J)
                
                ! Uppdate the presure for NX=2 untill there is zero pressure near the exit. at NX=1 BCs restric this pressure to be 0
                eps1=Hminimum                     !D_IN: /asp/ -> Hminimum
                eps2=Hminimum
                
                !Limiting the pressure updation
                if( t .GE. 0) Then              ! These limits are for when we get asperity contact. So should be indep of C1
                    limit  = 1.0/K_use          ! Limit dependent on C. D(I) is multiplied with C later. Scaled with KK so that on KK iteration the presure could change with 1
                    limit2 = 3.0*limit           ! release the pressure faster than what it can build up to decreas effect of H00 buildup 
                else                            ! These limits are for initial convergece. So can be dep of C1. Since if C1 islarge we're at a stable point. 
                    limit  = 1                  ! Could be dependent on number of nodes in contact at previous itteration
                    limit2 = 1
                endif
                
                clim   = clim*0+limit           ! Resetting the clim vector   
                IF(contact_alg .EQ. 1) THEN     ! Working contact routine     !D_IN: /Ref/ -> contact_alg
                    CAll    res(D2, D2p, MM, NN, J0, J1, JJ, J, kvot,Tcvot, DT2, SS, Pcenter)     !CALL: see arguments from code line
                ELSE
                    WRITE(4,*)'Warning! This contact_algorithm is not fully developed and is thus not working correctly. Use contact_alg=1 instead'
                    DO I=1+1*SS,MM,SS
                        I0=I-1*SS
                        I2=I-2*SS
                        IF (I2 .LE. 0) I2=I0
                        I1=I+1*SS
                        dHdx=(term4*H(I1,JJ)+term3*H(I,JJ)+term2*H(I0,JJ)+term1*H(I2,JJ))/(2*DX*SS)                                                           !D_IN: /Method/ -> term1,term2,term3,term4; /CurrentH/ -> H(I,J); /Grid/ -> Dx
                    
                        IF          (H(I,J) .GT. eps1) THEN
                            CAll    res1(I, D2, D2p, MM, NN, J0, J1, JJ, J, C1, C2, kvot, DT2, SS, Pcenter)       !CALL: See arguments from code line (This was not done)
                        Elseif      ((dhdx)   .GT. eps2) THEN
                            CALL    res2(I, D2, D2p, MM, NN, J0, J1, JJ, J, C1, C2, kvot, DT2, SS, Pcenter)       !CALL: no comment on the arguments (this was not done either)
                            clim(I) = 10*limit ! If contact, the pressure can increase clim times faster
                        else
                            CALL    res3(I, D2, D2p, MM, NN, J0, J1, JJ, J, C1, C2, kvot, DT2, SS, Pcenter)       !CALL: no comment on the arguments (this was not done either)
                            clim(I) = 10*limit ! If contact, the pressure can increase clim times faster
                        endif
                    
                    ENDDO
                ENDIF
                    
                ! Subrutien for adjusting the pressureupdates according to nearby pressureupdates in X-direction
                CALL TRA4(MM,SS)                                                                                  !CALL: N = MM, SS 
                
                DO I=1+SS,MM,SS
                    ! In the outlet region the pressure has a tendancy to diverge around the asperity. I thnik it has something to do with that the pressure goes from zero to something
                    lim_check=0
                    IF( I .GT. NX/2 .and. t .gt. Ntime*2.0/3.0) then                                              !D_IN: /Itt/ -> Ntime
                        i4=i+4
                        if ( i4 .GT. NX) i4=NX
                            IF( P(i4,j) .lt. 0.01 ) then !Then were at risky regions. 
                                clim(I) = limit/10
                                lim_check =  1
                            endif
                    endif

                    IF ( H(I,J) .LT. Hminimum)then
                        I0=I-1*SS
                        I1=I+1*SS

                        FC = 1.0
                        AK_sum = 1.0
                        IF (H(I0,J) .LT. Hminimum) AK_sum=AK_sum+AK_cvot 
                        IF (H(I1,J) .LT. Hminimum) AK_sum=AK_sum+AK_cvot 
                        IF (H(I,J0) .LT. Hminimum) AK_sum=AK_sum+AK_cvot 
                        IF (H(I,J1) .LT. Hminimum) AK_sum=AK_sum+AK_cvot 

                        incre=-FC/AK_sum*(H(I,J)-Hminimum)/(AK(0,0)*PAIAK*DX*SS)                                  !D_IN: /Ref/  -> PAIAK
                        IF( incre .LT. 0) then
                             WRITE(4,*)'check 1332 Bad incre in Single grid, at i,j,=', i,j
                             WRITE(4,*)'Incre was ', incre
                             incre=-incre
                             call Stop_to_large_out(t)
                        ENDIF
                        
                        IF( D(I) .LT. incre) D(I)=incre                                                           !D_OUT: D -> /CTRA4/
                        contact(i,j) = 2                                                                          !D_OUT: contact(I,J) -> /Contact_mat/
                    ENDIF
                    
                    ! limit the value of the increment since it's used for updating neibouring nodes as well. 
                    IF( abs(D(i)) .GT. 0.2) then
                        D(i) = sign(0.2, D(i))
                        Di_max = max(D(i) , Di_max)
                    endif

                ENDDO
                IF( C_Meth==0)then
                    sum_D=sum(abs(D))/(NX/SS)
                    if(sum_D .gt. 0.01) then
                        D=D*0.01/sum_D
                    endif
                endif

                DO 60 I=1+SS,MM,SS
                    if ( isnan(D(I) )) then
                        WRITE(4,*)'check 2222 Bad D at i,j,=', i,j
                        stop
                    endif
                    if( abs(D(I)) .gt. 1 .and. SS == 1) then
                        WRITE(4,*)'Warning, high D =',D(I),' at', i,j
                        if( abs(D(I)) .gt. 2) then
                            WRITE(4,*)'The components of A is:'
                            WRITE(4,*) A(:,I)                                                                     !D_IN: /CTRA4/ -> A
                        endif
                    endif
                    
                    IF(ID(I).EQ.2)GOTO 54                   ! If negativ lubrication height then skip the pressure uppdation      !D_IN: /ITERp/ -> ID(I)                         ! !!! Branch to target 54!
                    IF(ID(I).EQ.0)GOTO 52                   ! Choise of method                                                                                                    ! !!! Branch to target 52!
                    
                    ! Algorithm for ID=1
                    DD=D(I+1*SS)
                    IF(I.EQ.MM)DD=0                         ! If at the end in x-direction
                    
                    IF( C_meth .ge. 1) then
                        P_upd(I,J)                              = P_upd(I,J)  + C2*(D(I)-0.25*(D(I-1*SS)+DD))
                        IF(J0 .NE. 1)               P_upd(I,J0) = P_upd(I,J0) - 0.25*C2*D(I) ! Do not neet to update the pressure at P(I,1) since these are coppied later on. 
                        IF(P_upd(I,J0) .LT. 0.0)    P_upd(I,J0) = 0.0
                        
                        !IF(J == 1+ss)               P_upd(I,J)  = P_upd(I,J) - 0.25*C2*D(I) ! To account for the not performed investigation at J=1
                    else
                        P(I,J)                                  = P(I,J)+C2*(D(I)-0.25*(D(I-1*SS)+DD))
                        IF(J0 .NE. 1)               P(I,J0)     = P(I,J0)-0.25*C2*D(I) ! Do not neet to optade the pressure at P(I,1) since these are coppied later on. 
                        IF(P(I,J0) .LT. 0.0)        P(I,J0)     = 0.0
                    endif
                    
                    ! Update all nodes the same
                    ! IF(J1.GE.NN)GOTO 54                   ! Why treat middel node different than other nodes? so that middel node not based on uppdated values. The pressure at middel node is corrected in Res
                    IF( C_meth .ge. 1) Then
                        P_upd(I,J1)                             = P_upd(I,J1)-0.25*C2*D(I)        ! P(I,J1) is not updated
                        !IF(J1 == NN-SS)            P_upd(I,J1) = P_upd(I,J1)-0.25*C2*D(I)        ! To account for the update that would have been present if investigated the line J=NN+SS
                        IF(P_upd(I,J1) .LT. 0.0)    P_upd(I,J1) = 0.0
                    else
                        P(I,J1) = P(I,J1)-0.25*C2*D(I)        ! P(I,J1) is not updated
                        IF(P(I,J1) .LT. 0.0)    P(I,J1)=0.0
                    endif
                    
                    
                    GOTO 54
                    
                    ! Algorithme for ID=0
52                  IF(C_meth .ge. 1) then                                                                                                                                        ! !!! Branch target 52
                        P_upd(I,J)=P_upd(I,J)+C1*D(I)
                    else
                        P(I,J)=P(I,J)+C1*D(I)
                    endif
                    
                    ! Check so no negative pressures
54                  IF( C_meth .ge. 1 ) then                                                                                                                                      ! !!! Branch target 54
                        IF(P_upd(I,J).LT.0.0)   P_upd(I,J)=0.0
                    else
                        IF(P(I,J).LT.0.0)   P(I,J)=0.0
                    endif
                    
                    
60              CONTINUE
                
                IF( C_meth == 2 ) Then
                    CALL P_update_line(NYs, SS, k, C11, t, k_use, M_conv, Di_max, J, ID)                                                                                          !CALL: No comment on arguments 
                endif
                
70          CONTINUE
            
            IF( C_meth == 1) CALL P_update(NYs, SS, k, C11, t, k_use, M_conv, Di_max)                                                                                             !CALL: No comment on arguments


            ! Now when the pressure is updated, update the other parameters based on the new pressure
            CALL HREE(H00,t,SS, NYs, T40, Tcvot, k, k_use, M_conv)      ! The pressure will be updated if contact.And also if lite line is averaged                               !CALL: No comment on arguments
            K=K+1
            
    ENDDO
        
        IF(Geom .EQ. 5) THEN !Ball                                                                                                                                                !D_IN: /Ref/ -> Geom
            width=((NYs-1)/SS+1)*DX*SS
            G0=width*PAI/2 *(PH_new/PH)**2                                                                                                                                        !D_IN: /Geom5/ -> PH_new; /Visc/ -> PH;         !D_OUT: G0 -> /G0DT/
        ELSEIF(Geom .GE. 2) THEN !Ball
            CONTINUE
        ELSE
            ! Scale the load depending on the mesh size IF Cylinder. 
            width=((NYs-1)/SS+1)*DX*SS
            G0=width*PAI/2 
        ENDIF
        
        
        load=0.0
        DO J=1,NYs,SS
            DO I=1,NX,SS                        ! The load
                load=load+P(I,J)
            ENDDO
        ENDDO
        
        Sload=DX*DX*load/G0*SS*SS               ! Normalized load
        DW=Sload-1.0                            ! 1-scaled loadsum
        
             
        RETURN
    END