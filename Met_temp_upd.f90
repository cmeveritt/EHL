subroutine Met_temp_upd(dT_lim, nys, SS, t0, time)
    ! A subroutine to calculate the temperature inside the fluid film
    ! As a start based on isothermal surfaces 
    implicit none
    ! Input
    include     'inc_CurrentT.h'
	include     'inc_CurrentT_met.h' 
    include     'inc_Grid.h'
    include     'inc_Met_Temp_inc.h'
    include     'inc_Cores_used.h'
    integer Nys, ss
    real    dt_lim, t0, time
    ! Calculations
    integer i,j,k, jj, n_met
    integer NX_m, NY_m, NN_m
    real    test1, test2, test3, test4, test5
    ! Output
    save  /CurrentT_met/
    
    NX_m=(nx+1)/2      !D_IN: /Grid/ -> Nx
    NY_m=(nys+1)/2
    NN_m=(NY_m+1)/2
    n_met=39           !This has the same value in every metal nodes subroutine. Should be placed in an inputfile to increas stability

    ! If in the partial time step for temperature, only update top layer
    IF( dt_lim .LT. 0 .or. time .lt. 0) then
        WRITE(4,*)'Negative time increment'
        WRITE(4,*)'dt_lim = ', dt_lim, ',  time = ', time
        STOP
    endif
    
    
    if( time==0) then
        k=1
        !$OMP PARALLEL DO SHARED(SS,NN_m,NX_m,Temp_ma,Temp_mb,temp_incre,t0,NY_m,dt_lim,k)&
        !$OMP&            PRIVATE(J,JJ,I) &
        !$OMP&  IF(use_multiple_cores)
        do j=1,NN_m,SS
            JJ=NY_m-J+1
            IF( JJ .LE. 1)  JJ=1+SS
            do i=1,NX_m,SS
                    
                Temp_ma(i,j,k)=Temp_ma(i,j,k)+temp_incre(i,j,k,1)*dt_lim   !D_IN: /Met_temp_inc/ -> temp_incre(I,J) !D_OUT: Temp_ma(I,J) -> /CurrentT_met/;
                Temp_mb(i,j,k)=Temp_mb(i,j,k)+temp_incre(i,j,k,2)*dt_lim   !D_OUT: Temp_mb(I,J) -> /CurrentT_met/
                
                if( Temp_ma(i,j,k) .GT. 1200) Temp_ma(i,j,k)=1200
                if( Temp_mb(i,j,k) .GT. 1200) Temp_mb(i,j,k)=1200
                
                ! Validate the the metals downt cool of
                if( Temp_ma(i,j,k) .LT. t0) then 
                    if( Temp_ma(i,j,k) .LT. t0-2 ) WRITE(4,*)'Warning! Something wrong with metal temperature' ! Some numerics can reduce the temperature slightly below the initial temp
                    Temp_ma(i,j,k)=t0
                endif
                
                if( Temp_mb(i,j,k) .LT. t0) then
                    if( Temp_ma(i,j,k) .LT. t0-2 ) WRITE(4,*)'Warning! Something wrong with metal temperature' ! Some numerics can reduce the temperature slightly below the initial temp
                    Temp_mb(i,j,k) = t0               
                endif
                Temp_ma(i,jj,k)=Temp_ma(i,j,k)
                Temp_mb(i,jj,k)=Temp_mb(i,j,k)
            enddo
        enddo
        !$OMP END PARALLEL DO
    
    ! If finished with the partial time simulations then update down in the material
    else
        !$OMP PARALLEL DO SHARED(SS,N_met,NN_m,NX_m,Temp_ma,Temp_mb,temp_incre,time,t0,NY_m)&
        !$OMP&            PRIVATE(K,J,JJ,I) &
        !$OMP&  IF(use_multiple_cores)
        do k=1+ss,n_met,ss
            do j=1,NN_m,SS
                JJ=NY_m-J+1
                IF( JJ .LE. 1)  JJ=1+SS
                do i=1,NX_m,SS
                  
                    Temp_ma(i,j,k)=Temp_ma(i,j,k)+temp_incre(i,j,k,1)*time
                    Temp_mb(i,j,k)=Temp_mb(i,j,k)+temp_incre(i,j,k,2)*time
                
                    if( Temp_ma(i,j,k) .GT. 1200) Temp_ma(i,j,k)=1200
                    if( Temp_mb(i,j,k) .GT. 1200) Temp_mb(i,j,k)=1200
                
                    if( Temp_ma(i,j,k) .LT. t0) then
                        if( Temp_ma(i,j,k) .LT. t0-2 ) WRITE(4,*)'Warning! Something wrong with metal temperature' ! Some numerics can reduce the temperature slightly below the initial temp
                        Temp_ma(i,j,k)=t0
                    endif
                
                    if( Temp_mb(i,j,k) .LT. t0) then
                        if( Temp_ma(i,j,k) .LT. t0-2 ) WRITE(4,*)'Warning! Something wrong with metal temperature' ! Some numerics can reduce the temperature slightly below the initial temp
                        Temp_mb(i,j,k)=t0
                    endif
                
                    Temp_ma(i,jj,k)=Temp_ma(i,j,k)
                    Temp_mb(i,jj,k)=Temp_mb(i,j,k)
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
    endif
    
    return
    end
    
    