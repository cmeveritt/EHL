! A subroutine to calculate the increments to the residual from the Poiseuille term in the Reynlds equation based on a pressure, p, change
    !d(eps(i,j))/d(P(k,l))=incre
    ! The subroutine was initially setup with all the derivatives. This was then discarded since no increase in stability nor convergence  was obtained 
    ! Therefore is so much of the code commented out. Theoretically it should be a good ide to use all derivatives so it could be a good idea to look into in the future. 
    Subroutine Poiseuille_increment(i,j,i0,j0,i1,j1,i2, K00, K10, K20, K30, D1, D2, D3, control, Tcvot, ss, AII1dxx, AII2dxx, AII3dxx, AII4dxx, dRodP, d2ROdP2)
    implicit none
    include     'inc_AKparam.h'
    include     'inc_Contact_mat.h'
    include     'inc_CurrentP.h'
    include     'inc_Grid.h'
    !Input
    integer     i,j,i0,j0,i1,j1,i2, control
    real        K00, K10, K20, K30, D1, D2, D3
    real        Tcvot
    !Calculations
    integer     k, l, ii, jj, m, ss
    real        incre_i, incre_im1, incre_jm1, incre_ip1, incre_jp1
    real        AII(1:4)
    real        test
    !Output
    real        AII1dxx, AII2dxx, AII3dxx, AII4dxx
    real        dRodP, d2ROdP2
        
    
    ! To save calculations (DX*SS)**2 is added later as DX3
    
    if (control==1) then
        AII1dxx =   0
        AII2dxx =   D1+0.25*D3
        AII3dxx =   -1.25*D3
        AII4dxx =   D2+0.25*D3
        K00     =   BK00
        K10     =   BK10
        K20     =   BK20
        K30     =   BK30
    else
        AII1dxx =   0
        AII2dxx =   D1
        AII3dxx =   -D3
        AII4dxx =   D2
        K00     =   AK00
        K10     =   AK10
        K20     =   AK20
        K30     =   AK30
    endif
    
    d2ROdP2=0
    dRodP=0
    
    !If( i2 .NE. i0 .and. i1 .NE. i) then ! Do not take this derivatives into acocunt if to close to the border
    !    do m=1,4
    !        if(m==1) then
    !            k=i2
    !        elseif(m==2)then
    !            k=i0
    !        elseif(m==3)then
                k=i
    !        elseif(m==4)then
    !            k=i1
    !        endif
    !    
            l=j
            jj=j
    !        ii=i0
    !        CALL Epsilon_derivative(ii,jj,k,l, K00, K10, K20, K30,Tcvot, ss, incre_im1, test, test)
            ii=i
            CALL Epsilon_derivative(ii,jj,k,l, K00, K10, K20, K30,Tcvot, ss, incre_i, dRodP, d2ROdP2)  !CALL: incre_i is not initialized yet and is probably used to read an output from this function. Apparently incre dRodP and D2ROdP2 are outputs
    !        ii=i1
    !        CALL Epsilon_derivative(ii,jj,k,l, K00, K10, K20, K30,Tcvot, ss, incre_ip1, test, test) 
    !
    !        ii=i
    !        jj=j0
    !        CALL Epsilon_derivative(ii,jj,k,l, K00, K10, K20, K30,Tcvot, ss, incre_jm1, test, test)
    !        jj=j1
    !        CALL Epsilon_derivative(ii,jj,k,l, K00, K10, K20, K30,Tcvot, ss, incre_jp1, test, test)
    !        !incre_i=0
    !        !incre_ip1=0
    !        !incre_im1=0
    !        
    !        !incre_jm1=0
    !        !incre_jp1=0
    !        ! d(eps)/dx
    !        AII(m)  = ((incre_im1+incre_i)*P(i0,j) - (2*incre_i+(incre_im1+incre_ip1))*P(i,j) + (incre_ip1+incre_i)*P(i1,j))*0.5
    !    
    !        !d(eps)/dy
    !        AII(m)  = AII(m) + ((incre_jm1+incre_i)*P(i,j0) - (2*incre_i+(incre_jm1+incre_jp1))*P(i,j) + (incre_jp1+incre_i)*P(i,j1))*0.5
    !        
    !        !AII(m) = AII(m)*0.10
    !        
    !        if ( isnan( AII(m))  ) then
    !            WRITE(4,*)' check P_inc 2234 Bad AII(m). Therfore its reset to 0, at II=', II
    !            Call Stop_to_large_out
    !            AII(m)=0.0
    !        endif 
    !        
    !    enddo
    !endif
    
    
    !if( i .GT. 1 .and. j .gt. 1) then
    !    if(contact(i,j)==1 .or. contact(i+ss,j)==1 .or. contact(i-ss,j)==1  .or. contact(i,j-ss)==1 .or. contact(i,j+ss)==1 )then
            AII1dxx  = AII1dxx
            AII2dxx  = AII2dxx
            AII3dxx  = AII3dxx
            AII4dxx  = AII4dxx
    !    else
    !        AII1dxx  = 0.5*AII(1)+AII1dxx
    !        AII2dxx  = 0.5*AII(2)+AII2dxx
    !        AII3dxx  = 0.5*AII(3)+AII3dxx
    !        AII4dxx  = 0.5*AII(4)+AII4dxx
    !    endif
    !endif                                            
                
                
    return
    end
    