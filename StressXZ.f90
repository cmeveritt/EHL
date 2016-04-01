! Subrutine to calculate the stresses  in the x-z plane at the centerline of the contact
    SUBROUTINE StressXZ
        implicit none
        COMMON      /Grid/NX,NY,X0,XE,DX 
        COMMON      /Stress/sig_x, sig_y, sig_z
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep                      
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real         sig_x(300), sig_y(300), sig_z(300)
        real         Pmin(NX), Pside(NX)
        real         X0, XE, DX
        integer      NX, NY, NN
        integer      i, j, k, l
        real         ik, jl, ikk, jll
        real         abs, min, pi
        save        /Stress/
        
    pi=3.14159265
    NN=(NY+1)/2
    ! Finds the minimum lineload for each line
    DO i=1,NX
        Pmin(i)=minval(P(i,1:NY))
        sig_x(i)=0.0
        sig_y(i)=0.0
        sig_z(i)=0.0
        Pside(i)=P(i,2)
    enddo
    
    
    ! Calculats the stresses from the pressure at the area
    DO i=2,NX-1
        j=NN            ! Only interested in the centerline
            
                        !  from the pressure at the area
        Do l=1,NY
            Do k=1,NX
                ik=abs(i-k)*DX;
                jl=abs(j-l)*DX;
                ikk=(k-i)*DX;
                jll=(l-j)*DX;
                    
                ! If we're at the nod where the pressure is applied
                if ( ik .EQ. 0 .AND. jl .EQ. 0) THEN
                    sig_x(i)=sig_x(i)-(P(k,l)-Pmin(k))/(2)*(1+2*0.3);
                    sig_y(i)=sig_y(i)-(P(k,l)-Pmin(k))/(2)*(1+2*0.3);
                        
                else
                    sig_x(i)=sig_x(i)+(P(k,l)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(ik**2-jl**2)/(ik**2+jl**2)**2;
                    sig_y(i)=sig_y(i)+(P(k,l)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(jl**2-ik**2)/(ik**2+jl**2)**2;
                    
                    ! No shear stresses due to symmetry
                    ! tau_xy(i,j,t)=tau_xy(i,j,t)+(P(k,l,t)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(2*ikk*jll)/(ik**2+jl**2)**2;
                endif
            enddo
        enddo
            
        ! stresses from pressure outside the area
            
        !For the nodes on the -y side of the contact
        DO k=1,NX
            ik=abs(i-k)*DX;
            ikk=(k-i)*DX;
            do l=1,NX+1-j
                jl=(j+l-1)*DX;
                jll=-jl;
                
                sig_x(i)=sig_x(i)+(Pside(k)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(ik**2-jl**2)/(ik**2+jl**2)**2;
                sig_y(i)=sig_y(i)+(Pside(k)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(jl**2-ik**2)/(ik**2+jl**2)**2;
                !tau_xy(i,j,t)=tau_xy(i,j,t)+(Pside(k,l,t)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(2*ikk*jll)/(ik**2+jl**2)**2;
                    
            enddo
                
            !For the nodes on the +y side of the contact
            do l=1,NX+1+j-NY-1
                jl=(l+NY-j)*DX;
                jll=jl;
                sig_x(i)=sig_x(i)+(Pside(k)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(ik**2-jl**2)/(ik**2+jl**2)**2;
                sig_y(i)=sig_y(i)+(Pside(k)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(jl**2-ik**2)/(ik**2+jl**2)**2;
                !tau_xy(i,j,t)=tau_xy(i,j,t)+(Pside(k,l,t)-Pmin(k))/(2*pi)*(1-2*0.3)*DX**2*(2*ikk*jll)/(ik**2+jl**2)**2;
            enddo
        enddo
            
            !%-----------------------------------------
            !% Add shear stresses acc Hamrock page 465  Not for pure rolling
            !%----------------------------------------
            
            !EDA(i,j,t)=EDA0*exp(alpha*Pref/Z*(-1+(1+P(i,j,t)*10**6/Pref)**Z));   %Pref is in Pa

            
            !tau_xz(i,j,t)=-H(i,j,t)*a**2/(2*R)*(P(i+1,j,t)-P(i-1,j,t))/(2*DX*a)+(Ub-Ua)/H(i,j,t)*R/a**2*EDA(i,j,t)*10**-6; %Converted to MPa
            !tau_yz(i,j,t)=-H(i,j,t)*a**2/(2*R)*(P(i,j+1,t)-P(i,j-1,t))/(2*DX*a);
           
            !%Limiting shear value
            !if tau_xz(i,j,t)>500
            !    tau_xz(i,j,t)=500;
            !end
            
            !if tau_yz(i,j,t)>500
            !   tau_yz(i,j,t)=500;
            !end
            
            sig_z(i)=-P(i,j);
            
            ! Effective stress
            !sig_vm(i,j,t)=1/sqrt(2)*sqrt((sig(1)-sig(2))**2+(sig(1)-sig(3))**2+(sig(2)-sig(3))**2);
            
    enddo
end
    