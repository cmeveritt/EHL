! Subroutine for calculating the elastic deformation
	SUBROUTINE VIside
        implicit none
        COMMON      /Grid/NX,NY,X0,XE,DX     
        COMMON      /COMAK2D/AK(0:300,0:300)
        COMMON      /Ref/ref, PAIAK, tmeth                                          ! Coise of equations based on referense
        COMMON      /Current/P,H,RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        real         P(1:300,1:300),H(1:300,1:300),RO(1:300,1:300),EPSx(1:300,1:300),EPSy(1:300,1:300),EDAx(1:300,1:300),EDAy(1:300,1:300),xi(1:300,1:300),W(1:300,1:300),Wside(1:300,1:300)                      ! Current timestep
        real        AK
        integer     ref, I,J,K,L, IK,JL, NX,NY, tmeth
        real        PAIAK, H0, DX, X0, XE
        SAVE        /Current/

        
        DO 40 I=1,NX
            DO 40 J=1,NY
                H0=0.0
                !Defromation from pressure outside the area. All nodes deforms from pressure 4*2*NX number of nodes in X-directon. 
                    !For the nodes on the -y side of the contact
                    DO K=1,NX              
                        IK=IABS(I-K)
                        DO L=1,NX+1-J
                            JL=J+L-1
                            H0=H0+AK(IK,JL)*P(K,3)
                           
                        END DO
      
                        !For the nodes on the +y side of the contact 
                        DO L=1,NX+1+J-NY-1
                            JL=L+NY-J
                            H0=H0+AK(IK,JL)*P(K,3)
                           
                        END DO
                    END DO
                

        40      Wside(I,J)=H0*DX*PAIAK  ! See Eq 2.12 in the new book. The scaling factor between cylinder änd ball is pi/4. V_cyl=pi/4*V_ball due to the scaling R/a^2
        
        
        RETURN
    END