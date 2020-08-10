! Subroutine for calculating the elastic deformation based upon the fast fourier transforms
! This is much faster than the VI subroutine. For referenses:
! Liu et al. 2000 A versatile method of discrete convolution and FFT (DC-FFT) for contact analyses 
! Wang et al. 2003 A comparative study of the methods for calculation of surface elastic deformation
include "mkl_dfti.f90"
SUBROUTINE VIFFT(SS, NYs, NN, t)
    use MKL_DFTI
    implicit none
    include     'inc_COMAK2D.h'
    include     'inc_Contact_mat.h'
    include     'inc_Current.h'
	include     'inc_CurrentP.h'
    include     'inc_Grid.h'
    include     'inc_Itt.h'
    include     'inc_Ref.h'
    include     'inc_Pres_ave_param.h'
    include     'inc_P_line_side.h' 
    include     'inc_Cores_used.h'
    include     'inc_Output_control.h' 
    ! Input
    integer     SS, NYs, NN, iter, t, full
    ! Calculation parameters
    integer     K, L, I, J, IK, JL, JJ, ii
    integer     Js, Je, number_n
    integer     LL, JL1, JL2
    real        H0
    real, allocatable :: P_fft_real(:,:),W_fft_real(:,:),AK_fft_real(:,:),P_extension(:,:)
    complex, allocatable :: P_fft_cmplx(:,:),W_fft_cmplx(:,:),AK_fft_cmplx(:,:)
    integer                 cstrides(3), rstrides(3)
    Integer                 status, dimensions(2)
    real                 scale
    type(DFTI_DESCRIPTOR), pointer :: My_desc1_handle
    
    real ::               W_test(NX,NYs), err_array(NX,NYs), dummy(NX,NYs)
    ! Output
    save        /current/

    
    
    IF(Geom .EQ. 2 .OR. Geom .EQ. 3 .OR. Geom .EQ. 6) THEN                                     !D_IN: /Ref/ -> Geom
                                                 
        Allocate(P_fft_real(1:2*(NX/SS)-1+10,1:2*(NYs/SS)-1+3), W_fft_real(1:2*(NX/SS)-1+10,1:2*(NYs/SS)-1+3), AK_fft_real(1:2*(NX/SS)-1+10,1:2*(NYs/SS)-1+3))  !Allocate the arrays used in the calculation. ...
                                                                                                                                                                !The size is to accomodate the expanded AK matrix. ...
                                                                                                                                                                !The added 10 rows and 3 columns of zeros are to give the solution the space to appear and was found through trial and error.
        Allocate(P_fft_cmplx((2*(NX/SS)-1+10)/2 + 1, 2*(NYs/SS)-1 + 3), W_fft_cmplx((2*(NX/SS)-1+10)/2 + 1, 2*(NYs/SS)-1 + 3), AK_fft_cmplx((2*(NX/SS)-1+10)/2 + 1, 2*(NYs/SS)-1+3))
        rstrides = [0, 1, 2*(NX/SS)-1+10]
        cstrides = [0, 1, (2*(NX/SS)-1+10)/2 + 1]
        dimensions = (/2*(NX/SS)-1+10,2*(NYs/SS)-1+3/)
        
        AK_fft_real = 0
        AK_fft_real(10+1:2*(NX/SS)-1+10,1+3:2*(NYs/SS)-1+3) = AK((/ (Nx/SS)-1:1:-1 , 0:(Nx/SS)-1 /),(/ (Nys/SS)-1:1:-1 , 0:(Nys/SS)-1 /))
        P_fft_real = 0
        P_fft_real((Nx/SS)+10:2*(Nx/SS)-1+10,(NYs/SS)+3:2*(Nys/SS)-1+3) =  P((/ 1:NX:SS /) , (/ 1:NN:SS , NN-1:1:-1*SS /))
        
        status = DftiCreateDescriptor(My_desc1_handle, DFTI_SINGLE, DFTI_REAL,2,dimensions)         !Single 2D forward transformation 
        status = DftiSetValue(My_desc1_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)                    !Not in place produces the result in another array variable
        status = DftiSetValue(My_desc1_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)   !This works but if you want you can read more about the different complex storage scheme in MKL here: https://software.intel.com/en-us/mkl-developer-reference-c-dfti-complex-storage-dfti-real-storage-dfti-conjugate-even-storage
        status = DftiSetValue(My_desc1_handle, DFTI_INPUT_STRIDES, rstrides)                        !Strides are the distances between references and data, successive elements and successive rows.
        status = DftiSetValue(My_desc1_handle, DFTI_OUTPUT_STRIDES, cstrides)                       !The strides are set for both inputs and outputs
        status = DftiCommitDescriptor(My_desc1_handle)                                              !Creating and commiting the descriptor is necessary before FFT operations
    

        status = DftiComputeForward(My_desc1_handle, AK_fft_real(:,1), AK_fft_cmplx(:,1))           !Compute the forward transform of the deformation matrix
        status = DftiComputeForward(My_desc1_handle, P_fft_real(:,1), P_fft_cmplx(:,1))             !Compute the forward transform of the pressure profile

    
        W_fft_cmplx = AK_fft_cmplx * P_fft_cmplx                                        
        scale = 1.0/(dimensions(1)*dimensions(2))                                              

        status = DftiSetValue(My_desc1_handle, DFTI_INPUT_STRIDES, cstrides)                        !The strides of the output of the forward transformation are the strides of the input to the inverse transformation.
        status = DftiSetValue(My_desc1_handle, DFTI_OUTPUT_STRIDES, rstrides)                       !Same here. The strides of the input to the forward transformation are the strides of the output here.
        status = DftiSetValue(My_desc1_handle, DFTI_BACKWARD_SCALE, scale)                          !The scale equals the size of the arrays and needed to make the IFFT truly the inverse of the forwards FFT.
        status = DftiCommitDescriptor(My_desc1_handle)
        status = DftiComputeBackward(My_desc1_handle, W_fft_cmplx(:,1), W_fft_real(:,1))
        status = DftiFreeDescriptor(My_desc1_handle)

        W(1:NX:SS,1:NYs:SS) = W_fft_real(10:(NX/SS)+9,1+2:(Nys/SS)+2) * DX * SS * PAIAK             !The output is picked out from the expanded result. This index is exactly where the correct output shows up given how the problem is setup in this case.
        

    ELSE 
        Allocate(P_fft_real(1:2*(NX/SS)-1+10,1:2*(2*(NX/SS)+(NYs/SS))-1+1), W_fft_real(1:2*(NX/SS)-1+10,1:2*(2*(NX/SS)+(NYs/SS))-1+1), AK_fft_real(1:2*(NX/SS)-1+10,1:2*(2*(NX/SS)+(NYs/SS))-1+1)) !The only difference here is that the pressure profile is expanded first which makes matching the sizes harder. 
        Allocate(P_extension(1:(NX/SS),1:(NX/SS)))  !This stores p_line repeated to create a matrix.
        Allocate(P_fft_cmplx((2*(NX/SS)-1+10)/2 + 1, 2*(2*(NX/SS)+(NYs/SS))-1+1), W_fft_cmplx((2*(NX/SS)-1+10)/2 + 1, 2*(2*(NX/SS)+(NYs/SS))-1+1), AK_fft_cmplx((2*(NX/SS)-1+10)/2 + 1, 2*(2*(NX/SS)+(NYs/SS))-1+1))
        rstrides = [0, 1, 2*(NX/SS)-1+10]
        cstrides = [0, 1, (2*(NX/SS)-1+10)/2 + 1]
        dimensions = (/2*(NX/SS)-1+10,2*(2*(NX/SS)+(NYs/SS))-1+1/)
        
        AK_fft_real = 0
        AK_fft_real(10+1:2*(NX/SS)-1+10,(NX/SS)+(NYs/SS)+1:(NYs/SS)+3*(NX/SS)-1) = AK((/ (NX/SS)-1:1:-1 , 0:(NX/SS)-1 /),(/ (NX/SS)-1:1:-1 , 0:(NX/SS)-1 /))
        
        do j = 1,(NX/SS)
            k = 1
            do i = 1,(NX/SS)
                P_extension(i,j) = P_line(k) !the pressure extenstion is the the line pressure profile, outside the simulated area, expanded to an array.                
                k = k + SS
            end do
        end do
        
        P_fft_real = 0
        P_fft_real((NX/SS)+10:2*(NX/SS)-1+10, 2*(NX/SS)+(NYs/SS):3*(NX/SS)+(NYs/SS)-1) = P_extension
        P_fft_real((NX/SS)+10:2*(NX/SS)-1+10, 3*(NX/SS)+(NYs/SS):3*(NX/SS)+2*(NYs/SS)-1) = P((/ 1:NX:SS /) , (/ 1:NN:SS , NN-1:1:-1*SS /)) !P(1:NX:SS,1:NYs:SS) 
        P_fft_real((NX/SS)+10:2*(NX/SS)-1+10, 3*(NX/SS)+2*(NYs/SS):2*(2*(NX/SS)+(NYs/SS))-1) = P_extension

        status = DftiCreateDescriptor(My_desc1_handle, DFTI_SINGLE, DFTI_REAL,2,dimensions)
        status = DftiSetValue(My_desc1_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        status = DftiSetValue(My_desc1_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
        status = DftiSetValue(My_desc1_handle, DFTI_INPUT_STRIDES, rstrides)
        status = DftiSetValue(My_desc1_handle, DFTI_OUTPUT_STRIDES, cstrides)
        status = DftiCommitDescriptor(My_desc1_handle)

        status = DftiComputeForward(My_desc1_handle, AK_fft_real(:,1), AK_fft_cmplx(:,1))
        status = DftiComputeForward(My_desc1_handle, P_fft_real(:,1), P_fft_cmplx(:,1))

        W_fft_cmplx = AK_fft_cmplx * P_fft_cmplx                     
        scale = 1.0/(dimensions(1)*dimensions(2))

        status = DftiSetValue(My_desc1_handle, DFTI_INPUT_STRIDES, cstrides)
        status = DftiSetValue(My_desc1_handle, DFTI_OUTPUT_STRIDES, rstrides)
        status = DftiSetValue(My_desc1_handle, DFTI_BACKWARD_SCALE, scale)
        status = DftiCommitDescriptor(My_desc1_handle)
        status = DftiComputeBackward(My_desc1_handle, W_fft_cmplx(:,1), W_fft_real(:,1))
        status = DftiFreeDescriptor(My_desc1_handle)
        W(1:NX:SS,1:NYs:SS) = W_fft_real(10:(NX/SS)+9,(NX/SS)-1:(NX/SS)+(Nys/SS)-2) * DX * SS * PAIAK


    ENDIF

    
110     FORMAT(2001(E13.6,1X))        

    END
    