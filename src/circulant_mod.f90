module circulant_mod
    use iso_fortran_env
    use structblas_util_mod
    implicit none

    private

    type :: dcirculant
        real(kind=real64), dimension(:), allocatable :: c
        complex(kind=real64), dimension(:), allocatable :: lambda
        integer(kind=int32) :: n = -1
        logical :: isbuild = .false.
        contains
        procedure, pass(circ) :: setc => sb_d_circulant_setc
        procedure, pass(circ) :: setlambda => sb_d_circulant_setlambda
        procedure, pass(circ) :: build => sb_d_circulant_build
        procedure, pass(circ) :: tomat => sb_d_circulant_tomat
        procedure, pass(circ) :: free => sb_d_circulant_free
    end type dcirculant

    !! d submodule interfaces
    interface
        module subroutine sb_d_circulant_setc(circ,c,info,n)
            use iso_fortran_env
            use structblas_const_mod
            implicit none
            class(dcirculant), intent(inout) :: circ !! circulant matrix object
            real(kind=real64), intent(in), dimension(:) :: c !! input circulant coefficients
            integer(kind=int32), intent(out) :: info !! success of the routine
            integer(kind=int32), intent(in), optional :: n !! size of the circulant matrix
        end subroutine sb_d_circulant_setc
    end interface

    interface
        module subroutine sb_d_circulant_setlambda(circ,lambda,info)
            !! Implementation that sets the entries of the first column of a circulant matrix
            use iso_fortran_env
            use structblas_const_mod
            implicit none
            class(dcirculant), intent(inout) :: circ !! circulant matrix object
            complex(kind=real64), intent(in), dimension(:) :: lambda !! input eigenvalues of the circulant matrix
            integer(kind=int32), intent(out) :: info !! success of the routine
        end subroutine sb_d_circulant_setlambda
    end interface

    interface
        module subroutine sb_d_circulant_build(circ,info)
            !! If the circulant matrix has c and not lambda builds the latter, if it has lambda 
            !! and not c builds the former
            use iso_fortran_env
            use iso_c_binding
            use structblas_const_mod
            implicit none
            class(dcirculant), intent(inout) :: circ !! circulant matrix object
            integer(kind=int32), intent(out) :: info !! success of the routine
        end subroutine sb_d_circulant_build
    end interface

    interface
        module subroutine sb_d_circulant_tomat(circ,mat,info)
            !! Assembles the circulant matrix in the "full" representation
            use iso_fortran_env
            use structblas_const_mod
            implicit none
            class(dcirculant), intent(in) :: circ !! circulant matrix object
            real(kind=real64), allocatable, dimension(:,:), intent(inout) :: mat !! matrix to be filled
            integer(kind=int32), intent(out) :: info !! success of the routine
        end subroutine sb_d_circulant_tomat
    end interface

    interface 
        module subroutine sb_d_circulant_free(circ,info)
            !! Free memory occupated by circulant matrix
            use iso_fortran_env
            use structblas_const_mod
            implicit none
            class(dcirculant), intent(inout) :: circ !! circulant matrix object
            integer(kind=int32), intent(out) :: info !! success of the routine
        end subroutine sb_d_circulant_free
    end interface

    public :: dcirculant, d_circ_gemv

contains

        subroutine d_circ_gemv(trans,alpha,circm,x,beta,y,planf,planb)
            !! Computes \( \mathbf{y} = \alpha C \mathbf{x} + \beta \mathbf{y} \) or \( \mathbf{y} = \alpha C^T \mathbf{x} + \beta \mathbf{y} \)
            use iso_fortran_env
            use iso_c_binding
            use structblas_const_mod
            use structblas_util_mod
            use fftw_mod
            implicit none
            character(len=1), intent(in) :: trans
            real(kind=real64), intent(in) :: alpha
            type(dcirculant), intent(inout) :: circm
            real(kind=real64), intent(inout), dimension(:) :: x
            real(kind=real64), intent(in) :: beta
            real(kind=real64), intent(inout), dimension(:) :: y
            type(c_ptr), intent(inout), optional :: planf,planb

            ! Local variables
            integer(kind=int32) :: info,n,i
            complex(kind=c_double_complex), dimension(:), allocatable :: forwardvec
            real(kind=c_double), dimension(:), allocatable :: backvec

            logical :: wehaveaplan
            type(c_ptr) :: planf_,planb_

            ! Check inputs
            if (.not.circm%isbuild) then
                if (allocated(circm%c).or.allocated(circm%lambda)) then
                    call circm%build(info)
                else
                    write(error_unit,'("d_circ_genmv: circm is empty")')
                    return
                end if
            else
                n = circm%n
                if ((n /= size(x)).or.(n /= size(y))) then
                    write(error_unit,'("d_circ_genmv: circm, x, y sizes are not compatible")')
                    return
                end if
            end if
            if (present(planf).and.present(planb)) then
                !!
                !! Passing a fftw plan to this routine is useful for repeated matrix-vector products
                !! if you pass it, and the plan is not associated to something after the first
                !! matrix vector product the plan is created and can be reused
                !!
                if (c_associated(planf).and.c_associated(planb)) then
                    planf_ = planf
                    planb_ = planb
                    wehaveaplan = .true.
                else
                    wehaveaplan = .false.
                end if
            else
                wehaveaplan = .false.
            end if

            !! 
            !! The function checks in case of quick returns: \( n = 0 \) or \( \alpha =0 \land \beta = 0 \).
            !!
            if ((n.eq.0).or.((alpha.eq.dzero).and.(beta.eq.done))) return

            allocate(forwardvec(n/2+1),stat=info)
            allocate(backvec(n),stat=info)

            if ((beta /= done).and.(beta /= dzero)) y = beta*y
            select case (sb_toupper(trans))
            case ('T')
                !!
                !! If `trans = 'T'` then \( \mathbf{y} = \alpha C^T \mathbf{x} + \beta \mathbf{y} \)
                !!

                ! if (.not.wehaveaplan) then
                !     !! We don't have a plan and we have to generate them
                !     planf_ = fftw_plan_dft_r2c_1d(n,x,forwardvec,FFTW_MEASURE)
                !     planb_ = fftw_plan_dft_r2c_1d(n,forwardvec,backvec,FFTW_MEASURE)
                ! end if
                ! !! We have FFTWs plan and we can just perform multiplications
                ! call fftw_execute_dft_r2c(planf_,x,forwardvec)
                ! do i=1,n/2+1
                !     forwardvec(i) = alpha*forwardvec(i)*circm%lambda(i)
                ! end do
                ! call fftw_execute_dft_c2r(planb_,forwardvec,backvec)
                ! y = backvec + x        
            case ('N')
                !!
                !! If `trans = 'N'` then \( \mathbf{y} = \alpha C \mathbf{x} + \beta \mathbf{y} \)
                !!
                if (.not.wehaveaplan) then
                    ! We don't have a plan and we have to generate them
                    planf_ = fftw_plan_dft_r2c_1d(n,backvec,forwardvec,FFTW_MEASURE)
                    planb_ = fftw_plan_dft_c2r_1d(n,forwardvec,backvec,FFTW_MEASURE)
                end if
                ! We have FFTWs plan and we can just perform multiplications
                call fftw_execute_dft_r2c(planf_,x,forwardvec)
                ! write(output_unit,*) "fft(x) = ", forwardvec
                do i=1,n/2+1
                    forwardvec(i) = forwardvec(i)*circm%lambda(i)
                end do
                call fftw_execute_dft_c2r(planb_,forwardvec,backvec)
                y = alpha*backvec/real(n) + y
            case default
                write(error_unit,'("d_circ_genmv: unknown transpose/non-transpose charachter")')
                deallocate(forwardvec,stat=info)
                deallocate(backvec,stat=info)
                return
            end select

            deallocate(forwardvec,stat=info)
            deallocate(backvec,stat=info)

        end subroutine
    
end module circulant_mod
