submodule (circulant_mod) dcirculant_mod
    use fftw_mod
    implicit none
    
contains
    
    subroutine sb_d_circulant_setc(circ,c,info,n)
        !! Implementation that sets the entries of the first column of a circulant matrix
        use iso_fortran_env
        use structblas_const_mod
        implicit none
        class(dcirculant), intent(inout) :: circ !! circulant matrix object
        real(kind=real64), intent(in), dimension(:) :: c !! input circulant coefficients
        integer(kind=int32), intent(out) :: info !! success of the routine
        integer(kind=int32), intent(in), optional :: n !! size of the circulant matrix

        !! Local variables
        integer(kind=int32) :: n_ !! size of the circulant matrix
        integer(kind=int32) :: nonzero !! number of nonzero elements of the circulant matrix

        info = structblas_success_
        
        if (present(n)) then
            n_ = n
        else
            n_ = size(c)
        end if

        if ( (n_.ne.circ%n).and.(allocated(circ%lambda))) then
            deallocate(circ%lambda, stat=info)
        end if

        if (allocated(circ%c)) then
            deallocate(circ%c, stat=info)           
        end if

        allocate(circ%c(n_), stat=info, source=dzero)

        nonzero = min(n_,size(c))

        circ%c(1:nonzero) = c(1:nonzero) 
        circ%n = n_

    end subroutine sb_d_circulant_setc

    subroutine sb_d_circulant_setlambda(circ,lambda,info)
        !! Implementation that sets the entries of the first column of a circulant matrix
        use iso_fortran_env
        use structblas_const_mod
        implicit none
        class(dcirculant), intent(inout) :: circ !! circulant matrix object
        complex(kind=real64), intent(in), dimension(:) :: lambda !! input eigenvalues of the circulant matrix
        integer(kind=int32), intent(out) :: info !! success of the routine
        
        !! Local variables
        integer(kind=int32) :: n_ !! size of the circulant matrix

        info = structblas_success_
        n_ = size(lambda)
        
        if ( (n_.ne.circ%n).and.(allocated(circ%c))) then
            deallocate(circ%c, stat=info)
        end if

        if (allocated(circ%lambda)) then
            deallocate(circ%lambda, stat=info)           
        end if

        allocate(circ%lambda(n_), stat=info)

        circ%lambda = lambda 
        circ%n = n_

    end subroutine sb_d_circulant_setlambda

    subroutine sb_d_circulant_build(circ,info)
        !! If the circulant matrix has c and not lambda builds the latter, if it has lambda 
        !! and not c builds the former
        use iso_fortran_env
        use fftw_mod
        use structblas_const_mod
        implicit none
        class(dcirculant), intent(inout) :: circ !! circulant matrix object
        integer(kind=int32), intent(out) :: info !! success of the routine

        !! Local variables
        type(C_PTR) :: plan
        complex(kind=c_double_complex), allocatable, dimension(:) :: templambda
        real(kind=c_double), allocatable, dimension(:) :: tempc
        integer(kind=int32) :: n,ntransform

        info = structblas_success_
        n = circ%n
        ntransform = n/2+1

        if (allocated(circ%c).and..not.allocated(circ%lambda)) then
            !! we have the first column, but we don't have the eigenvalues
            allocate(circ%lambda(n),templambda(ntransform),tempc(n),stat=info)
            plan = fftw_plan_dft_r2c_1d(n,tempc,templambda,FFTW_ESTIMATE)
            tempc(1:n) = circ%c(1:n) 
            call fftw_execute_dft_r2c(plan,tempc,templambda)
            circ%lambda(1:ntransform) = templambda!/sqrt(real(n))
            circ%lambda(ntransform+1:n) = conjg(templambda(ntransform:2:-1))/sqrt(real(n))
            call fftw_destroy_plan(plan)
            deallocate(templambda,tempc,stat=info)
            circ%isbuild = .true.
        else if (allocated(circ%lambda).and..not.allocated(circ%c)) then
            !! we have the eigenvalues, but we don't have the first column
            allocate(circ%c(n),templambda(ntransform),tempc(n),stat=info)
            plan = fftw_plan_dft_c2r_1d(n,templambda,tempc,FFTW_ESTIMATE)
            templambda(1:ntransform) = circ%lambda(1:ntransform)
            call fftw_execute_dft_c2r(plan,templambda,tempc)
            circ%c(1:n) = tempc(1:n)/real(n)
            call fftw_destroy_plan(plan)
            deallocate(templambda,tempc,stat=info)
            circ%isbuild = .true.
        end if

    end subroutine sb_d_circulant_build

    subroutine sb_d_circulant_tomat(circ,mat,info)
        !! Assembles the circulant matrix in the "full" representation
        use iso_fortran_env
        use structblas_const_mod
        implicit none
        class(dcirculant), intent(in) :: circ !! circulant matrix object
        real(kind=real64), allocatable, dimension(:,:), intent(inout) :: mat !! matrix to be filled
        integer(kind=int32), intent(out) :: info !! success of the routine

        !! Local variables
        integer(kind=int32) :: i,j

        info = structblas_success_

        if (.not.allocated(circ%c)) then
            info = structblas_fail_
            return
        end if

        allocate(mat(circ%n,circ%n),stat=info)

        do i=0,circ%n-1
            do j=0,circ%n-1
                mat(i+1,j+1) = circ%c(modulo(i-j,circ%n)+1)
            end do
        end do

    end subroutine sb_d_circulant_tomat

    subroutine sb_d_circulant_free(circ,info)
        !! Free memory occupated by circulant matrix
        use iso_fortran_env
        use structblas_const_mod
        implicit none
        class(dcirculant), intent(inout) :: circ !! circulant matrix object
        integer(kind=int32), intent(out) :: info !! success of the routine

        info = structblas_success_

        if (allocated(circ%c)) deallocate(circ%c,stat=info)
        if (allocated(circ%lambda)) deallocate(circ%lambda,stat=info)
        circ%n = -1
        circ%isbuild = .false.

    end subroutine

end submodule dcirculant_mod