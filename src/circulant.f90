program circulant
    !! This test program verifies all the routine and functionalities for the circulant matrices. It uses the
    !! Fortuno library for performing unit tests
    use structblas_mod
    use structblas_const_mod
    use iso_fortran_env
    use iso_c_binding
    use fftw_mod
    use fortuno_serial, only : execute_serial_cmd_app, is_equal, test => serial_case_item,&
        & check => serial_check
    implicit none

    call execute_serial_cmd_app(&
    testitems=[&
        test("setc_1", test_setc_1),&
        test("setc_2", test_setc_2),&
        test("setc_3", test_setc_3),&
        test("setlambda_1", test_setlambda_1),&
        test("build_1", test_build_1),&
        test("build_2", test_build_2),&
        test("tomat", test_tomat),&
        test("dgemv",test_dgemv)&
    ]&
    )

contains

        subroutine test_setc_1()
            implicit none
            type(dcirculant) :: c1
            real(kind=real64), dimension(:), allocatable :: cvec
            integer(kind=int32) :: info
            

            allocate(cvec(10),stat=info)
            call check(is_equal(info,0), msg="allocate failed")

            cvec = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            call c1%setc(cvec,info)
            call check(is_equal(info, structblas_success_), msg="setc failed")

            deallocate(cvec,stat=info)
            call check(is_equal(info,0), msg="deallocate failed")

            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
        end subroutine test_setc_1

        subroutine test_setc_2()
            implicit none
            type(dcirculant) :: c1
            real(kind=real64), dimension(:), allocatable :: cvec
            integer(kind=int32) :: info,n
            
            n = 5
            allocate(cvec(10),stat=info)
            call check(is_equal(info,0), msg="allocate failed")

            cvec = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            call c1%setc(cvec,info,n=n)
            call check(is_equal(info, structblas_success_), msg="setc with n failed")

            call check(is_equal(size(c1%c),n), msg="setc%c wrong size")

            deallocate(cvec,stat=info)
            call check(is_equal(info,0), msg="deallocate failed")
            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
        end subroutine test_setc_2

        subroutine test_setc_3()
            implicit none
            type(dcirculant) :: c1
            real(kind=real64), dimension(:), allocatable :: cvec
            integer(kind=int32) :: info,n
            
            n = 15
            allocate(cvec(10),stat=info)
            call check(is_equal(info,0), msg="allocate failed")

            cvec = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            call c1%setc(cvec,info,n=n)
            call check(is_equal(info, structblas_success_), msg="setc with n failed")

            call check(is_equal(size(c1%c),n), msg="setc%c wrong size")

            deallocate(cvec,stat=info)
            call check(is_equal(info,0), msg="deallocate failed")
            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
        end subroutine test_setc_3

        subroutine test_setlambda_1()
            implicit none
            type(dcirculant) :: c1
            complex(kind=real64), dimension(:), allocatable :: lambda
            integer(kind=int32) :: info,i
            
            allocate(lambda(10),stat=info)
            call check(is_equal(info,0), msg="allocate failed")

            do i=1,10
                lambda(i) = exp(zunit*i*dpi/10)
            end do

            call c1%setlambda(lambda,info)
            call check(is_equal(info, structblas_success_), msg="setc with n failed")

            call check(is_equal(size(c1%lambda),10), msg="setc%lambda wrong size")

            deallocate(lambda,stat=info)
            call check(is_equal(info,0), msg="deallocate failed")
            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
        end subroutine test_setlambda_1

        subroutine test_build_1()
            implicit none
            type(dcirculant) :: c1
            real(kind=real64), dimension(:), allocatable :: cvec
            integer(kind=int32) :: info,n
            
            n = 10
            allocate(cvec(2),stat=info)
            call check(is_equal(info,0), msg="allocate failed")

            cvec = [1, -1]
            call c1%setc(cvec,info,n=n)
            call check(is_equal(info, structblas_success_), msg="setc with n failed")
            call check(is_equal(size(c1%c),n), msg="setc%c wrong size")

            call c1%build(info)
            call check(is_equal(info, structblas_success_), msg="build from c failed")
            call check(is_equal(size(c1%lambda),n), msg="setc%lambda wrong size")
            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
        end subroutine test_build_1

        subroutine test_build_2()
            implicit none
            type(dcirculant) :: c1
            complex(kind=real64), dimension(:), allocatable :: lambda
            real(kind=real64) :: norm2err
            integer(kind=int32) :: info,n
            real(c_double), allocatable :: data_in(:)
            complex(c_double_complex), allocatable :: data_out(:)
            integer(int32) :: i
            type(c_ptr) :: planf
                  
            n = 10
            allocate(data_in(n),stat=info)
            call check(is_equal(info,0), msg="allocate data_in failed")
            allocate(data_out(n/2+1),stat=info)
            call check(is_equal(info,0), msg="allocate data_out failed")
            allocate(lambda(n),stat=info)
            call check(is_equal(info,0), msg="allocate lambda failed")

            data_in(1:2) = [1.0d0, -1.0d0]
            planf = fftw_plan_dft_r2c_1d(size(data_in), data_in, data_out, FFTW_ESTIMATE)
            call fftw_execute_dft_r2c(planf, data_in, data_out)
            lambda(1:n/2+1) = data_out(1:n/2+1)
            lambda(n/2+2:n) = data_out(n/2+1:-1:2)
            
            call c1%setlambda(lambda,info)
            call check(is_equal(info, structblas_success_), msg="setlambda failed")
            call check(is_equal(size(c1%lambda),n), msg="setlambda%c wrong size")

            call c1%build(info)
            call check(is_equal(info, structblas_success_), msg="build from lambda failed")
            call check(is_equal(size(c1%c),n), msg="setlambda%c wrong size")
            
            norm2err = dzero
            do i = 1,n
                norm2err = norm2err + (c1%c(i)-data_in(i))**2
            end do
            norm2err = sqrt(norm2err)
            if (norm2err.le.1d1*epsilon(data_in(1))) then
                info = structblas_success_
            else
                info = structblas_fail_
            end if
            call check(is_equal(info,structblas_success_), msg="wrong vector reconstruction")
            
            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
            deallocate(data_in,data_out,lambda, stat=info)
            call check(is_equal(info,0), msg="auxiliary vector free failed")
        end subroutine test_build_2

        subroutine test_tomat()
            implicit none
            type(dcirculant) :: c1
            real(kind=real64), dimension(:), allocatable :: cvec
            real(kind=real64), dimension(:,:), allocatable :: mat,controlmat
            real(kind=real32) :: normfroerr
            integer(kind=int32) :: info,n,i,j
            
            n = 5
            allocate(cvec(2),stat=info)
            call check(is_equal(info,0), msg="cvec allocate failed")
            allocate(controlmat(n,n),stat=info)
            call check(is_equal(info,0), msg="controlmat allocate failed")
            
            do i=1,n
                do j=1,n
                    if (i-j== 1) then
                        controlmat(i,j) = -done
                    else if (i == j) then
                        controlmat(i,j) = done
                    else
                        controlmat(i,j) = dzero
                    end if
                end do
            end do
            controlmat(1,n) = -done


            cvec = [1, -1]
            call c1%setc(cvec,info,n=n)
            call check(is_equal(info, structblas_success_), msg="setc with n failed")
            call check(is_equal(size(c1%c),n), msg="setc%c wrong size")
            
            call c1%build(info)
            call check(is_equal(info, structblas_success_), msg="build from c failed")
            call check(is_equal(size(c1%lambda),n), msg="setc%lambda wrong size")

            call c1%tomat(mat,info)
            call check(is_equal(info, structblas_success_), msg="tomat failed")            

            normfroerr = dzero
            do i=1,n
                do j=1,n
                    normfroerr = normfroerr + (mat(i,j)-controlmat(i,j))**2
                end do
            end do
            normfroerr = sqrt(normfroerr)
            if (normfroerr.le.1d1*epsilon(cvec(1))) then
                info = structblas_success_
            else
                info = structblas_fail_
                print *,"mat = ",mat
                print *,"controlmat = ",controlmat
            end if
            call check(is_equal(info,structblas_success_), msg="wrong matrix reconstruction")


            deallocate(cvec,stat=info)
            call check(is_equal(info,0), msg="deallocate failed")

            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")
        end subroutine test_tomat

        subroutine test_dgemv()
            implicit none
            type(dcirculant) :: c1
            real(kind=real64) :: alpha, beta
            real(kind=real64), dimension(:), allocatable :: cvec,x,y,z
            real(kind=real64), dimension(:,:), allocatable :: mat,controlmat
            real(kind=real32) :: normfroerr
            integer(kind=int32) :: info,n,i,j
            
            n = 20
            allocate(cvec(2),stat=info)
            call check(is_equal(info,0), msg="cvec allocate failed")
            allocate(controlmat(n,n),stat=info)
            call check(is_equal(info,0), msg="controlmat allocate failed")
            allocate(x(n),y(n),z(n),stat=info)
            call check(is_equal(info,0), msg="x/y/z allocate failed")

            call random_number(x)
            call random_number(y)
            !call random_number(alpha)
            !call random_number(beta)
            alpha = 1.0
            beta = 1.0
             
            do i=1,n
                do j=1,n
                    if (i-j== 1) then
                        controlmat(i,j) = -done
                    else if (i == j) then
                        controlmat(i,j) = done
                    else
                        controlmat(i,j) = dzero
                    end if
                end do
            end do
            controlmat(1,n) = -done

            cvec = [1, -1]
            call c1%setc(cvec,info,n=n)
            call check(is_equal(info, structblas_success_), msg="setc with n failed")
            call check(is_equal(size(c1%c),n), msg="setc%c wrong size")
            
            call c1%build(info)
            call check(is_equal(info, structblas_success_), msg="build from c failed")
            call check(is_equal(size(c1%lambda),n), msg="setc%lambda wrong size")

            call c1%tomat(mat,info)
            call check(is_equal(info, structblas_success_), msg="tomat failed")            

            normfroerr = dzero
            do i=1,n
                do j=1,n
                    normfroerr = normfroerr + (mat(i,j)-controlmat(i,j))**2
                end do
            end do
            normfroerr = sqrt(normfroerr)
            if (normfroerr.le.1d1*epsilon(cvec(1))) then
                info = structblas_success_
            else
                info = structblas_fail_
                write(output_unit,'("(Matrix) Error is ",E16.5)') normfroerr
            end if
            call check(is_equal(info,structblas_success_), msg="wrong matrix reconstruction")

            ! Compute control product
            z = alpha*matmul(mat,x) + beta*y
            ! Compute real product
            call sb_gemv('N', alpha, c1, x, beta, y)

            ! Check the error;
            normfroerr = dzero
            do i=1,n
                normfroerr = normfroerr + (y(i) - z(i))**2                
            end do
            normfroerr = sqrt(normfroerr)
            if (normfroerr.le.1d1*epsilon(cvec(1))) then
                info = structblas_success_
            else
                info = structblas_fail_
                write(output_unit,'("(Vector) Error is ",E16.5)') normfroerr
            end if
            call check(is_equal(info,structblas_success_), msg="wrong vector reconstruction")

            deallocate(cvec,controlmat,x,y,z,stat=info)
            call check(is_equal(info,0), msg="deallocate failed")

            call c1%free(info)
            call check(is_equal(info,0), msg="circulant free failed")

        end subroutine test_dgemv

    

end program circulant
