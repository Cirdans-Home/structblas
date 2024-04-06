module structblas_const_mod
    !! This module contains the constants  used everywhere in the code.
    use iso_fortran_env 
    implicit none

    !! Integer constants
    integer(kind=int32), parameter :: structblas_success_ = 0
    integer(kind=int32), parameter :: structblas_fail_ = 1

    !! Short-hands for numerical constants
    real(kind=real64), parameter :: dzero = 0.0d0
    real(kind=real64), parameter :: done = 1.0d0
    complex(kind=real64), parameter :: zunit = (0.0d0,1.0d0)
    real(kind=real64), parameter :: dpi = 3.14159265358979323846

end module