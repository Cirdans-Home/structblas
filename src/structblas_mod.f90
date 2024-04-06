module structblas_mod
  use iso_fortran_env  !! Intrinsic module providing constants, derived types, and intrinsic procedures relating to the Fortran environment.
  use circulant_mod
  implicit none
  
  !! general operation interface
  interface sb_gemv
    module procedure d_circ_gemv
  end interface sb_gemv
  
contains
  subroutine libraryinfo
    use iso_fortran_env, only: output_unit,compiler_options,compiler_version
    implicit none

    write(output_unit,'("Welcome to the STRUCTBLAS library")')
    write(output_unit,'("This library has been compiled with: ",a)')  compiler_version()
    write(output_unit,'("Compiler options are: ",a)') compiler_options()
    
  end subroutine libraryinfo
end module structblas_mod
