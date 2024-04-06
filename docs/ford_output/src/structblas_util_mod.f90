module structblas_util_mod

    private
    public :: sb_toupper
    
    contains

pure Function sb_toupper (str) result(string)
!! Changes a string to upper case
    implicit None
    
    character(*), Intent(In) :: str
    character(LEN(str))      :: string

    integer :: ic, i

    character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

end function sb_toupper


end module structblas_util_mod