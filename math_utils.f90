
module math_utils
  implicit none
contains
  subroutine check_allocation(error, N)
    integer, intent(in) :: error, N
    if (error .ne. 0) then
      print *, "error: couldn't allocate memory, N = ", N
      stop
    end if
  end subroutine check_allocation
end module math_utils
