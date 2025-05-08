
module io_utils
  implicit none
contains
  subroutine get_input_args(filename, N)
    character(len=100), intent(out) :: filename
    integer, intent(out) :: N
    character(len=100) :: option1, option2

    call get_command_argument(1, option1)
    call get_command_argument(2, option2)

    filename = option1
    read(option2, '(I9)') N
  end subroutine get_input_args
end module io_utils
