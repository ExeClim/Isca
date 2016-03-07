subroutine _SUBNAME_(errortype, errormsg1, array, errormsg2, scalar, errormsg3)
  integer,            intent(in) :: errortype
  _ARRAY1TYPE_, dimension(:), intent(in) :: array
  _ARRAY2TYPE_,               intent(in) :: scalar
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, array, errormsg2, (/scalar/), errormsg3)

end subroutine _SUBNAME_
