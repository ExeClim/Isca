subroutine _SUBNAME_(errortype, errormsg1, scalar2, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  _ARRAY1TYPE_,               intent(in) :: scalar2
  _ARRAY2TYPE_, dimension(:), intent(in) :: array2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar2/), errormsg2, array2, errormsg3)

end subroutine _SUBNAME_
