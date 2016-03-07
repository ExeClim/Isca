subroutine _SUBNAME_(errortype, errormsg1, scalar1, errormsg2, scalar2, errormsg3)
  integer,            intent(in) :: errortype
  _ARRAY1TYPE_, intent(in) :: scalar1
  _ARRAY2TYPE_, intent(in) :: scalar2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar1/), errormsg2, (/scalar2/), errormsg3)

end subroutine _SUBNAME_
