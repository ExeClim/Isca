subroutine _SUBNAME_(errortype, errormsg1, array1, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  _ARRAY1TYPE_, dimension(:), intent(in) :: array1
  _ARRAY2TYPE_, dimension(:), intent(in) :: array2
  character(len=*),      intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array1)) 
  string = trim(string)//errormsg2//trim(array_to_char(array2)) 
  if(present(errormsg3)) string = trim(string)//errormsg3
  call mpp_error_basic( errortype, trim(string))

end subroutine _SUBNAME_
