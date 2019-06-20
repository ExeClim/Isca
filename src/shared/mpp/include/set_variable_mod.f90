module set_variable_mod
  implicit none
  interface func
     module procedure func_r4, func_l8
  end interface

  
  public    :: func
  private   :: func_r4, func_l8

contains

  function func_r4(val) result(d)
    real, intent(in) :: val
    real             :: d
    d = 0
  end function func_r4

  function func_l8(val) result(d)
    logical, intent(in) :: val
    logical             :: d
    d = .F.
  end function func_l8
    
end module set_variable_mod
