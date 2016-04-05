    subroutine MPP_WRITE_( unit, field, data, tstamp)
      use mpp_parameter_mod, only : NULLUNIT
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      MPP_TYPE_, intent(in) :: data MPP_RANK_
      real(DOUBLE_KIND), intent(in), optional :: tstamp

      if (unit == NULLUNIT) return
      MPP_WRITE_RECORD_
      return
    end subroutine MPP_WRITE_
