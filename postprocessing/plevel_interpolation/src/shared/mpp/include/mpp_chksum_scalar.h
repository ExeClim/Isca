    function MPP_CHKSUM_( var, pelist )
!mold is a dummy array to be used by TRANSFER()
!must be same TYPE as result
!result is LONG_KIND, which will actually be int ifdef no_8byte_integers
      integer(LONG_KIND) :: MPP_CHKSUM_
      MPP_TYPE_, intent(in) :: var
      integer, intent(in), optional :: pelist(:)
      integer(LONG_KIND) :: mold(1)
      pointer( p, mold )

      p = LOC(var)
      MPP_CHKSUM_ = mpp_chksum( mold, pelist )
      return
    end function MPP_CHKSUM_
