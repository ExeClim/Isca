    function MPP_CHKSUM_( var, pelist )
!mold is a dummy array to be used by TRANSFER()
!must be same TYPE as result
!result is LONG_KIND, which will actually be int ifdef no_8byte_integers
      integer(LONG_KIND) :: MPP_CHKSUM_, mold(1)
      MPP_TYPE_, intent(in) :: var MPP_RANK_
      integer, intent(in), optional :: pelist(:)

      MPP_CHKSUM_ = mpp_chksum( TRANSFER(var,mold), pelist )
      return
    end function MPP_CHKSUM_
