    function MPP_CHKSUM_INT_( var, pelist )
      MPP_TYPE_ :: MPP_CHKSUM_INT_
      MPP_TYPE_, intent(in) :: var MPP_RANK_
      integer, optional :: pelist(:)
      MPP_CHKSUM_INT_ = sum(var)
      call mpp_sum( MPP_CHKSUM_INT_, pelist )
      return
    end function MPP_CHKSUM_INT_
