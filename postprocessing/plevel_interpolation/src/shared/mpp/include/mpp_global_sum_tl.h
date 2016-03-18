  subroutine MPP_GLOBAL_SUM_TL_( domain, field, field_tl, gsum, gsum_tl, flags, position, tile_count )
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(inout) :: field(:,: MPP_EXTRA_INDICES_ )
    MPP_TYPE_, intent(inout) :: field_tl(:,: MPP_EXTRA_INDICES_ )
    MPP_TYPE_, intent(inout) :: gsum
    MPP_TYPE_, intent(inout) :: gsum_tl
    integer, intent(in), optional :: position
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: tile_count

    gsum = mpp_global_sum(domain, field, flags, position, tile_count )
    gsum_tl = mpp_global_sum(domain, field_tl, flags, position, tile_count )

    return
  end subroutine MPP_GLOBAL_SUM_TL_
