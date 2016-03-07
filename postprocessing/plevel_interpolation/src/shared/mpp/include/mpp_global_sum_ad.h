  subroutine MPP_GLOBAL_SUM_AD_( domain, field, gsum, position )
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(inout) :: field(:,: MPP_EXTRA_INDICES_ )
    MPP_TYPE_, intent(inout) :: gsum
    integer, intent(in), optional :: position
    
    integer :: i,j, ioff,joff
    type(domain2d), pointer :: Dom => NULL()

    Dom => get_domain(domain, position)

    if( size(field,1).EQ.Dom%x(1)%compute%size .AND. size(field,2).EQ.Dom%y(1)%compute%size )then
!field is on compute domain
        ioff = -Dom%x(1)%compute%begin + 1
        joff = -Dom%y(1)%compute%begin + 1
    else if( size(field,1).EQ.Dom%x(1)%data%size .AND. size(field,2).EQ.Dom%y(1)%data%size )then
!field is on data domain
        ioff = -Dom%x(1)%data%begin + 1
        joff = -Dom%y(1)%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_SUM_AD_: incoming field array must match either compute domain or data domain.' )
    end if

        do j = Dom%y(1)%compute%begin, Dom%y(1)%compute%end
           do i = Dom%x(1)%compute%begin, Dom%x(1)%compute%end
              field(i+ioff:i+ioff,j+joff:j+joff MPP_EXTRA_INDICES_)= gsum
           end do
        end do
      
        gsum = 0.
    return
  end subroutine MPP_GLOBAL_SUM_AD_
