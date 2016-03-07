module platform_mod
!platform-dependent settings
#include <fms_platform.h>
  public
  integer, parameter :: r8_kind=DOUBLE_KIND, r4_kind=FLOAT_KIND, &
                        c8_kind=DOUBLE_KIND, c4_kind=FLOAT_KIND, &
                        l8_kind=LONG_KIND, l4_kind=INT_KIND, &
                        i8_kind=LONG_KIND, i4_kind=INT_KIND, i2_kind=SHORT_KIND
!could additionally define things like OS, compiler...: useful?
end module platform_mod
