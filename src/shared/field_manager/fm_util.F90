!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fm_util_mod  !{
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Utility routines for the field manager
!</OVERVIEW>
!
!<DESCRIPTION>
! This module provides utility routines for the field manager.
! Basically, it provides for error catching, reporting and
! termination while interfacing with the field manager.
!</DESCRIPTION>
!
! <INFO>
! </INFO>
!

use field_manager_mod, only: fm_string_len, fm_path_name_len, fm_field_name_len, fm_type_name_len
use field_manager_mod, only: fm_get_type, fm_get_index, fm_get_length
use field_manager_mod, only: fm_get_current_list, fm_new_list, fm_change_list, fm_loop_over_list
use field_manager_mod, only: fm_new_value, fm_get_value
use field_manager_mod, only: fm_exists, fm_dump_list
use fms_mod,           only: FATAL, stdout
use mpp_mod,           only: mpp_error

implicit none

private

public  fm_util_start_namelist
public  fm_util_end_namelist
public  fm_util_check_for_bad_fields
public  fm_util_set_caller
public  fm_util_reset_caller
public  fm_util_set_no_overwrite
public  fm_util_reset_no_overwrite
public  fm_util_set_good_name_list
public  fm_util_reset_good_name_list
public  fm_util_get_length
public  fm_util_get_integer
public  fm_util_get_logical
public  fm_util_get_real
public  fm_util_get_string
public  fm_util_get_integer_array
public  fm_util_get_logical_array
public  fm_util_get_real_array
public  fm_util_get_string_array
public  fm_util_set_value
public  fm_util_set_value_integer_array
public  fm_util_set_value_logical_array
public  fm_util_set_value_real_array
public  fm_util_set_value_string_array
public  fm_util_set_value_integer
public  fm_util_set_value_logical
public  fm_util_set_value_real
public  fm_util_set_value_string
!public  fm_util_get_index
public  fm_util_get_index_list
public  fm_util_get_index_string

!
!       Public variables
!

character(len=128), public      :: fm_util_default_caller = ' '

!
!       private parameters
!

character(len=48), parameter    :: mod_name = 'fm_util_mod'

!
!       Private variables
!

character(len=128)              :: save_default_caller = ' '
character(len=128)              :: default_good_name_list = ' '
character(len=128)              :: save_default_good_name_list = ' '
logical                         :: default_no_overwrite = .false.
logical                         :: save_default_no_overwrite = .false.
character(len=fm_path_name_len) :: save_current_list
character(len=fm_path_name_len) :: save_path
character(len=fm_path_name_len) :: save_name
character(len=128) :: version = '$Id: fm_util.F90,v 17.0 2009/07/21 03:19:16 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'

!
!        Interface definitions for overloaded routines
!

!interface  fm_util_get_value  !{
  !module procedure  fm_util_get_value_integer
  !module procedure  fm_util_get_value_logical
  !module procedure  fm_util_get_value_real
  !module procedure  fm_util_get_value_string
  !module procedure  fm_util_get_value_integer_array
  !module procedure  fm_util_get_value_logical_array
  !module procedure  fm_util_get_value_real_array
  !module procedure  fm_util_get_value_string_array
!end interface  !}

interface  fm_util_set_value  !{
  module procedure  fm_util_set_value_integer_array
  module procedure  fm_util_set_value_logical_array
  module procedure  fm_util_set_value_real_array
  module procedure  fm_util_set_value_string_array
  module procedure  fm_util_set_value_integer
  module procedure  fm_util_set_value_logical
  module procedure  fm_util_set_value_real
  module procedure  fm_util_set_value_string
end interface  !}

!interface  fm_util_get_index  !{
  !module procedure  fm_util_get_index_list
  !module procedure  fm_util_get_index_string
!end interface  !}


contains


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_caller">
!
! <DESCRIPTION>
! Set the default value for the optional "caller" variable used in many of these
! subroutines. If the argument is blank, then set the default to blank, otherwise
! the deault will have brackets placed around the argument.
!
! </DESCRIPTION>
!

subroutine fm_util_set_caller(caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)          :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_set_caller'

!
!       Local variables
!

!
!       save the default caller string
!

save_default_caller = fm_util_default_caller

!
!       set the default caller string
!

if (caller .eq. ' ') then  !{
  fm_util_default_caller = ' '
else  !}{
  fm_util_default_caller = '[' // trim(caller) // ']'
endif  !}

return

end subroutine fm_util_set_caller  !}
! </SUBROUTINE> NAME="fm_util_set_caller"


!#######################################################################
! <SUBROUTINE NAME="fm_util_reset_caller">
!
! <DESCRIPTION>
! Reset the default value for the optional "caller" variable used in many of these
! subroutines to blank.
!
! </DESCRIPTION>
!

subroutine fm_util_reset_caller  !{

implicit none

!
!       arguments
!

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_reset_caller'

!
!       Local variables
!

!
!       reset the default caller string
!

fm_util_default_caller = save_default_caller
save_default_caller = ' '

return

end subroutine fm_util_reset_caller  !}
! </SUBROUTINE> NAME="fm_util_reset_caller"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_good_name_list">
!
! <DESCRIPTION>
! Set the default value for the optional "good_name_list" variable used in many of these
! subroutines.
!
! </DESCRIPTION>
!

subroutine fm_util_set_good_name_list(good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)          :: good_name_list

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_set_good_name_list'

!
!       Local variables
!

!
!       save the default good_name_list string
!

save_default_good_name_list = default_good_name_list

!
!       set the default good_name_list string
!

default_good_name_list = good_name_list

return

end subroutine fm_util_set_good_name_list  !}
! </SUBROUTINE> NAME="fm_util_set_good_name_list"


!#######################################################################
! <SUBROUTINE NAME="fm_util_reset_good_name_list">
!
! <DESCRIPTION>
! Reset the default value for the optional "good_name_list" variable used in many of these
! subroutines to the saved value.
!
! </DESCRIPTION>
!

subroutine fm_util_reset_good_name_list  !{

implicit none

!
!       arguments
!

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_reset_good_name_list'

!
!       Local variables
!

!
!       reset the default good_name_list string
!

default_good_name_list = save_default_good_name_list
save_default_good_name_list = ' '

return

end subroutine fm_util_reset_good_name_list  !}
! </SUBROUTINE> NAME="fm_util_reset_good_name_list"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_no_overwrite">
!
! <DESCRIPTION>
! Set the default value for the optional "no_overwrite" variable used in some of these
! subroutines.
!
! </DESCRIPTION>
!

subroutine fm_util_set_no_overwrite(no_overwrite)  !{

implicit none

!
!       arguments
!

logical, intent(in)          :: no_overwrite

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_set_no_overwrite'

!
!       Local variables
!

!
!       save the default no_overwrite string
!

save_default_no_overwrite = default_no_overwrite

!
!       set the default no_overwrite value
!

default_no_overwrite = no_overwrite

return

end subroutine fm_util_set_no_overwrite  !}
! </SUBROUTINE> NAME="fm_util_set_no_overwrite"


!#######################################################################
! <SUBROUTINE NAME="fm_util_reset_no_overwrite">
!
! <DESCRIPTION>
! Reset the default value for the optional "no_overwrite" variable used in some of these
! subroutines to false.
!
! </DESCRIPTION>
!

subroutine fm_util_reset_no_overwrite  !{

implicit none

!
!       arguments
!

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_reset_no_overwrite'

!
!       Local variables
!

!
!       reset the default no_overwrite value
!

default_no_overwrite = save_default_no_overwrite
save_default_no_overwrite = .false.

return

end subroutine fm_util_reset_no_overwrite  !}
! </SUBROUTINE> NAME="fm_util_reset_no_overwrite"


!#######################################################################
! <SUBROUTINE NAME="fm_util_check_for_bad_fields">
!
! <DESCRIPTION>
! Check for unrecognized fields in a list
!
! </DESCRIPTION>
!

subroutine fm_util_check_for_bad_fields(list, good_fields, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                    :: list
character(len=*), intent(in), dimension(:)      :: good_fields
character(len=*), intent(in), optional          :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_check_for_bad_fields'

!
!       Local variables
!

logical                                 :: fm_success
integer                                 :: i
integer                                 :: ind
integer                                 :: list_length
integer                                 :: good_length
character(len=fm_type_name_len)         :: typ
character(len=fm_field_name_len)        :: name
logical                                 :: found
character(len=256)                      :: error_header
character(len=256)                      :: warn_header
character(len=256)                      :: note_header
character(len=128)                      :: caller_str
integer                         :: out_unit

out_unit = stdout()

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a list is given (fatal if not)
!

if (list .eq. ' ') then  !{
  write (out_unit,*) trim(error_header) // ' Empty list given'
  call mpp_error(FATAL, trim(error_header) // ' Empty list given')
endif  !}

!
!       Check that we have been given a list
!

if (fm_get_type(list) .ne. 'list') then  !{
  write (out_unit,*) trim(error_header) // ' Not given a list: ' // trim(list)
  call mpp_error(FATAL, trim(error_header) // ' Not given a list: ' // trim(list))
endif  !}

!
!       Get the list length
!

list_length = fm_get_length(list)
if (list_length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(list))
endif  !}

!
!       Get the number of good fields
!

good_length = size(good_fields)

if (list_length .lt. good_length) then  !{

!
!       If the list length is less than the number of good fields this is an error 
!       as the list should be fully populated and we'll check which extra fields
!       are given in good_fields
!

  write (out_unit,*) trim(error_header), ' List length < number of good fields (',       &
       list_length, ' < ', good_length, ') in list ', trim(list)

  write (out_unit,*)
  write (out_unit,*) 'The list contains the following fields:'
  fm_success= fm_dump_list(list, .false.)
  write (out_unit,*)
  write (out_unit,*) 'The supposed list of good fields is:'
  do i = 1, good_length  !{
    if (fm_exists(trim(list) // '/' // good_fields(i))) then  !{
      write (out_unit,*) 'List field: "', trim(good_fields(i)), '"'
    else  !}{
      write (out_unit,*) 'EXTRA good field: "', trim(good_fields(i)), '"'
    endif  !}
  enddo  !} i
  write (out_unit,*)

  call mpp_error(FATAL, trim(error_header) //                                           &
       ' List length < number of good fields for list: ' // trim(list))

elseif (list_length .gt. good_length) then  !}{

!
!       If the list length is greater than the number of good fields this is an error
!       as the there should not be any more fields than those given in the good fields list
!       and we'll check which extra fields are given in the list
!

  write (out_unit,*) trim(warn_header), 'List length > number of good fields (',        &
       list_length, ' > ', good_length, ') in list ', trim(list)

  write (out_unit,*) trim(error_header), ' Start of list of fields'
  do while (fm_loop_over_list(list, name, typ, ind))  !{
    found = .false.
    do i = 1, good_length  !{
      found = found .or. (name .eq. good_fields(i))
    enddo  !} i
    if (found) then  !{
      write (out_unit,*) 'Good list field: "', trim(name), '"'
    else  !}{
      write (out_unit,*) 'EXTRA list field: "', trim(name), '"'
    endif  !}
  enddo  !}
  write (out_unit,*) trim(error_header), ' End of list of fields'

  call mpp_error(FATAL, trim(error_header) //                                           &
       ' List length > number of good fields for list: ' // trim(list))

endif  !}

!
!       If the list length equals the number of good fields then all is good
!

return

end subroutine fm_util_check_for_bad_fields  !}
! </SUBROUTINE> NAME="fm_util_check_for_bad_fields"


!#######################################################################
! <FUNCTION NAME="fm_util_get_length">
!
! <DESCRIPTION>
! Get the length of an element of the Field Manager tree
! </DESCRIPTION>
!
function fm_util_get_length(name, caller)       &
         result (field_length)  !{

implicit none

!
!       Return type
!

integer :: field_length

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_length'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Get the field's length
!

field_length = fm_get_length(name)
if (field_length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
endif  !}

return

end function fm_util_get_length  !}
! </FUNCTION> NAME="fm_util_get_length"


!#######################################################################
! <FUNCTION NAME="fm_util_get_index_string">
!
! <DESCRIPTION>
! Get the index of an element of a string in the Field Manager tree
! </DESCRIPTION>
!
function fm_util_get_index_string(name, string, caller)       &
         result (fm_index)  !{

implicit none

!
!       Return type
!

integer :: fm_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: string
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_index_string'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
character(len=fm_string_len)    :: fm_string
integer                         :: i
integer                         :: length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check the field's type and get the index
!

fm_index = 0
fm_type = fm_get_type(name)
if (fm_type .eq. 'string') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    do i = 1, length  !{
      if (.not. fm_get_value(name, fm_string, index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
      if (fm_string .eq. string) then  !{
        fm_index = i
        exit
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

!if (fm_index .eq. 0) then  !{
  !call mpp_error(FATAL, trim(error_header) // ' "' // trim(string) // '" does not exist in ' // trim(name))
!endif  !}

return

end function fm_util_get_index_string  !}
! </FUNCTION> NAME="fm_util_get_index_string"


!#######################################################################
! <FUNCTION NAME="fm_util_get_index_list">
!
! <DESCRIPTION>
! Get the length of an element of the Field Manager tree
! </DESCRIPTION>
!
function fm_util_get_index_list(name, caller)       &
         result (fm_index)  !{

implicit none

!
!       Return type
!

integer :: fm_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_index_list'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=fm_type_name_len) :: fm_type

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check the field's type and get the index
!

fm_index = 0
fm_type = fm_get_type(name)
if (fm_type .eq. 'list') then  !{
  fm_index = fm_get_index(name)
  if (fm_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' List does not exist: ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' List does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}


return

end function fm_util_get_index_list  !}
! </FUNCTION> NAME="fm_util_get_index_list"


!#######################################################################
! <FUNCTION NAME="fm_util_get_integer_array">
!
! <DESCRIPTION>
! Get an integer value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_integer_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

integer, pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_integer_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'integer') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_integer_array  !}
! </FUNCTION> NAME="fm_util_get_integer_array"


!#######################################################################
! <FUNCTION NAME="fm_util_get_logical_array">
!
! <DESCRIPTION>
! Get a logical value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_logical_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

logical, pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_logical_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'logical') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_logical_array  !}
! </FUNCTION> NAME="fm_util_get_logical_array"


!#######################################################################
! <FUNCTION NAME="fm_util_get_real_array">
!
! <DESCRIPTION>
! Get a real value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_real_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

real, pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_real_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'real') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_real_array  !}
! </FUNCTION> NAME="fm_util_get_real_array"


!#######################################################################
! <FUNCTION NAME="fm_util_get_string_array">
!
! <DESCRIPTION>
! Get a string value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_string_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

character(len=fm_string_len), pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_string_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'string') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_string_array  !}
! </FUNCTION> NAME="fm_util_get_string_array"


!#######################################################################
! <FUNCTION NAME="fm_util_get_integer">
!
! <DESCRIPTION>
! Get an integer value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_integer(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

integer :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
integer, intent(in), optional           :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_integer'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'integer') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Field does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_integer  !}
! </FUNCTION> NAME="fm_util_get_integer"


!#######################################################################
! <FUNCTION NAME="fm_util_get_logical">
!
! <DESCRIPTION>
! Get a logical value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_logical(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

logical :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_logical'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'logical') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Field does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_logical  !}
! </FUNCTION> NAME="fm_util_get_logical"


!#######################################################################
! <FUNCTION NAME="fm_util_get_real">
!
! <DESCRIPTION>
! Get a real value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_real(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

real :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
real, intent(in), optional              :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_real'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'real') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Field does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_real  !}
! </FUNCTION> NAME="fm_util_get_real"


!#######################################################################
! <FUNCTION NAME="fm_util_get_string">
!
! <DESCRIPTION>
! Get a string value from the Field Manager tree.
! </DESCRIPTION>
!
function fm_util_get_string(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

character(len=fm_string_len) :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
character(len=*), intent(in), optional  :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_get_string'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'string') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Field does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function fm_util_get_string  !}
! </FUNCTION> NAME="fm_util_get_string"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_integer_array">
!
! <DESCRIPTION>
! Set an integer array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_integer_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
integer, intent(in)                                     :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_integer_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, 0, index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_integer_array  !}
! </SUBROUTINE> NAME="fm_util_set_value_integer_array"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_logical_array">
!
! <DESCRIPTION>
! Set a logical array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_logical_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
logical, intent(in)                                     :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_logical_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, .false., index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_logical_array  !}
! </SUBROUTINE> NAME="fm_util_set_value_logical_array"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_real_array">
!
! <DESCRIPTION>
! Set a real array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_real_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
real, intent(in)                                        :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_real_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, 0.0, index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_real_array  !}
! </SUBROUTINE> NAME="fm_util_set_value_real_array"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_string_array">
!
! <DESCRIPTION>
! Set a string array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_string_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
character(len=*), intent(in)                            :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_string_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, ' ', index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_string_array  !}
! </SUBROUTINE> NAME="fm_util_set_value_string_array"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_integer">
!
! <DESCRIPTION>
! Set an integer value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_integer(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
integer, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_integer'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_integer  !}
! </SUBROUTINE> NAME="fm_util_set_value_integer"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_logical">
!
! <DESCRIPTION>
! Set a logical value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_logical(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
logical, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_logical'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_logical  !}
! </SUBROUTINE> NAME="fm_util_set_value_logical"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_real">
!
! <DESCRIPTION>
! Set a real value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_real(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
real, intent(in)                        :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_real'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_real  !}
! </SUBROUTINE> NAME="fm_util_set_value_real"


!#######################################################################
! <SUBROUTINE NAME="fm_util_set_value_string">
!
! <DESCRIPTION>
! Set a string value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine fm_util_set_value_string(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'fm_util_set_value_string'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_exists(good_name_list_use)) then  !{
    add_name = fm_util_get_index_string(good_name_list_use, name,               &
       caller = caller_str) .le. 0              ! true if name does not exist in string array
  else  !}{
    add_name = .true.                           ! always add to new list
  endif  !}
  if (add_name .and. fm_exists(name)) then  !{
    if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) //                               &
           ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
    endif  !}
  endif  !}
endif  !}

return

end subroutine fm_util_set_value_string  !}
! </SUBROUTINE> NAME="fm_util_set_value_string"


!#######################################################################
! <SUBROUTINE NAME="fm_util_start_namelist">
!
! <DESCRIPTION>
! Start processing a namelist
! </DESCRIPTION>
!
subroutine fm_util_start_namelist(path, name, caller, no_overwrite, check)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: path
character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
logical,          intent(in), optional  :: no_overwrite
logical,          intent(in), optional  :: check

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_start_namelist'

!
!       Local variables
!

integer                         :: namelist_index
character(len=fm_path_name_len) :: path_name
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: out_unit

out_unit = stdout()

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Concatenate the path and name
!

if (path .eq. ' ') then  !{
  path_name = name
else  !}{
  path_name = trim(path) // '/' // name
endif  !}
save_path = path
save_name = name

!
!       set the default caller string, if desired
!

if (present(caller)) then  !{
  call fm_util_set_caller(caller)
else  !}{
  call fm_util_reset_caller
endif  !}

!
!       set the default no_overwrite flag, if desired
!

if (present(no_overwrite)) then  !{
  call fm_util_set_no_overwrite(no_overwrite)
else  !}{
  call fm_util_reset_no_overwrite
endif  !}

!
!       set the default good_name_list string, if desired
!

if (present(check)) then  !{
  if (check) then  !{
    call fm_util_set_good_name_list('/ocean_mod/GOOD/namelists/' // trim(path_name) // '/good_list')
  else  !}{
    call fm_util_reset_good_name_list
  endif  !}
else  !}{
  call fm_util_reset_good_name_list
endif  !}

!
!       Process the namelist
!

write (out_unit,*)
write (out_unit,*) trim(note_header), ' Processing namelist ', trim(path_name)

!
!       Check whether the namelist already exists. If so, then use that one
!

namelist_index = fm_get_index('/ocean_mod/namelists/' // trim(path_name))
if (namelist_index .gt. 0) then  !{

  !write (out_unit,*) trim(note_header), ' Namelist already set with index ', namelist_index

else  !}{

!
!       Set a new namelist and get its index
!

  namelist_index = fm_new_list('/ocean_mod/namelists/' // trim(path_name), create = .true.)
  if (namelist_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set namelist ' // trim(path_name))
  endif  !}

endif  !}

!
!       Add the namelist name to the list of good namelists, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/namelists/' // trim(path) // '/good_values',    &
                 name, append = .true., create = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(name) // ' to "' // trim(path) // '/good_values" list')
endif  !}

!
!       Change to the new namelist, first saving the current list
!

save_current_list = fm_get_current_list()
if (save_current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list('/ocean_mod/namelists/' // trim(path_name))) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the namelist ' // trim(path_name))
endif  !}

return

end subroutine fm_util_start_namelist  !}
! </SUBROUTINE> NAME="fm_util_start_namelist"


!#######################################################################
! <SUBROUTINE NAME="fm_util_end_namelist">
!
! <DESCRIPTION>
! Finish up processing a namelist
! </DESCRIPTION>
!
subroutine fm_util_end_namelist(path, name, caller, check)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: path
character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
logical,          intent(in), optional  :: check

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'fm_util_end_namelist'

!
!       Local variables
!

character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
character(len=fm_path_name_len)                         :: path_name
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a path is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check that the path ane name match the preceding call to
!       fm_util_start_namelist
!

if (path .ne. save_path) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Path "' // trim(path) // '" does not match saved path "' // trim(save_path) // '"')
elseif (name .ne. save_name) then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Name "' // trim(name) // '" does not match saved name "' // trim(save_name) // '"')
endif  !}

!
!       Concatenate the path and name
!

if (path .eq. ' ') then  !{
  path_name = name
else  !}{
  path_name = trim(path) // '/' // name
endif  !}
save_path = ' '
save_name = ' '

!
!       Check for any errors in the number of fields in this list
!

if (present(check)) then  !{
  if (check) then  !{
    if (caller_str .eq. ' ') then  !{
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
    endif  !}
    good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(path_name) // '/good_list',            &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    if (associated(good_list)) then  !{
      call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(path_name), good_list, caller = caller_str)
      deallocate(good_list)
    else  !}{
      call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(path_name) // '" list')
    endif  !}
  endif  !}
endif  !}

!
!       Change back to the saved list
!

if (save_current_list .ne. ' ') then  !{
  if (.not. fm_change_list(save_current_list)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not change to the saved list: ' // trim(save_current_list))
  endif  !}
endif  !}
save_current_list = ' '

!
!       reset the default caller string
!

call fm_util_reset_caller

!
!       reset the default no_overwrite string
!

call fm_util_reset_no_overwrite

!
!       reset the default good_name_list string
!

call fm_util_reset_good_name_list

return

end subroutine fm_util_end_namelist  !}
! </SUBROUTINE> NAME="fm_util_end_namelist"


end module fm_util_mod  !}
