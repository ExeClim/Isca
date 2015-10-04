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

module field_manager_mod
#ifndef MAXFIELDS_ 
#define MAXFIELDS_ 150
#endif

#ifndef MAXFIELDMETHODS_
#define MAXFIELDMETHODS_ 150
#endif

!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> William Cooke
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
! </REVIEWER>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matthew Harrison
! </REVIEWER>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
! </REVIEWER>
!
! <HISTORY
!  SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/field_manager/field_manager.F90"/>

! <OVERVIEW>

! The field manager reads entries from a field table and stores this
! information along with the type  of field it belongs to. This allows
! the component models to query the field manager to see if  non-default
! methods of operation are desired. In essence the field table is a
! powerful type of namelist. Default values can be provided for all the
! fields through a namelist, individual fields can be modified  through
! the field table however.

!</OVERVIEW>

! <DESCRIPTION>
!
! An example of field table entries could be
! <PRE>
!"tracer","atmos_mod","sphum"/
!
!"tracer","atmos_mod","sf6"
!"longname","sulf_hex"
!"advection_scheme_horiz","2nd_order"
!"Profile_type","Fixed","surface_value = 0.0E+00"/
!
!"prog_tracers","ocean_mod","age_global"
!horizontal-advection-scheme = mdfl_sweby
!vertical-advection-scheme = mdfl_sweby
!restart_file = ocean_age.res.nc
! </PRE>
! 
! The field table consists of entries in the following format.
!
! The first line of an entry should consist of three quoted strings.
!
! The first quoted string will tell the field manager what type of 
! field it is.
! 
! The second quoted string will tell the field manager which model the 
! field is being applied to.
! The supported types at present are
!<PRE>
!      "coupler_mod" for the coupler,
!      "atmos_mod" for the atmosphere model,
!      "ocean_mod" for the ocean model,
!      "land_mod" for the land model, and,
!      "ice_mod" for the ice model.
!</PRE>
! The third quoted string should be a unique name that can be used as a
! query.
!
! The second and following lines of each entry are called methods in
! this context. Methods can be developed within any module and these
! modules can query the field manager to find any methods that are
! supplied in the field table.
!
! These lines can be coded quite flexibly.
!
! The line can consist of two or three quoted strings or a simple unquoted 
! string.
!
! If the line consists two or three quoted strings, then the first string will 
! be an identifier that the querying module will ask for.
!
! The second string will be a name that the querying module can use to
! set up values for the module. 
!
! The third string, if present, can supply parameters to the calling module that can be
! parsed and used to further modify values.
!
! If the line consists of a simple unquoted string then quotes are not allowed 
! in any part of the line.
!
! An entry is ended with a backslash (/) as the final character in a
! row.
! 
! Comments can be inserted in the field table by having a # as the
! first character in the line.
! 
! In the example above we have three field entries. 
! 
! The first is a simple declaration of a tracer called "sphum". 
!
! The second is for a tracer called "sf6". In this case a field named
! "longname" will be given the value "sulf_hex". A field named 
! "advection_scheme_horiz" will be given the value "2nd_order". Finally a field
! name "Profile_type" will be given a child field called "Fixed", and that field
! will be given a field called "surface_value" with a real value of 0.0E+00.
!
! The third entry is an example of a oceanic age tracer. Note that the 
! method lines are formatted differently here. This is the flexibility mentioned 
! above.
! 
! With these formats, a number of restrictions are required. 
!
! The following formats are equally valid.
!<PRE>
!      "longname","sulf_hex"
!      "longname = sulf_hex"
!      longname = sulf_hex
!</PRE>
! However the following is not valid.
!<PRE>
!      longname = "sulf_hex"
!</PRE>
!
! In the SF6 example above the last line of the entry could be written in the 
! following ways.
!<PRE>
!      "Profile_type","Fixed","surface_value = 0.0E+00"/
!      Profile_type/Fixed/surface_value = 0.0E+00/
!</PRE>
!
! Values supplied with fields are converted to the various types with the
! following assumptions.
!<PRE>
! Real values : These values contain a decimal point or are in exponential format.
!    These values only support e or E format for exponentials.
!    e.g. 10.0, 1e10 and 1E10 are considered to be real numbers.
!
! Integer values : These values only contain numbers. 
!    e.g 10 is an integer. 10.0 and 1e10 are not.
!
! Logical values : These values are supplied as one of the following formats.
!    T, .T., TRUE, .TRUE.
!    t, .t., true, .true.
!    F, .F., FALSE, .FALSE.
!    f, .f., false, .false.
!    These will be converted to T or F in a dump of the field.
!
! Character strings : These values are assumed to be strings if a character 
!    other than an e (or E) is in the value. Numbers can be suppled in the value.
!    If the value does not meet the criteria for a real, integer or logical type,
!    it is assumed to be a character type.
!</PRE>
! The entries within the field table can be designed by the individual
! authors of code to allow modification of their routines.
!
! </DESCRIPTION>

use    mpp_mod, only : mpp_error,   &
                       FATAL,       &
                       NOTE,        &
                       WARNING,     &
                       mpp_pe,      &
                       mpp_root_pe, &
                       stdlog,      &
                       stdout
use mpp_io_mod, only : mpp_io_init, &
                       mpp_open,    &
                       mpp_close,   &
                       MPP_ASCII,   &
                       MPP_RDONLY
use    fms_mod, only : lowercase,   &
                       file_exist,  &
                       write_version_number

implicit none
private


character(len=128) :: version = '$Id: field_manager.F90,v 19.0.6.1 2012/09/13 15:23:54 Seth.Underwood Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'
logical            :: module_is_initialized  = .false.

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Public routines
!        Interface definitions (optional arguments are in [brackets]):
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
public :: field_manager_init   ! (nfields, [table_name]) returns number of fields
public :: field_manager_end    ! ()
public :: find_field_index     ! (model, field_name) or (list_path)
public :: find_field_index_old ! (model, field_name) returns index of field_name in 
public :: find_field_index_new ! (list_path) returns index of field_name in 
                               ! component model model
public :: get_field_info       ! (n,fld_type,fld_name,model,num_methods)
                               ! Returns parameters relating to field n.
public :: get_field_method     ! (n, m, method) Returns the m-th method of field n
public :: get_field_methods    ! (n, methods) Returns the methods related to field n
public :: parse                ! (text, label, values) Overloaded function to parse integer,
                               ! real or character. Parse returns the number of values 
                               ! decoded (> 1 => an array of values)
public :: fm_change_list       ! (list) return success
public :: fm_change_root       ! (list) return success
public :: fm_dump_list         ! (list [, recursive]) return success
public :: fm_exists            ! (field) return success
public :: fm_get_index         ! (field) return index
public :: fm_get_current_list  ! () return path
public :: fm_get_length        ! (list) return length
public :: fm_get_type          ! (field) return string
public :: fm_get_value         ! (entry, value [, index]) return success !! generic
public :: fm_get_value_integer !   as above (overloaded function)
public :: fm_get_value_logical !   as above (overloaded function)
public :: fm_get_value_real    !   as above (overloaded function)
public :: fm_get_value_string  !   as above (overloaded function)
public :: fm_intersection      ! (lists, num_lists) return fm_array_list pointer
public :: fm_loop_over_list    ! (list, name, type, index) return success
public :: fm_new_list          ! (list [, create] [, keep]) return index
public :: fm_new_value         ! (entry, value [, create] [, index]) return index !! generic
public :: fm_new_value_integer !   as above (overloaded function)
public :: fm_new_value_logical !   as above (overloaded function)
public :: fm_new_value_real    !   as above (overloaded function)
public :: fm_new_value_string  !   as above (overloaded function)
public :: fm_reset_loop        ! ()
public :: fm_return_root       ! () return success
public :: fm_modify_name       ! (oldname, newname) return success
public :: fm_query_method      ! (name, method_name, method_control) return success and 
                               ! name and control strings
public :: fm_find_methods      ! (list, methods, control) return success and name and 
                               ! control strings.
public :: fm_copy_list         ! (list, suffix, [create]) return index
public :: fm_set_verbosity     ! ([verbosity])

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!   Private routines
!   Interface definitions (optional arguments are in [brackets]):
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

private :: create_field        ! (list_p, name) return field pointer
private :: dump_list           ! (list_p, recursive, depth) return success
private :: find_base           ! (field, path, base)
private :: find_field          ! (field, list_p) return field pointer
private :: find_head           ! (field, head, rest)
private :: find_list           ! (list, list_p, create) return field pointer
private :: get_field           ! (field, list_p) return field pointer
private :: initialize          ! ()
private :: make_list           ! (list_p, name) return field pointer

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Public parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer, parameter, public :: fm_field_name_len = 48
! <DATA NAME="fm_field_name_len" TYPE="integer, parameter" DEFAULT="48">
!   The length of a character string representing the field name.
! </DATA>
integer, parameter, public :: fm_path_name_len  = 512
! <DATA NAME="fm_path_name_len" TYPE="integer, parameter" DEFAULT="512">
!   The length of a character string representing the field path.
! </DATA>
integer, parameter, public :: fm_string_len     = 128
! <DATA NAME="fm_string_len" TYPE="integer, parameter" DEFAULT="128">
!   The length of a character string representing character values for the field.
! </DATA>
integer, parameter, public :: fm_type_name_len  = 8
! <DATA NAME="fm_type_name_len" TYPE="integer, parameter" DEFAULT="8">
!   The length of a character string representing the various types that the values of the field can take.
! </DATA>
integer, parameter, public :: NUM_MODELS        = 5
! <DATA NAME="NUM_MODELS" TYPE="integer, parameter" DEFAULT="5">
!   Number of models (ATMOS, OCEAN, LAND, ICE, COUPLER).
! </DATA>
integer, parameter, public :: NO_FIELD          = -1
! <DATA NAME="NO_FIELD" TYPE="integer, parameter" DEFAULT="-1">
!   The value returned if a field is not defined.
! </DATA>! 
integer, parameter, public :: MODEL_ATMOS       = 1
! <DATA NAME="MODEL_ATMOS" TYPE="integer, parameter" DEFAULT="1">
!   Atmospheric model.
! </DATA>! 
integer, parameter, public :: MODEL_OCEAN       = 2
! <DATA NAME="MODEL_OCEAN" TYPE="integer, parameter" DEFAULT="2">
!   Ocean model.
! </DATA>
integer, parameter, public :: MODEL_LAND        = 3
! <DATA NAME="MODEL_LAND" TYPE="integer, parameter" DEFAULT="3">
!   Land model.
! </DATA>
integer, parameter, public :: MODEL_ICE         = 4
! <DATA NAME="MODEL_ICE" TYPE="integer, parameter" DEFAULT="4">
!   Ice model.
! </DATA>
integer, parameter, public :: MODEL_COUPLER     = 5
! <DATA NAME="MODEL_COUPLER" TYPE="integer, parameter" DEFAULT="5">
!   Ice model.
! </DATA>
character(len=11), parameter, public, dimension(NUM_MODELS) :: &
   MODEL_NAMES=(/'atmospheric','oceanic    ','land       ','ice        ','coupler    '/)
! <DATA NAME="MODEL_NAMES" TYPE="character(len=11), parameter">
!   Model names, e.g. MODEL_NAMES(MODEL_OCEAN) is 'oceanic'
! </DATA>

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Public type definitions
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

type, public :: fm_array_list_def  !{
  character (len=fm_field_name_len), dimension(:), pointer :: names => NULL()
  integer                                                  :: length
end type  fm_array_list_def  !}

!
! <TYPE NAME="method_type">
! <DESCRIPTION>

! This method_type is a way to allow a component module to alter the parameters it needs
! for various tracers. In essence this is a way to modify a namelist. A namelist can supply
! default parameters for all tracers. This  method will allow the user to modify these
! default parameters for an individual tracer. An example could be that  the user wishes to
! use second order advection on a tracer and also use fourth order advection on a second
! tracer  within the same model run. The default advection could be second order and the
! field table would then indicate  that the second tracer requires fourth order advection.
! This would be parsed by the advection routine.

!
! </DESCRIPTION>
type, public :: method_type

  ! <DATA NAME="method_type :: method_type" TYPE="character" DIM="(128)">
  !
  !   This string represents a tag that a module using this method can
  !   key on. Typically this should contain some reference to the module
  !   that is calling it.
  ! </DATA>
  !
  ! <DATA NAME="method_type :: method_name" TYPE="character" DIM="(128)">
  !   This is the name of a method which the module can parse and use
  !   to assign different default values to a field method.
  ! </DATA> 
  !
  ! <DATA NAME="method_type :: method_control" TYPE="character" DIM="(256)">
  !   This is the string containing parameters that the module can use
  !   as values  for a field method. These should override default
  !   values within the module.
  ! </DATA>
  character(len=fm_string_len) :: method_type
  character(len=fm_string_len) :: method_name
  character(len=fm_string_len) :: method_control
end type
! </TYPE> NAME="method_type"

! <TYPE NAME="method_type_short">
! <DESCRIPTION>
!   This method_type is the same as method_type except that the
!   method_control string is not present. This is used when you wish to
!   change to a scheme within a module but do not need to pass 
!   parameters.
! </DESCRIPTION>
type, public :: method_type_short
  ! <DATA NAME="method_type_short :: method_type" TYPE="character" DIM="(128)">
  !   see method_type :: method_type above.
  ! </DATA>
  !
  ! <DATA NAME="method_type_short :: method_name" TYPE="character" DIM="(128)">
  !   see method_type :: method_name above.
  ! </DATA> 
  character(len=fm_string_len) :: method_type
  character(len=fm_string_len) :: method_name
end type
! </TYPE> NAME="method_type_short"

! <TYPE NAME="method_type_very_short">
! <DESCRIPTION>
!   This method_type is the same as method_type except that the
!   method_control and method_name strings are not present. This is used
!   when you wish to change to a scheme within a module but do not need
!   to pass  parameters.
! </DESCRIPTION>
type, public :: method_type_very_short
  ! <DATA NAME="method_type_short :: method_type" TYPE="character" DIM="(128)">
  !   see method_type :: method_type above.
  ! </DATA>
  character(len=fm_string_len) :: method_type
end type
! </TYPE> NAME="method_type_very_short"


!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Public types
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

type(method_type), public :: default_method


!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Public variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Interface definitions for overloaded routines
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

interface find_field_index
  module procedure  find_field_index_old
  module procedure  find_field_index_new
end interface

interface parse
  module procedure  parse_real
  module procedure  parse_reals
  module procedure  parse_integer
  module procedure  parse_integers
  module procedure  parse_string
  module procedure  parse_strings
end interface

interface  fm_new_value  !{
  module procedure  fm_new_value_integer
  module procedure  fm_new_value_logical
  module procedure  fm_new_value_real
  module procedure  fm_new_value_string
end interface  !}

interface  fm_get_value  !{
  module procedure  fm_get_value_integer
  module procedure  fm_get_value_logical
  module procedure  fm_get_value_real
  module procedure  fm_get_value_string
end interface  !}

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Private parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=17), parameter :: module_name       = 'field_manager_mod'
character(len=1),  parameter :: bracket_left      = '['
character(len=1),  parameter :: bracket_right     = ']'
character(len=1),  parameter :: comma             = ","
character(len=1),  parameter :: comment           = '#'
character(len=1),  parameter :: dquote            = '"'
character(len=1),  parameter :: equal             = '='
character(len=1),  parameter :: list_sep          = '/'
character(len=1),  parameter :: space             = ' '
character(len=1),  parameter :: squote            = "'"
character(len=1),  parameter :: tab               = char(9) ! ASCII

integer,           parameter :: null_type         = 0
integer,           parameter :: integer_type      = 1
integer,           parameter :: list_type         = 2
integer,           parameter :: logical_type      = 3
integer,           parameter :: real_type         = 4
integer,           parameter :: string_type       = 5
integer,           parameter :: num_types         = 5
integer,           parameter :: line_len          = 256
integer,           parameter :: array_increment   = 10
integer,           parameter :: MAX_FIELDS        = MAXFIELDS_
integer,           parameter :: MAX_FIELD_METHODS = MAXFIELDMETHODS_


!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Private type definitions
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

type, private :: field_mgr_type !{
  character(len=fm_field_name_len)                    :: field_type
  character(len=fm_string_len)                    :: field_name
  integer                                             :: model, num_methods
  type(method_type)                                   :: methods(MAX_FIELD_METHODS)
end type field_mgr_type !}

type, private :: field_names_type !{
  character(len=fm_field_name_len)                    :: fld_type
  character(len=fm_field_name_len)                    :: mod_name
  character(len=fm_string_len)                    :: fld_name
end  type field_names_type !}

type, private :: field_names_type_short !{
  character(len=fm_field_name_len)                    :: fld_type
  character(len=fm_field_name_len)                    :: mod_name
end type field_names_type_short !}

type, private :: field_def  !{
  character (len=fm_field_name_len)                   :: name
  integer                                             :: index
  type (field_def), pointer                           :: parent => NULL()
  integer                                             :: field_type
  integer                                             :: length
  integer                                             :: array_dim
  integer                                             :: max_index
  type (field_def), pointer                           :: first_field => NULL()
  type (field_def), pointer                           :: last_field => NULL()
  integer, pointer, dimension(:)                      :: i_value => NULL()
  logical, pointer, dimension(:)                      :: l_value => NULL()
  real, pointer, dimension(:)                         :: r_value => NULL()
  character(len=fm_string_len), pointer, dimension(:) :: s_value => NULL()
  type (field_def), pointer                           :: next => NULL()
  type (field_def), pointer                           :: prev => NULL()
end type field_def  !}

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Private types
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

type(field_mgr_type), private :: fields(MAX_FIELDS)


!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        Private variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=fm_path_name_len)  :: loop_list
character(len=fm_type_name_len)  :: field_type_name(num_types)
character(len=fm_field_name_len) :: save_root_name
! The string set is the set of characters. 
character(len=52)                :: set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
! If a character in the string being parsed matches a character within
! the string set_nonexp then the string being parsed cannot be a number.
character(len=50)                :: set_nonexp = "ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz"
! If a character in the string being parsed matches a character within
! the string setnum then the string may be a number.
character(len=13)                :: setnum     = "0123456789+-."
integer                          :: num_fields         = 0
integer                          :: verb               = 0
integer                          :: verb_level_warn    = 0
integer                          :: verb_level_note    = 0
integer                          :: default_verbosity  = 0
integer                          :: max_verbosity      = 1
type (field_def), pointer        :: loop_list_p        => NULL()
type (field_def), pointer        :: current_list_p     => NULL()
type (field_def), pointer        :: root_p             => NULL()
type (field_def), pointer        :: save_root_parent_p => NULL()
type (field_def), target, save   :: root 


contains

! <SUBROUTINE NAME="field_manager_init">
!   <OVERVIEW>
!     Routine to initialize the field manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine reads from a file containing formatted strings. 
!     These formatted strings contain information on which schemes are
!     needed within various modules. The field manager does not
!     initialize any of those schemes however. It simply holds the
!     information and is queried by the appropriate  module.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call field_manager_init(nfields, table_name)
!   </TEMPLATE>

subroutine field_manager_init(nfields, table_name)

! <OUT NAME="nfields" TYPE="integer">
!   The number of fields.
! </OUT>

integer,                      intent(out), optional :: nfields

! <IN NAME="table_name" TYPE="character, optional"
!     DIM="(len=128)" DEFAULT="field_table">
!   The name of the field table. The default name is field_table.
! </IN>

character(len=fm_string_len), intent(in), optional :: table_name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=18), parameter :: sub_name     = 'field_manager_init'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=1024)              :: record
character(len=fm_path_name_len)  :: control_str
character(len=fm_path_name_len)  :: list_name
character(len=fm_path_name_len)  :: method_name
character(len=fm_path_name_len)  :: name_str
character(len=fm_path_name_len)  :: type_str
character(len=fm_path_name_len)  :: val_name
character(len=fm_string_len)     :: tbl_name
integer                          :: control_array(MAX_FIELDS,3)
integer                          :: endcont
integer                          :: icount
integer                          :: index_list_name
integer                          :: iunit
integer                          :: l
integer                          :: log_unit
integer                          :: ltrec
integer                          :: m
integer                          :: midcont
integer                          :: model
integer                          :: startcont
logical                          :: flag_method
logical                          :: fm_success
type(field_names_type_short)     :: text_names_short
type(field_names_type)           :: text_names
type(method_type_short)          :: text_method_short
type(method_type)                :: text_method
type(method_type_very_short)     :: text_method_very_short



if (module_is_initialized) then
   if(present(nfields)) nfields = num_fields
   return
endif
num_fields = 0
call initialize

call mpp_io_init()

if (.not.PRESENT(table_name)) then
   tbl_name = 'field_table'
else
   tbl_name = trim(table_name)
endif

if (.not. file_exist(trim(tbl_name))) then
!   <ERROR MSG="No field table available, so no fields are being registered." STATUS="NOTE">
!      The field table does not exist.
!   </ERROR>
if (mpp_pe() == mpp_root_pe()) then
  if (verb .gt. verb_level_warn) then
    call mpp_error(NOTE, trim(warn_header)//                       &
         'No field table ('//trim(tbl_name)//') available, so no fields are being registered.')
  endif
endif
if(present(nfields)) nfields = 0
return
endif


call mpp_open(iunit,file=trim(tbl_name), form=MPP_ASCII, action=MPP_RDONLY)
!write_version_number should precede all writes to stdlog from field_manager
call write_version_number (version, tagname)
log_unit = stdlog()
do while (.TRUE.)
   read(iunit,'(a)',end=89,err=99) record
   write( log_unit,'(a)' )record
   if (record(1:1) == "#" ) cycle
   ltrec =  LEN_TRIM(record)
   if (ltrec .le. 0 ) cycle ! Blank line


         icount = 0
         do l= 1, ltrec
            if (record(l:l) == '"' ) then
               icount = icount + 1
            endif
         enddo
!     <ERROR MSG="Too many fields in field table header entry." STATUS="FATAL">
!       There are more that 3 fields in the field table header entry. 
!       The entry should look like <BR/>
!       "Field_Type","Model_Type","Field_Name" <BR/>
!        or<BR/>
!       "Field_Type","Model_Type"
!     </ERROR>
      if (icount > 6 ) then
        call mpp_error(FATAL,trim(error_header)//'Too many fields in field table header entry.'//trim(record))
      endif

         select case (icount)
           case (6)
             read(record,*,end=79,err=79) text_names
             text_names%fld_type = lowercase(trim(text_names%fld_type))
             text_names%mod_name = lowercase(trim(text_names%mod_name))
             text_names%fld_name = lowercase(trim(text_names%fld_name))
           case(4)
! If there is no control string then the last string can be omitted and there are only 4 '"' in the record.
             read(record,*,end=79,err=79) text_names_short
             text_names%fld_type = lowercase(trim(text_names_short%fld_type))
             text_names%mod_name = lowercase(trim(text_names_short%mod_name))
             text_names%fld_name = lowercase(trim(text_names_short%mod_name))
           case(2)
! If there is only the method_type string then the last 2 strings need to be blank and there are only 2 '"' in the record.
             read(record,*,end=79,err=79) text_names_short
             text_names%fld_type = lowercase(trim(text_names_short%fld_type))
             text_names%mod_name = lowercase(trim(text_names_short%mod_name))
             text_names%fld_name = lowercase(trim(text_names_short%mod_name))
           case default
!     <ERROR MSG="Unterminated field in field table header entry." STATUS="FATAL">
!       There is an unterminated or unquoted string in the field table entry.
             text_names%fld_type = " "
             text_names%mod_name = lowercase(trim(record))
             text_names%fld_name = " "
!             call mpp_error(FATAL,trim(error_header)//'Unterminated field in field_table header entry.'//trim(record))
!     </ERROR>
         end select    

! Create a list with Rick Slaters field manager code

   list_name = list_sep//trim(text_names%mod_name)//list_sep//trim(text_names%fld_type)//&
               list_sep//trim(text_names%fld_name)
   if (mpp_pe() == mpp_root_pe() ) then
     if (verb .gt. verb_level_note) then
!   <ERROR MSG="Creating list name = list_name." STATUS="NOTE">
!      A field is being created called list_name.
!   </ERROR>
       call mpp_error(NOTE, trim(note_header)//'Creating list name = '//trim(list_name))
     endif
   endif

   index_list_name = fm_new_list(list_name, create = .true.)
!   <ERROR MSG="Could not set field list for list_name." STATUS="FATAL">
!      A field called list_name could not be created.
!   </ERROR>
   if ( index_list_name == NO_FIELD ) &
     call mpp_error(FATAL, trim(error_header)//'Could not set field list for '//trim(list_name))

   fm_success = fm_change_list(list_name)  
   select case (text_names%mod_name)
   case ('coupler_mod')
      model = MODEL_COUPLER
   case ('atmos_mod')
      model = MODEL_ATMOS
   case ('ocean_mod')
      model = MODEL_OCEAN
   case ('land_mod')
      model = MODEL_LAND
   case ('ice_mod')
      model = MODEL_ICE
   case default
!   <ERROR MSG="The model name is unrecognised : model_name" STATUS="FATAL">
!      The model name being supplied in the field entry is unrecognised.
!      This should be the second string in the first line of the field entry.
!      Recognised names are atmos_mod, ice_mod, land_mod and ocean_mod.
!   </ERROR>
     call mpp_error(FATAL, trim(error_header)//'The model name is unrecognised : '//trim(text_names%mod_name))
   end select
   if (find_field_index(list_name) > 0) then
      num_fields = num_fields + 1


!     <ERROR MSG="max fields exceeded" STATUS="FATAL">
!       Maximum number of fields for this module has been exceeded.
!     </ERROR>
      if (num_fields > MAX_FIELDS) call mpp_error(FATAL,trim(error_header)//'max fields exceeded')
      fields(num_fields)%model       = model
      fields(num_fields)%field_name  = lowercase(trim(text_names%fld_name))
      fields(num_fields)%field_type  = lowercase(trim(text_names%fld_type))
      fields(num_fields)%num_methods = 0
      call check_for_name_duplication

! Check to see that the first line is not the only line
      if ( record(LEN_TRIM(record):LEN_TRIM(record)) == list_sep) cycle

      flag_method = .TRUE.
      m = 1
      do while (flag_method)
         read(iunit,'(a)',end=99,err=99) record
! If the line is blank then fetch the next line.
         if (LEN_TRIM(record) .le. 0) cycle
! If the last character in the line is / then this is the end of the field methods
         if ( record(LEN_TRIM(record):LEN_TRIM(record)) == list_sep) then
            flag_method = .FALSE.
            if (LEN_TRIM(record) == 1) cycle
            record = record(:LEN_TRIM(record)-1) ! Remove the end of field method marker
         endif
! If the line is now blank, after removing the field separator marker, then fetch the next line.
         if (LEN_TRIM(record) .le. 0) cycle
! If the first character in the line is # then it is treated as a comment
         if (record(1:1) == comment ) cycle

         icount = 0
         do l= 1, LEN_TRIM(record)
            if (record(l:l) == dquote ) then
               icount = icount + 1
            endif
         enddo     
!     <ERROR MSG="Too many fields in field entry." STATUS="FATAL">
!       There are more that 3 fields in the tracer entry. This is probably due
!       to separating the parameters entry into multiple strings. 
!       The entry should look like <BR/>       
!       "Type","Name","Control1=XXX,Control2=YYY" <BR/>
!        and not like<BR/>
!       "Type","Name","Control1=XXX","Control2=YYY"
!     </ERROR>
      if (icount > 6 ) call mpp_error(FATAL,trim(error_header)//'Too many fields in field entry.'//trim(record))

      if (.not. fm_change_list ( list_name)) &
         call mpp_error(FATAL, trim(error_header)//'Could not change to '//trim(list_name)//' list')

      select case (icount)
        case (6)
          read(record,*,end=99,err=99) text_method
          fields(num_fields)%methods(m)%method_type = lowercase(trim(text_method%method_type))
          fields(num_fields)%methods(m)%method_name = lowercase(trim(text_method%method_name))
          fields(num_fields)%methods(m)%method_control = lowercase(trim(text_method%method_control))

          type_str    = text_method%method_type
          name_str    = text_method%method_name
          control_str = text_method%method_control

        case(4)
! If there is no control string then the last string can be omitted and there are only 4 '"' in the record.
          read(record,*,end=99,err=99) text_method_short
          fields(num_fields)%methods(m)%method_type = lowercase(trim(text_method_short%method_type))
          fields(num_fields)%methods(m)%method_name = lowercase(trim(text_method_short%method_name))
          fields(num_fields)%methods(m)%method_control = " "

          type_str    = text_method_short%method_type
          name_str    = ""
          control_str = text_method_short%method_name

        case(2)
! If there is only the method_type string then the last 2 strings need to be blank and there are only 2 '"' in the record.
          read(record,*,end=99,err=99) text_method_very_short
          fields(num_fields)%methods(m)%method_type = lowercase(trim(text_method_very_short%method_type))
          fields(num_fields)%methods(m)%method_name = " "
          fields(num_fields)%methods(m)%method_control = " "

          type_str    = ""
          name_str    = ""
          control_str = text_method_very_short%method_type

        case(0)
          read(record,'(A)',end=99,err=99) control_str
          type_str = ""
          name_str = ""

        case default
!     <ERROR MSG="Unterminated field in field entry." STATUS="FATAL">
!       There is an unterminated or unquoted string in the field table entry.
          call mpp_error(FATAL,trim(error_header)//'Unterminated field in field entry.'//trim(record))
!     </ERROR>
      end select

! This section of code breaks the control string into separate strings. 
! The array control_array contains the following parameters.
! control_array(:,1) = index within control_str of the first character of the name.
! control_array(:,2) = index within control_str of the equal sign
! control_array(:,3) = index within control_str of the last character of the value.
! 
! control_array(:,1)   -> control_array(:,2) -1 = name of the parameter.
! control_array(:,2)+1 -> control_array(:,3)    = value of the parameter.

      ltrec= len_trim(control_str)
      control_array(:,1) = 1
      control_array(:,2:3) = ltrec
      icount = 0
      do l= 1, ltrec
         if (control_str(l:l) == equal ) then
            icount = icount + 1
            control_array(icount,2) = l ! Middle of string
         elseif (control_str(l:l) == comma ) then
            if (icount .eq. 0) then

!     <ERROR MSG="Unterminated field in field entry." STATUS="FATAL">
!       Bad format for field entry (comma without equals sign)
              call mpp_error(FATAL,trim(error_header) //                                &
                   ' Bad format for field entry (comma without equals sign): ''' //     &
                   trim(control_str) // '''')
!     </ERROR>

            elseif (icount .gt. MAX_FIELDS) then

!     <ERROR MSG="Unterminated field in field entry." STATUS="FATAL">
!       Too many fields in field entry
              call mpp_error(FATAL,trim(error_header) //        &
                   ' Too many fields in field entry: ''' //     &
                   trim(control_str) // '''')
!     </ERROR>

            else

              control_array(icount,3) = l-1   !End of previous string
              control_array(min(MAX_FIELDS,icount+1),1) = l+1 !Start of next string

            endif
         endif
      enddo     

      ! Make sure that we point to the end of the string (minus any trailing comma)
      ! for the last set of values. This fixes the case where the last set of values
      ! is a comma separated list

      if (control_str(ltrec:ltrec) .ne. comma) then
        control_array(max(1,icount),3) = ltrec
      endif


      if ( icount == 0 ) then
        method_name = type_str
        if (len_trim(method_name) > 0 ) then
          method_name = trim(method_name)//list_sep// trim(name_str)
        else
          method_name = trim(name_str)
        endif
        val_name = control_str
        
        call new_name(list_name, method_name, val_name )
      
      else
      
        do l = 1,icount
          startcont = control_array(l,1)
          midcont   = control_array(l,2)
          endcont   = control_array(l,3)
          
          method_name = trim(type_str)
          if (len_trim(method_name) > 0 ) then
            method_name = trim(method_name)//list_sep// trim(name_str)
          else
            method_name = trim(name_str)
          endif
          
          if (len_trim(method_name) > 0 ) then
            method_name = trim(method_name)//list_sep//&
                          trim(control_str(startcont:midcont-1))
          else
            method_name = trim(control_str(startcont:midcont-1))
          endif
          val_name =    trim(control_str(midcont+1:endcont))
        
          call new_name(list_name, method_name, val_name )
        enddo

      endif

      fields(num_fields)%num_methods = fields(num_fields)%num_methods + 1
!     <ERROR MSG="Maximum number of methods for field exceeded" STATUS="FATAL">
!       Maximum number of methods allowed for entries in the field table has been exceeded.
!     </ERROR>
      if (fields(num_fields)%num_methods > MAX_FIELD_METHODS) &
         call mpp_error(FATAL,trim(error_header)//'Maximum number of methods for field exceeded')
         m = m + 1
      enddo
   else

!     <ERROR MSG="Field with identical name and model name duplicate found, skipping" STATUS="NOTE">
!       The name of the field and the model name are identical. Skipping that field.
!     </ERROR>
      if (mpp_pe() == 0) then
         if (verb .gt. verb_level_warn) then
           call mpp_error(WARNING, trim(warn_header)//                              &
                'Field with identical name and model name duplicate found, skipping')
          endif
      endif
      flag_method = .TRUE.
      do while (flag_method)
         read(iunit,'(A)',end=99,err=99) record
         if ( record(LEN_TRIM(record):LEN_TRIM(record)) == list_sep) then
            flag_method = .FALSE.
         endif
      enddo
   endif
79 continue
enddo
         
89 continue
close(iunit)

if(present(nfields)) nfields = num_fields
if (verb .gt. verb_level_warn) &
  fm_success= fm_dump_list("/", .true.)
  
default_method%method_type = 'none'
default_method%method_name = 'none'
default_method%method_control = 'none'
return

99 continue

!     <ERROR MSG="error reading field table" STATUS="FATAL">
!       There is an error in reading the field table.
!     </ERROR>
call mpp_error(FATAL,trim(error_header)//' Error reading field table. Record = '//trim(record))

end subroutine field_manager_init
! </SUBROUTINE>

subroutine check_for_name_duplication
integer :: i

! Check that name is unique amoung fields of the same field_type and model.
do i=1,num_fields-1
  if ( fields(i)%field_type == fields(num_fields)%field_type .and. &
       fields(i)%model      == fields(num_fields)%model      .and. &
       fields(i)%field_name == fields(num_fields)%field_name ) then
    if (mpp_pe() .eq. mpp_root_pe()) then
      call mpp_error(WARNING,'Error in field_manager_mod. Duplicate field name: Field type='//trim(fields(i)%field_type)// &
         ',  Model='//trim(MODEL_NAMES(fields(i)%model))// &
         ',  Duplicated name='//trim(fields(i)%field_name))
    endif
  endif
enddo

end subroutine check_for_name_duplication

!#######################################################################
!#######################################################################

! <PRIVATE><SUBROUTINE NAME="new_name">
!   <OVERVIEW>
!     Subroutine to add new values to list parameters.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine uses input strings list_name, method_name
!     and val_name_in to add new values to the list. Given
!     list_name a new list item is created that is named
!     method_name and is given the value or values in
!     val_name_in. If there is more than 1 value in
!     val_name_in, these values should be  comma-separated.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call new_name ( list_name, method_name , val_name_in)
!   </TEMPLATE>
subroutine new_name ( list_name, method_name_in , val_name_in)
!   <IN NAME="list_name" TYPE="character(len=*)">
!     The name of the field that is of interest here.
!   </IN>
!   <IN NAME="method_name" TYPE="character(len=*)">
!     The name of the method that values are being supplied for.
!   </IN>
character(len=*), intent(in)    :: list_name
character(len=*), intent(in)    :: method_name_in
!   <INOUT NAME="val_name_in" TYPE="character(len=*)">
!     The value or values that will be parsed and used as the value when 
!     creating a new field or fields.
!   </INOUT>
character(len=*), intent(inout) :: val_name_in

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=8),  parameter :: sub_name     = 'new_name'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_string_len)   :: method_name
character(len=fm_string_len)   :: val_list
character(len=fm_string_len)   :: val_name
integer, dimension(MAX_FIELDS) :: end_val
integer, dimension(MAX_FIELDS) :: start_val
integer                        :: i
integer                        :: index_t
integer                        :: left_br
integer                        :: num_elem
integer                        :: out_unit
integer                        :: right_br
integer                        :: val_int
integer                        :: val_type
logical                        :: append_new
logical                        :: val_logic
real                           :: val_real
integer                        :: length

call strip_front_blanks(val_name_in)
method_name = trim (method_name_in)
call strip_front_blanks(method_name)

index_t  = 1
num_elem = 1
append_new = .false.
start_val(1) = 1
end_val(:) = len_trim(val_name_in)

! If the array of values being passed in is a comma delimited list then count 
! the number of elements.

do i = 1, len_trim(val_name_in)
  if ( val_name_in(i:i) == comma ) then
    end_val(num_elem) = i-1
    start_val(num_elem+1) = i+1
    num_elem = num_elem + 1
  endif
enddo

! Check to see if this is an array element of form array[x] = value
left_br  = scan(method_name,'[')
right_br = scan(method_name,']')
if ( num_elem .eq. 1 ) then 
!     <ERROR MSG="Left bracket present without right bracket in method_name" STATUS="FATAL">
!       When using an array element an unpaired bracket was found.
!     </ERROR>
  if ( left_br > 0 .and. right_br == 0 ) &
    call mpp_error(FATAL, trim(error_header)//"Left bracket present without right bracket in "//trim(method_name))
!     <ERROR MSG="Right bracket present without left bracket in method_name" STATUS="FATAL">
!       When using an array element an unpaired bracket was found.
!     </ERROR>
  if ( left_br== 0 .and. right_br > 0 ) &
    call mpp_error(FATAL, trim(error_header)//"Right bracket present without left bracket in "//trim(method_name))
  

  if ( left_br > 0 .and. right_br > 0 ) then 
!     <ERROR MSG="Using a non-numeric value for index in method_name" STATUS="FATAL">
!       An array assignment was requested but a non-numeric value was found. i.e. array[a] = 1
!     </ERROR>
    if ( scan( method_name(left_br+1:right_br -1), set ) > 0 ) &
       call mpp_error(FATAL, trim(error_header)//"Using a non-numeric value for index in "//trim(method_name))
    read(method_name(left_br+1:right_br -1), *) index_t
    method_name = method_name(:left_br -1)
  endif
else
! If there are multiple values then there cannot be a bracket in method_name.
!     <ERROR MSG="Using a comma delimited list with an indexed array element in method_name" STATUS="FATAL">
!       When supplying multiple values an index was found. i.e array[3] = 4,5,6 is invalid.
!     </ERROR>
  if ( left_br > 0 .or. right_br > 0 ) &
    call mpp_error(FATAL, &
      trim(error_header)//"Using a comma delimited list with an indexed array element in "//trim(method_name))

endif

do i = 1, num_elem

  if ( i .gt. 1 .or. index_t .eq. 0 ) then
    append_new = .true.
    index_t = 0 ! If append is true then index must be <= 0
  endif  
  val_type = string_type  ! Assume it is a string
  val_name = val_name_in(start_val(i):end_val(i))
  call strip_front_blanks(val_name)


!
!       if the string starts and ends with matching single quotes, then this is a string
!       if there are quotes which do not match, then this is an error
!

  length = len_trim(val_name)
  if (val_name(1:1) .eq. squote) then  !{

    if (val_name(length:length) .eq. squote) then
      val_name = val_name(2:length-1)
      val_type = string_type
    elseif (val_name(length:length) .eq. dquote) then
      call mpp_error(FATAL, trim(error_header) // ' Quotes do not match in ' // trim(val_name) //       &
           ' for ' // trim(method_name) // ' of ' // trim(list_name))
    else
      call mpp_error(FATAL, trim(error_header) // ' No trailing quote in ' // trim(val_name) //         &
           ' for ' // trim(method_name) // ' of ' // trim(list_name))
    endif

  elseif (val_name(1:1) .eq. dquote .or. val_name(length:length) .eq. dquote) then  !}{

    call mpp_error(FATAL, trim(error_header) // ' Double quotes not allowed in ' // trim(val_name) //   &
         ' for ' // trim(method_name) // ' of ' // trim(list_name))

  elseif (val_name(length:length) .eq. squote) then  !}{

    call mpp_error(FATAL, trim(error_header) // ' No leading quote in ' // trim(val_name) //            &
         ' for ' // trim(method_name) // ' of ' // trim(list_name))

  else  !}{
! If the string to be parsed is a real then all the characters must be numeric, 
! be a plus/minus, be a decimal point or, for exponentials, be e or E.

! If a string is an integer, then all the characters must be numeric.

  if ( scan(val_name(1:1), setnum ) > 0 ) then  

! If there is a letter in the name it may only be e or E

      if ( scan(val_name, set_nonexp ) > 0 ) then
        if (verb .gt. verb_level_warn) then
!     <ERROR MSG="First character of value is numerical but the value does not appear to be numerical." STATUS="WARNING">
!       The value may not be numerical. This is a warning as the user may wish to use a value of 2nd_order.
!     </ERROR>
          call mpp_error(WARNING, trim(warn_header)//                                  &
               'First character of value is numerical but the value does not appear to be numerical.')
          call mpp_error(WARNING, 'Name = '// trim(list_name)// list_sep//                &
               trim(method_name)// ' Value = '// trim(val_name))
        endif

      else
! It is real if there is a . in the name or the value appears exponential
        if ( scan(val_name, '.') > 0 .or. scan(val_name, 'e') > 0 .or. scan(val_name, 'E') > 0) then 
          read(val_name, *) val_real
          val_type = real_type
        else
          read(val_name, *) val_int
          val_type = integer_type
        endif   
      endif

    endif

! If val_name is t/T or f/F then this is a logical flag.
    if ( len_trim(val_name) == 1 .or. len_trim(val_name) == 3) then
       if ( val_name == 't' .or. val_name == 'T' .or. val_name == '.t.' .or. val_name == '.T.' ) then
         val_logic = .TRUE.
         val_type = logical_type
       endif
       if ( val_name == 'f' .or. val_name == 'F' .or. val_name == '.f.' .or. val_name == '.F.' ) then
         val_logic = .FALSE.
         val_type = logical_type
       endif
    endif
    if ( trim(lowercase(val_name)) == 'true' .or. trim(lowercase(val_name)) == '.true.' ) then
      val_logic = .TRUE.
      val_type = logical_type
    endif
    if ( trim(lowercase(val_name)) == 'false' .or. trim(lowercase(val_name)) == '.false.' ) then
      val_logic = .FALSE.
      val_type = logical_type
    endif
  endif  !}

  select case(val_type) 

    case (integer_type)
      if ( fm_new_value( method_name, val_int, create = .true., index = index_t, append = append_new ) < 0 ) &
        call mpp_error(FATAL, trim(error_header)//'Could not set "' // trim(val_name) // '" for '//trim(method_name)//&
                              ' (I) for '//trim(list_name))

    case (logical_type)
      if ( fm_new_value( method_name, val_logic, create = .true., index = index_t, append = append_new) < 0 ) &
        call mpp_error(FATAL, trim(error_header)//'Could not set "' // trim(val_name) // '" for '//trim(method_name)//&
                              ' (L) for '//trim(list_name))

    case (real_type)
      if ( fm_new_value( method_name, val_real, create = .true., index = index_t, append = append_new) < 0 ) &
        call mpp_error(FATAL, trim(error_header)//'Could not set "' // trim(val_name) // '" for '//trim(method_name)//&
                              ' (R) for '//trim(list_name))

    case (string_type)
      if ( fm_new_value( method_name, val_name, create = .true., index = index_t, append = append_new) < 0 ) &
        call mpp_error(FATAL, trim(error_header)//'Could not set "' // trim(val_name) // '" for '//trim(method_name)//&
                              ' (S) for '//trim(list_name))
    case default
      call mpp_error(FATAL, trim(error_header)//'Could not find a valid type to set the '//trim(method_name)//&
                            ' for '//trim(list_name))
    
  end select

  if (mpp_pe() == mpp_root_pe() ) then
    if (verb .gt. verb_level_note) then
      out_unit = stdout()
      write (out_unit,*) trim(note_header), 'Creating new value = ', trim(method_name), ' ', trim(val_name)
    endif
  endif

enddo

end subroutine new_name 
!</SUBROUTINE>
!</PRIVATE>
!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="field_manager_end">
!   <OVERVIEW>
!     Destructor for field manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine writes to the logfile that the user is exiting field_manager and 
!     changes the initialized flag to false.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call field_manager_end
!   </TEMPLATE>
subroutine field_manager_end

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=17), parameter :: sub_name     = 'field_manager_end'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

integer :: unit 

call write_version_number (version, tagname)
if ( mpp_pe() == mpp_root_pe() ) then
   unit = stdlog()
   write (unit,'(/,(a))') trim(note_header), 'Exiting field_manager, have a nice day ...'
   unit = stdout()
   write (unit,'(/,(a))') trim(note_header), 'Exiting field_manager, have a nice day ...'
endif

module_is_initialized = .false.

end subroutine field_manager_end
! </SUBROUTINE>

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="strip_front_blanks">
!   <OVERVIEW>
!     A routine to strip whitespace from the start of character strings.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine removes spaces and tabs from the start of a character string.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call strip_front_blanks(name)
!   </TEMPLATE>
subroutine strip_front_blanks(name)

character(len=*), intent(inout) :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=18), parameter :: sub_name     = 'strip_front_blanks'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

integer :: i, j

j = 1
do i = 1,len_trim(name) !{
   if ( .not. (name(i:i) .eq. space .or.                        &
               name(i:i) .eq. tab)) then  !{
    j = i
    exit
  endif !}
enddo !}
name = name(j:)
end subroutine strip_front_blanks
!</SUBROUTINE>

!#######################################################################
!#######################################################################

! <FUNCTION NAME="find_field_index">
!   <OVERVIEW>
!     Function to return the index of the field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This function when passed a model number and a field name will 
!     return the index of the field within the field manager. This index 
!     can be used to access other information from the field manager.
!   </DESCRIPTION>
!   <TEMPLATE>
!     value=find_field_index( model, field_name )
!     value=find_field_index( field_name )
!   </TEMPLATE>

function find_field_index_old(model, field_name)
! 
!   <IN NAME="model" TYPE="integer">
!     The number indicating which model is used.
!   </IN>
!   <IN NAME="field_name" TYPE="character">
!     The name of the field that an index is being requested for.
!   </IN>
!   <OUT NAME="find_field_index" TYPE="integer">
!     The index of the field corresponding to field_name.
!   </OUT>

integer                      :: find_field_index_old
integer,          intent(in) :: model
character(len=*), intent(in) :: field_name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=16), parameter :: sub_name     = 'find_field_index'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
integer :: i

find_field_index_old = NO_FIELD

do i=1,num_fields
   if (fields(i)%model == model .and. fields(i)%field_name == lowercase(field_name)) then
      find_field_index_old = i
      return
   endif
enddo

end function find_field_index_old

function find_field_index_new(field_name)
! 
!   <IN NAME="field_name" TYPE="character">
!     The path to the name of the field that an index is being requested for.
!   </IN>
!   <OUT NAME="find_field_index" TYPE="integer">
!     The index of the field corresponding to field_name.
!   </OUT>

integer                      :: find_field_index_new
character(len=*), intent(in) :: field_name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=16), parameter :: sub_name     = 'find_field_index'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
integer :: i

find_field_index_new = NO_FIELD

find_field_index_new = fm_get_index(field_name)

end function find_field_index_new
! </FUNCTION>

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="get_field_info">
!   <OVERVIEW>
!     This routine allows access to field information given an index.
!   </OVERVIEW>
!   <DESCRIPTION>
!     When passed an index, this routine will return the type of field, 
!     the name of the field, the model which the field is associated and 
!     the number of methods associated with the field.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_field_info( n,fld_type,fld_name,model,num_methods )
!   </TEMPLATE>
subroutine get_field_info(n,fld_type,fld_name,model,num_methods)
!
!   <IN NAME="n" TYPE="integer">
!     The field index.
!   </IN>
integer,          intent(in)  :: n

!   <OUT NAME="fld_type" TYPE="character" DIM="(*)">
!     The field type.
!   </OUT>

!   <OUT NAME="fld_name" TYPE="character" DIM="(*)">
!     The name of the field.
!   </OUT>

!   <OUT NAME="model" TYPE="integer">
!     The number indicating which model is used.
!   </OUT>

!   <OUT NAME="num_methods" TYPE="integer">
!     The number of methods.
!   </OUT>
character (len=*),intent(out) :: fld_type, fld_name
integer, intent(out) :: model, num_methods

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=14), parameter :: sub_name     = 'get_field_info'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

!   <ERROR MSG="invalid field index" STATUS="FATAL">
!     The field index is invalid because it is less than 1 or greater than the 
!     number of fields.
!   </ERROR>
if (n < 1 .or. n > num_fields) call mpp_error(FATAL,trim(error_header)//'Invalid field index')

fld_type    = fields(n)%field_type
fld_name    = fields(n)%field_name
model       = fields(n)%model
num_methods = fields(n)%num_methods

end subroutine get_field_info
! </SUBROUTINE>

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="get_field_method">
!   <OVERVIEW>
!     A routine to get a specified method.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine, when passed a field index and a method index will 
!     return the method text associated with the field(n) method(m).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_field_method( n,m,method )
!   </TEMPLATE>
subroutine get_field_method(n,m,method)
!
!   <IN NAME="n" TYPE="integer">
!     The field index.
!   </IN>
!   <IN NAME="m" TYPE="integer">
!     The method index.
!   </IN>
!   <OUT NAME="method" TYPE="type(method_type)">
!     The m-th method of field with index n.
!   </OUT>
integer,           intent(in)    :: n
integer,           intent(in)    :: m
type(method_type) ,intent(inout) :: method

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=16), parameter :: sub_name     = 'get_field_method'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

!   <ERROR MSG="invalid field index" STATUS="FATAL">
!     The field index is invalid because it is less than 1 or greater than the 
!     number of fields.
!   </ERROR>
if (n < 1 .or. n > num_fields) call mpp_error(FATAL,trim(error_header)//'Invalid field index')

!   <ERROR MSG="invalid method index" STATUS="FATAL">
!     The method index is invalid because it is less than 1 or greater than 
!     the number of methods.
!   </ERROR>
if (m < 1 .or. m > fields(n)%num_methods) call mpp_error(FATAL,trim(error_header)//'Invalid method index')

  method = fields(n)%methods(m)

end subroutine get_field_method
! </SUBROUTINE>

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="get_field_methods">
!   <OVERVIEW>
!     A routine to obtain all the methods associated with a field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     When passed a field index, this routine will return the text 
!     associated with all the methods attached to the field.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_field_methods( n,methods )
!   </TEMPLATE>
subroutine get_field_methods(n,methods)
!
!   <IN NAME="n" TYPE="integer">
!     The field index.
!   </IN>
!   <OUT NAME="method" TYPE="type(method_type)" DIM="(:)">
!     An array of methods for field with index n.
!   </OUT>
integer,          intent(in)  :: n

type(method_type),intent(inout) :: methods(:)

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=17), parameter :: sub_name     = 'get_field_methods'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len), dimension(size(methods(:))) :: control
character(len=fm_path_name_len), dimension(size(methods(:))) :: method
logical                                                   :: found_methods
!   <ERROR MSG="invalid field index" STATUS="FATAL">
!     The field index is invalid because it is less than 1 or greater than the 
!     number of fields.
!   </ERROR>
  if (n < 1 .or. n > num_fields) &
    call mpp_error(FATAL,trim(error_header)//'Invalid field index')

!   <ERROR MSG="method array too small" STATUS="FATAL">
!     The method array is smaller than the number of methods.
!   </ERROR>
  if (size(methods(:)) <  fields(n)%num_methods) &
    call mpp_error(FATAL,trim(error_header)//'Method array too small')

  methods = default_method
  methods(1:fields(n)%num_methods) = fields(n)%methods(1:fields(n)%num_methods)

end subroutine get_field_methods
! </SUBROUTINE>

!#######################################################################
!#######################################################################
  
! <FUNCTION NAME="parse">
!   <OVERVIEW>
!     A function to parse an integer or an array of integers, 
!     a real or an array of reals, a string or an array of strings.
!   </OVERVIEW>
!   <DESCRIPTION>
!  Parse is an integer function that decodes values from a text string.
!  The text string has the form: "label=list" where "label" is an
!  arbitrary user defined label describing the values being decoded,
!  and "list" is a list of one or more values separated by commas.
!  The values may be integer, real, or character.
!  Parse returns the number of values decoded.
!   </DESCRIPTION>
!   <TEMPLATE>
!     number = parse(text, label, value)
!   </TEMPLATE>


function parse_reals ( text, label, values ) result (parse)
!
!   <IN NAME="text" TYPE="character(len=*)">
!     The text string from which the values will be parsed.
!   </IN>
!   <IN NAME="label" TYPE="character(len=*)">
!     A label which describes the values being decoded. 
!   </IN>
!   <OUT NAME="value" TYPE="integer, real, character(len=*)">
!     The value or values that have been decoded.
!   </OUT>
!   <OUT NAME="parse" TYPE="integer">
!     The number of values that have been decoded. This allows 
!     a user to define a large array and fill it partially with 
!     values from a list. This should be the size of the value array.
!   </OUT>
character(len=*), intent(in)  :: text, label
real,             intent(out) :: values(:)

include 'parse.inc'
end function parse_reals
! </FUNCTION>

!#######################################################################
!#######################################################################

function parse_integers ( text, label, values ) result (parse)
character(len=*), intent(in)  :: text, label
integer,          intent(out) :: values(:)

include 'parse.inc'
end function parse_integers

!#######################################################################
!#######################################################################

function parse_strings ( text, label, values ) result (parse)
character(len=*), intent(in)  :: text, label
character(len=*), intent(out) :: values(:)

include 'parse.inc'
end function parse_strings

!#######################################################################
!#######################################################################

!---- scalar overloads -----

function parse_real ( text, label, value ) result (parse)
character(len=*), intent(in)  :: text, label
real,             intent(out) :: value
integer :: parse

real :: values(1)

   parse = parse_reals ( text, label, values )
   if (parse > 0) value = values(1)
end function parse_real

!#######################################################################
!#######################################################################

function parse_integer ( text, label, value ) result (parse)
character(len=*), intent(in)  :: text, label
integer,          intent(out) :: value
integer :: parse

integer :: values(1)

   parse = parse_integers ( text, label, values )
   if (parse > 0) value = values(1)
end function parse_integer

!#######################################################################
!#######################################################################

function parse_string ( text, label, value ) result (parse)
character(len=*), intent(in)  :: text, label
character(len=*), intent(out) :: value
integer :: parse

character(len=len(value)) :: values(1)

   parse = parse_strings ( text, label, values )
   if (parse > 0) value = values(1)
end function parse_string

!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="create_field">
!
! <OVERVIEW>
!    A function to create a field as a child of parent_p. This will return
!    a pointer to a field_def type.
! </OVERVIEW>
! <DESCRIPTION>
!    Allocate and initialize a new field in parent_p list.
!    Return a pointer to the field on success, or a null pointer
!    on failure.
! </DESCRIPTION>
!   <TEMPLATE>
!     list_p => create_field(parent_p, name)
!   </TEMPLATE>
!
!
function  create_field(parent_p, name)                        &
          result (list_p)  !{
!
!   <IN NAME="parent_p" TYPE="type(field_def), pointer">
!     A pointer to the parent of the field that is to be created.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     The name of the field that is to be created.
!   </IN>
!   <OUT NAME="list_p" TYPE="type(field_def), pointer">
!     A pointer to the field that has been created.
!   </OUT>
!
!        Function definition
!
type (field_def), pointer    :: list_p
!
!        arguments
!
type (field_def), pointer    :: parent_p
character(len=*), intent(in) :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=12), parameter :: sub_name     = 'create_field'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
integer                      :: ier
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                      :: error, out_unit
!
!        Check for fatal errors which should never arise
!
out_unit = stdout()
if (.not. associated(parent_p)) then  !{

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Unnassociated pointer'  &
                   , ' for ', trim(name)
  endif  !}
  nullify(list_p)
  return
endif  !}

if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Empty name for '        &
                   , trim(name)
  endif  !}
  nullify(list_p)
  return
endif  !}
!
!        Allocate space for the new list
!
allocate(list_p, stat = error)
if (error .ne. 0) then !{
  write (out_unit,*) trim(error_header), 'Error ', error,       &
       ' allocating memory for list ', trim(name)
  nullify(list_p)
  return
endif  !}
!
!        Initialize the new field
!
list_p%name = name

nullify(list_p%next)
list_p%prev => parent_p%last_field
nullify(list_p%first_field)
nullify(list_p%last_field)
list_p%length = 0
list_p%field_type = null_type
list_p%max_index = 0
list_p%array_dim = 0
if (associated(list_p%i_value)) deallocate(list_p%i_value)
if (associated(list_p%l_value)) deallocate(list_p%l_value)
if (associated(list_p%r_value)) deallocate(list_p%r_value)
if (associated(list_p%s_value)) deallocate(list_p%s_value)
!
!        If this is the first field in the parent, then set the pointer
!        to it, otherwise, update the "next" pointer for the last list
!
if (parent_p%length .le. 0) then  !{
  parent_p%first_field => list_p
else  !}{
  parent_p%last_field%next => list_p
endif  !}
!
!        Update the pointer for the last list in the parent
!
parent_p%last_field => list_p
!
!        Update the length for the parent
!
parent_p%length = parent_p%length + 1
!
!        Set the new index as the return value
!
list_p%index = parent_p%length
!
!        set the pointer to the parent list
!
list_p%parent => parent_p

end function  create_field  !}
! </FUNCTION> NAME="create_field"
!</PRIVATE>
!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="dump_list">
!
! <OVERVIEW>
!    This is a function that lists the parameters of a field.
! </OVERVIEW>
! <DESCRIPTION>
!    Given a pointer to a list, this function prints out the fields, and 
!    subfields, if recursive is true, associated with the list.
!
!    This is most likely to be used through fm_dump_list.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = dump_list(list_p, recursive= .true., depth=0)
!   </TEMPLATE>
!
recursive function dump_list(list_p, recursive, depth)                &
          result (success)  !{
!
!   <IN NAME="list_p" TYPE="type(field_def), pointer">
!     A pointer to the field, the contents of which will be printed out.
!   </IN>
!   <IN NAME="recursive" TYPE="logical">
!     A flag to make the function recursively print all the sub-fields 
!     of the field pointed to by list_p.
!   </IN>
!   <IN NAME="depth" TYPE="integer">
!     The listing will be padded so that 'depth' spaces appear before 
!     the field being printed.
!   </IN>
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!
!        Function definition
!
logical                             :: success
!
!        arguments
!
type (field_def), pointer           :: list_p
logical, intent(in)                 :: recursive
integer, intent(in)                 :: depth

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer,                  parameter :: max_depth    = 128
character(len=max_depth), parameter :: blank        = '    '
character(len=9),  parameter :: sub_name     = 'dump_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                             :: depthp1
integer                             :: first
integer                             :: i
integer                             :: j
integer                             :: last
integer                             :: nf
integer                             :: nl
integer                             :: out_unit
character(len=fm_field_name_len)    :: num
character(len=fm_field_name_len)    :: scratch
type (field_def), pointer           :: this_field_p
!
!        Check for a valid list
!

out_unit = stdout()
this_field_p => NULL()

if (.not. associated(list_p)) then  !{

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Invalid list pointer'
  endif  !}
  success = .false.
elseif (list_p%field_type .ne. list_type) then  !}{

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),               &
                       trim(list_p%name), ' is not a list'
  endif  !}
  success = .false.
else  !}{
!
!        set the default return value
!
  success = .true.
!
!        Print the name of this list
!
  write (out_unit,'(a,a,a)') blank(1:depth), trim(list_p%name), list_sep
!
!        Increment the indentation depth
!
  if (depth .eq. max_depth) then  !{
    if (verb .gt. verb_level_note) then  !{
      write (out_unit,*) trim(note_header),                        &
          'Indentation depth exceeded'
    endif  !}
  else  !}{
    ! The following max function is to work around an error in the IBM compiler for len_trim
    depthp1 = depth + max(len_trim(list_p%name),0) + len_trim(list_sep)
  endif  !}

  this_field_p => list_p%first_field

  do while (associated(this_field_p))  !{

    select case(this_field_p%field_type)
    case(list_type)
!
!        If this is a list, then call dump_list
!
      if (recursive) then  !{
! If recursive is true, then this routine will find and dump sub-fields.
        if (.not. dump_list(this_field_p, .true., depthp1)) then  !{
          success = .false.
          exit
        endif  !}
      else  !}{ ! Otherwise it will print out the name of this field.
        write (out_unit,'(a,a,a)') blank(1:depthp1),               &
                trim(this_field_p%name), list_sep
      endif  !}

    case(integer_type)

         if (this_field_p%max_index .eq. 0) then  !{
         ! Write out the solitary value for this field.
          write (out_unit,'(a,a,a)') blank(1:depthp1),             &
               trim(this_field_p%name), ' = NULL'
        elseif (this_field_p%max_index .eq. 1) then  !}{
          write (scratch,*) this_field_p%i_value(1)
          call strip_front_blanks(scratch)
          write (out_unit,'(a,a,a,a)') blank(1:depthp1),           &
                trim(this_field_p%name), ' = ', trim(scratch)

        else  !}{ Write out the array of values for this field.
          do j = 1, this_field_p%max_index - 1  !{
            write (scratch,*) this_field_p%i_value(j)
            call strip_front_blanks(scratch)
            write (num,*) j
            call strip_front_blanks(num)
            write (out_unit,'(a,a,a,a,a,a)') blank(1:depthp1),     &
                 trim(this_field_p%name), '[', trim(num),          &
                 '] = ', trim(scratch)
          enddo  !} j
          write (scratch,*) this_field_p%i_value(this_field_p%max_index)
          call strip_front_blanks(scratch)
          write (num,*) this_field_p%max_index
          call strip_front_blanks(num)
          write (out_unit,'(a,a,a,a,a,a)') blank(1:depthp1),       &
               trim(this_field_p%name), '[', trim(num),            &
               '] = ', trim(scratch)
        endif  !}



    case(logical_type)

        if (this_field_p%max_index .eq. 0) then  !{
         ! Write out the solitary value for this field.
          write (out_unit,'(a,a,a)') blank(1:depthp1),             &
               trim(this_field_p%name), ' = NULL'
        elseif (this_field_p%max_index .eq. 1) then  !}{
          write (out_unit,'(a,a,a,l1)') blank(1:depthp1),          &
               trim(this_field_p%name), ' = ',                     &
               this_field_p%l_value(1)
        else  !}{ Write out the array of values for this field.
          do j = 1, this_field_p%max_index - 1  !{
            write (num,*) j
            call strip_front_blanks(num)
            write (out_unit,'(a,a,a,a,a,l1)') blank(1:depthp1),    &
                 trim(this_field_p%name), '[', trim(num),          &
                 '] = ', this_field_p%l_value(j)
          enddo  !} j
          write (num,*) this_field_p%max_index
          call strip_front_blanks(num)

       write (out_unit,'(a,a,a,a,a,l1)') blank(1:depthp1),         &
               trim(this_field_p%name), '[', trim(num),            &
               '] = ', this_field_p%l_value(this_field_p%max_index)
        endif  !}


    case(real_type)

        if (this_field_p%max_index .eq. 0) then  !{
         ! Write out the solitary value for this field.
          write (out_unit,'(a,a,a)') blank(1:depthp1),             &
               trim(this_field_p%name), ' = NULL'
        elseif (this_field_p%max_index .eq. 1) then  !}{
          write (scratch,*) this_field_p%r_value(1)
          call strip_front_blanks(scratch)
          write (out_unit,'(a,a,a,a)') blank(1:depthp1),           &
                  trim(this_field_p%name), ' = ', trim(scratch)
        else  !}{ Write out the array of values for this field.
          do j = 1, this_field_p%max_index - 1  !{
            write (scratch,*) this_field_p%r_value(j)
            call strip_front_blanks(scratch)
            write (num,*) j
            call strip_front_blanks(num)
            write (out_unit,'(a,a,a,a,a,a)') blank(1:depthp1),     &
                 trim(this_field_p%name), '[', trim(num),          &
                 '] = ', trim(scratch)
          enddo  !} j
          write (scratch,*) this_field_p%r_value(this_field_p%max_index)
          call strip_front_blanks(scratch)
          write (num,*) this_field_p%max_index
          call strip_front_blanks(num)
          write (out_unit,'(a,a,a,a,a,a)') blank(1:depthp1),       &
               trim(this_field_p%name), '[', trim(num),            &
               '] = ', trim(scratch)
        endif  !}

    case(string_type)
        if (this_field_p%max_index .eq. 0) then  !{
         ! Write out the solitary value for this field.
          write (out_unit,'(a,a,a)') blank(1:depthp1),             &
               trim(this_field_p%name), ' = NULL'
        elseif (this_field_p%max_index .eq. 1) then  !}{
        write (out_unit,'(a,a,a,a,a)') blank(1:depthp1),           &
                trim(this_field_p%name), ' = ''',                  &
               trim(this_field_p%s_value(1)), ''''
        else  !}{ Write out the array of values for this field.
          do j = 1, this_field_p%max_index - 1  !{
            write (num,*) j
            call strip_front_blanks(num)
            write (out_unit,'(a,a,a,a,a,a,a)') blank(1:depthp1),   &
                 trim(this_field_p%name), '[', trim(num),          &
                 '] = ''', trim(this_field_p%s_value(j)), ''''
          enddo  !} j
          write (num,*) this_field_p%max_index
          call strip_front_blanks(num)
          write (out_unit,'(a,a,a,a,a,a,a)') blank(1:depthp1),     &
               trim(this_field_p%name), '[', trim(num),            &
               '] = ''',                                           &
               trim(this_field_p%s_value(this_field_p%max_index)), &
               ''''
        endif  !}

    case default

        if (verb .gt. verb_level_warn) then  !{
          write (out_unit,*) trim(warn_header),                    &
                  'Undefined type for ',                           &
                  trim(this_field_p%name)
        endif  !}
        success = .false.
        exit

    end select

    this_field_p => this_field_p%next
  enddo  !}
endif  !}

end function dump_list  !}
! </FUNCTION> NAME="dump_list"
!</PRIVATE>
!#######################################################################
!#######################################################################

! <PRIVATE><SUBROUTINE NAME="find_base">
!
! <OVERVIEW>
!    A subroutine that splits a listname into a path and a base.
! </OVERVIEW>
! <DESCRIPTION>
!    Find the base name for a list by splitting the list name into
!    a path and base. The base is the last field within name, while the
!    path is the preceding section of name. The base string can then be 
!    used to query for values associated with name.
! </DESCRIPTION>
!   <TEMPLATE>
!     call find_base(name, path, base)
!   </TEMPLATE>
!
subroutine find_base(name, path, base)  !{
!
!   <IN NAME="name" TYPE="character(len=*)">
!   </IN>
!   <OUT NAME="path" TYPE="character(len=*)">
!      A string containing the path of the base field.
!   </OUT>
!   <OUT NAME="base" TYPE="character(len=*)">
!      A string which can be used to query for values associated with name.
!   </OUT>
!
!        arguments
!
character(len=*), intent(in)  :: name
character(len=*), intent(out) :: path
character(len=*), intent(out) :: base

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=9),  parameter :: sub_name     = 'find_base'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

integer :: i
integer :: length

!
!        Check for the last occurrence of the list separator in name
!
! The following max function is to work around an error in the IBM compiler for len_trim
length = max(len_trim(name),0)

if (length .eq. 0) then  !{

   !
   !       Empty name, so return empty path and base
   !
   path = ' '
   base = ' '
else  !}{
   !
   !       Remove trailing list separators
   !
   do while (name(length:length) .eq. list_sep)  !{
      length = length - 1
      if (length .eq. 0) then  !{
         exit
      endif  !}
   enddo  !}
   if (length .eq. 0) then  !{

      !
      !       Name only list separators, so return empty path and base
      !
      path = ' '
      base = ' '
   else  !}{
      !
      !       Check for the last occurrence of the list separator in name
      !
      i = index(name(1:length), list_sep, back = .true.)
      if (i .eq. 0) then  !{
         !
         !       no list separators in the path, so return an empty path
         !       and name as the base
         !
         path = ' '
         base = name(1:length)
      else  !}{
         !
         !       Found a list separator, so return the part up to the last
         !       list separator in path, and the remainder in base
         !
         path = name(1:i)
         base = name(i+1:length)
      endif  !}
   endif  !}
endif  !}

end subroutine find_base  !}
! </SUBROUTINE> NAME="find_base"
!</PRIVATE>
!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="find_field">
!
! <OVERVIEW>
!    Find and return a pointer to the field in the specified
!    list. Return a null pointer on error.
! </OVERVIEW>
! <DESCRIPTION>
!    Find and return a pointer to the field in the specified
!    list. Return a null pointer on error. Given a pointer to a field, 
!    this function searchs for "name" as a sub field.
! </DESCRIPTION>
!   <TEMPLATE>
!     field_p => find_field(name, this_list_p)
!   </TEMPLATE>
!
function find_field(name, this_list_p)                                &
        result (field_p)  !{
!  <OUT NAME="field_p" TYPE="type(field_def), pointer">
!    A pointer to the field corresponding to "name" or an unassociated 
!    pointer if the field name does not exist.
!  </OUT>
!  <IN NAME="name" TYPE="character(len=*)">
!    The name of a field that the user wishes to find.
!  </IN>
!  <IN NAME="this_list_p" TYPE="type(field_def), pointer">
!    A pointer to a list which the user wishes to search for a field "name".
!  </IN>
!
!        Function definition
!
type (field_def), pointer    :: field_p
!
!        arguments
!
character(len=*), intent(in) :: name
type (field_def), pointer    :: this_list_p

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=10), parameter :: sub_name     = 'find_field'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save    :: temp_p 


nullify (field_p)

if (name .eq. '.') then  !{

!
!        If the field is '.' then return this list
!
  field_p => this_list_p
elseif (name .eq. '..') then  !}{
!
!        If the field is '..' then return the parent list
!
  field_p => this_list_p%parent
else  !}{
!
!        Loop over each field in this list
!
  temp_p => this_list_p%first_field

  do while (associated(temp_p))  !{
!
!        If the name matches, then set the return pointer and exit
!        the loop
!
    if (temp_p%name .eq. name) then  !{
      field_p => temp_p
      exit
    endif  !}

    temp_p => temp_p%next

  enddo  !}
endif  !}

end function find_field  !}
! </FUNCTION> NAME="find_field"
!</PRIVATE>

!#######################################################################
!#######################################################################

! <PRIVATE><SUBROUTINE NAME="find_head">
!
! <OVERVIEW>
!    Find the first list for a name by splitting the name into
!    a head and the rest.
! </OVERVIEW>
! <DESCRIPTION>
!   Find the first list for a name by splitting the name into a head and the
! rest. The head is the first field within name, while rest is the remaining
! section of name. The head string can then be used to find other fields that
! may be associated with name.
! </DESCRIPTION>
!   <TEMPLATE>
!     call find_head(name, head, rest)
!   </TEMPLATE>
!
subroutine find_head(name, head, rest)  !{
!
!   <IN NAME="name" TYPE="character(len=*)">
!      The name of a field of interest.
!   </IN>
!   <OUT NAME="head" TYPE="character(len=*)">
!      head is the first field within name.
!   </OUT>
!   <OUT NAME="rest" TYPE="character(len=*)">
!      rest is the remaining section of name.
!   </OUT>
!
!        arguments
!
character(len=*), intent(in)  :: name
character(len=*), intent(out) :: head
character(len=*), intent(out) :: rest

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=9),  parameter :: sub_name     = 'find_head'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer        :: i
!
!        Check for the first occurrence of the list separator in name
!
i = index(name, list_sep)
!
!        Check for additional consecutive list separators and return
!        those also
!
do while (i .le. len(name))  !{
  if (name(i+1:i+1) .eq. list_sep) then  !{
    i = i + 1
  else  !}{
    exit
  endif  !}
enddo  !}

if (i .eq. 0) then  !{
!
!        no list separators in the path, so return an empty head and
!        name as the rest
!
  head = ' '
  rest = name
elseif (i .eq. len(name)) then  !}{
!
!        The last character in name is a list separator, so return name
!        as head and an empty rest
!
  head = name
  rest = ' '
else  !}{
!
!        Found a list separator, so return the part up to the list
!        separator in head, and the remainder in rest
!
  head = name(1:i)
  rest = name(i+1:)
endif  !}

end subroutine find_head  !}
! </SUBROUTINE> NAME="find_head"
!</PRIVATE>

!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="find_list">
!
! <OVERVIEW>
!    Find and return a pointer to the specified list, relative to
!    relative_p. Return a null pointer on error.
! </OVERVIEW>
! <DESCRIPTION>
!    This function, when supplied a pointer to a field and a name of a second
!    field relative to that pointer, will find a list and return the pointer to 
!    the second field. If create is .true. and the second field does not exist,
!    it will be created.
! </DESCRIPTION>
!   <TEMPLATE>
!     list_p => find_list(path, relative_p, create)
!   </TEMPLATE>
!
function find_list(path, relative_p, create)                    &
        result (list_p)  !{
!
!   <OUT NAME="list_p" TYPE="type(field_def), pointer">
!     A pointer to the list to be returned.
!   </OUT>
!   <IN NAME="path" TYPE="character(len=*)">
!     A path to the list of interest.
!   </IN>
!   <IN NAME="list_p" TYPE="type(field_def), pointer">
!     A pointer to the list to which "path" is relative to.
!   </IN>
!   <IN NAME="create" TYPE="logical">
!     If the list does not exist, having create = .true. will create it.
!   </IN>
!
!        Function definition
!
type (field_def), pointer        :: list_p
!
!        arguments
!
character(len=*), intent(in)     :: path
type (field_def), pointer        :: relative_p
logical,          intent(in)     :: create

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=9),  parameter :: sub_name     = 'find_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len)  :: working_path
character(len=fm_path_name_len)  :: rest
character(len=fm_field_name_len) :: this_list
integer                          :: i, out_unit
type (field_def), pointer, save  :: working_path_p 
type (field_def), pointer, save  :: this_list_p 


out_unit = stdout()
nullify(list_p)
!
!        If the path is empty, then return the relative list
!
if (path .eq. ' ') then  !{

  list_p => relative_p

else  !}{
!
!        If a fully qualified path is given (i.e., starts with the
!        list separator) then do everything relative to root,
!        otherwise, do everything relative to relative list.
!
  if (path(1:1) .eq. list_sep) then  !{
    working_path_p => root_p
    working_path = path(2:)
  else  !}{
    working_path_p => relative_p
    working_path = path
  endif  !}
!
!        Loop over each field in the path
!
  do while (working_path .ne. ' ')  !{
!
!        Get the first list in the working path
!
    call find_head(working_path, this_list, rest)
!
!        If the first list is empty, then the 'rest' should hold the
!        final field in the path
!
    if (this_list .eq. ' ') then  !{
      this_list = rest
      rest = ' '
    endif  !}
!
!        Strip off trailing list separators
!
    i = len_trim(this_list)
    do while (i .gt. 0 .and. this_list(i:i) .eq. list_sep)  !{
      this_list(i:i) = ' '
      i = i - 1
    enddo  !}
!
!        Find a pointer to this field in the working list
!
    this_list_p => find_field(this_list, working_path_p)

    if (.not. associated(this_list_p)) then  !{
      if (create) then  !{
!
!        Create the list if so requested
!
        this_list_p => make_list(working_path_p, this_list)
        if (.not. associated(this_list_p)) then  !{
          if (verb .gt. verb_level_warn) then  !{
            write (out_unit,*) trim(warn_header), 'List "',       &
                 trim(this_list), '" could not be created in ',   &
                 trim(path)
          endif  !}
          nullify(list_p)
          return
        endif  !}
      else  !}{
!
!        Otherwise, return an error
!

        if (verb .gt. verb_level_note) then  !{
          write (out_unit,*) trim(note_header), 'List "',         &
               trim(this_list), '" does not exist in ', trim(path)
        endif  !}
        nullify(list_p)
        return
      endif  !}
    endif  !}
!
!        Make sure that the field found is a list, and if so, proceed to
!        the next field in the path, otherwise, return an error
!
    if (this_list_p%field_type .eq. list_type) then  !{
      working_path_p => this_list_p
      working_path = rest
    else  !}{
      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header), '"',                &
             trim(this_list), '" is not a list in ', trim(path)
      endif  !}
      nullify(list_p)
      return
    endif  !}
  enddo  !}
  list_p => working_path_p
endif  !}

end function find_list  !}
! </FUNCTION> NAME="find_list"
!</PRIVATE>

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_change_list">
!
! <OVERVIEW>
!    Change the current list. Return true on success,
!    false otherwise
! </OVERVIEW>
! <DESCRIPTION>
!    This function changes the currect list to correspond to the list named name.
!    If the first character of name is the list separator (/) then the list will 
!    search for "name" starting from the root of the field tree. Otherwise it 
!    will search for name starting from the current list.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_change_list(name)
!   </TEMPLATE>
!
function fm_change_list(name)                                        &
        result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list that the user wishes to change to.
!   </IN>
!
!        Function definition
!
logical        :: success
!
!        arguments
!
character(len=*), intent(in)  :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=14), parameter :: sub_name     = 'fm_change_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_p 
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Find the list if path is not empty
!
temp_p => find_list(name, current_list_p, .false.)

if (associated(temp_p)) then  !{
  current_list_p => temp_p
  success = .true.
else  !}{
  success = .false.
endif  !}

end function fm_change_list  !}
! </FUNCTION> NAME="fm_change_list"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_change_root">
!
! <OVERVIEW>
!    Change the root list
! </OVERVIEW>
! <DESCRIPTION>
!    This function changes the root of the field tree to correspond to the 
!    field named name. An example of a use of this would be if code is 
!    interested in a subset of fields with a common base. This common base 
!    could be set using fm_change_root and fields could be referenced using 
!    this root. 
!    
!    This function should be used in conjunction with fm_return_root.
!    
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_change_root(name)
!   </TEMPLATE>
!
function  fm_change_root(name)                                        &
          result (success)  !{
!
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of the field which the user wishes to become the root.
!   </IN>
!
!        Function definition
!
logical        :: success
!
!        arguments
!
character(len=*), intent(in)  :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=14), parameter :: sub_name     = 'fm_change_root'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_list_p 
integer :: out_unit
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
out_unit = stdout()
!
!        Must supply a field field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  success = .false.
  return
endif  !}
!
!        Get a pointer to the list
!
temp_list_p => find_list(name, current_list_p, .false.)

if (associated(temp_list_p)) then  !{
!
!        restore the saved root values if we've already changed root
!
  if (save_root_name .ne. ' ') then  !{
    root_p%name = save_root_name
    root_p%parent => save_root_parent_p
  endif  !}
!
!        set the pointer for the new root field
!
  root_p => temp_list_p
!
!        save the new root field's name and parent
!
  save_root_name = root_p%name
  save_root_parent_p => root_p%parent
!
!        set the new root name and parent fields to appropriate values
!
  root_p%name = ' '
  nullify(root_p%parent)
!
!        set the current list to the new root as it likely is not
!        going to be meaningful anymore
!
  current_list_p => root_p
  success = .true.
else  !}{
!
!        Couldn't find the list
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                      &
         'Could not find list ', trim(name)
  endif  !}
  success = .false.
endif  !}

end function  fm_change_root  !}
! </FUNCTION> NAME="fm_change_root"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_dump_list">
!
! <OVERVIEW>
!    A function to list properties associated with a field.
! </OVERVIEW>
! <DESCRIPTION>
!    This function writes the contents of the field named "name" to stdout.
!    If recursive is present and .true., then this function writes out the 
!    contents of any subfields associated with the field named "name".
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_dump_list(name, recursive = .true.) 
!   </TEMPLATE>
!
function  fm_dump_list(name, recursive)                        &
          result (success)  !{
!
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of the field for which output is requested.
!   </IN>
!   <IN NAME="recursive" TYPE="logical, optional">
!     If present and .true., then a recursive listing of fields will be
!     performed.
!   </IN>
!
!        Function definition
!
logical        :: success
!
!        arguments
!
character(len=*), intent(in)           :: name
logical,          intent(in), optional :: recursive

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=12), parameter :: sub_name     = 'fm_dump_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
logical                         :: recursive_t
type (field_def), pointer, save :: temp_list_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Check whether to do things recursively
!
if (present(recursive)) then  !{
  recursive_t = recursive
else  !}{
  recursive_t = .false.
endif  !}
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}

if (name .eq. ' ') then  !{
!
!        If list is empty, then dump the current list
!
  temp_list_p => current_list_p
  success = .true.
else  !}{
!
!        Get a pointer to the list
!
  temp_list_p => find_list(name, current_list_p, .false.)
  if (associated(temp_list_p)) then  !{
    success = .true.
  else  !}{
!
!        Error following the path
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                        &
           'Could not follow path for ', trim(name)
    endif  !}
    success = .false.
  endif  !}
endif  !}
!
!        Dump the list
!
if (success) then  !{
  success = dump_list(temp_list_p, recursive_t, 0)
endif  !}

end function  fm_dump_list  !}
! </FUNCTION> NAME="fm_dump_list"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_exists">
!
! <OVERVIEW>
!   A function to test whether a named field exists.
! </OVERVIEW>
! <DESCRIPTION>
!   This function determines is a field exists, relative to the current list,
!   and returns true if the list exists, false otherwise.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_exists(name)
!   </TEMPLATE>
!
function fm_exists(name)                                                &
        result (success)  !{
!
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of the field that is being queried.
!   </IN>
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!
!        Function definition
!
logical        :: success
!
!        arguments
!
character(len=*), intent(in) :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=9),  parameter :: sub_name     = 'fm_exists'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: dummy_p 
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Determine whether the field exists
!
dummy_p => get_field(name, current_list_p)
success = associated(dummy_p)

end function fm_exists  !}
! </FUNCTION> NAME="fm_exists"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_get_index">
!
! <OVERVIEW>
!    A function to return the index of a named field.
! </OVERVIEW>
! <DESCRIPTION>
!    Returns the index for name, returns the parameter NO_FIELD if it does not
!    exist. If the first character of the named field is the list peparator, 
!    then the named field will be relative to the root of the field tree. 
!    Otherwise the named field will be relative to the current list. 
! </DESCRIPTION>
!   <TEMPLATE>
!     index = fm_get_index(name)
!   </TEMPLATE>
!
function  fm_get_index(name)                        &
          result (index)  !{
!   <OUT NAME="index" TYPE="index">
!     The index of the named field if it exists. 
!     Otherwise the parameter NO_FIELD.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a field that the user wishes to get an index for.
!   </IN>
!
!        Function definition
!
integer        :: index
!
!        arguments
!
character(len=*), intent(in) :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=12), parameter :: sub_name     = 'fm_get_index'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  index = NO_FIELD
  return
endif  !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)
if (associated(temp_field_p)) then  !{
!
!        Set the index
!
  index = temp_field_p%index
else  !}{
!
!        Error following the path
!
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Could not follow path for ', trim(name)
  endif  !}
  index = NO_FIELD
endif  !}

end function  fm_get_index  !}
! </FUNCTION> NAME="fm_get_index"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_get_current_list">
!
! <OVERVIEW>
!    A function to return the full path of the current list.
! </OVERVIEW>
! <DESCRIPTION>
!    This function returns the full path for the current list. A blank 
!    path indicates an error condition has occurred.
! </DESCRIPTION>
!   <TEMPLATE>
!     path = fm_get_current_list()
!   </TEMPLATE>
!
function  fm_get_current_list()                                        &
          result (path)  !{
!
!   <OUT NAME="path" TYPE="character(len=fm_path_name_len)">
!     The path corresponding to the current list.
!   </OUT>
!
!        Function definition
!
character(len=fm_path_name_len) :: path
!
!        arguments
!

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=19), parameter :: sub_name     = 'fm_get_current_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_list_p 
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Set a pointer to the current list and proceed
!        up the tree, filling in the name as we go
!
temp_list_p => current_list_p
path = ' '

do while (associated(temp_list_p))  !{
!
!        Check whether we are at the root field--it is the
!        only field with a blank name
!
  if (temp_list_p%name .eq. ' ') then  !{
    exit
  endif  !}
!
!        Append the name to the path
!
  path = list_sep // trim(temp_list_p%name) // path
!
!        Point to the next field
!
  temp_list_p => temp_list_p%parent
enddo  !}

if (.not. associated(temp_list_p)) then  !{
!
!        The pointer is not associated, indicating an error has
!        occurred, so set the path accordingly
!
  path = ' '
elseif (path .eq. ' ') then  !}{
!
!        If path is empty, then the current list must be root,
!        so set path accordingly
!
  path = list_sep
endif  !}

end function  fm_get_current_list  !}
! </FUNCTION> NAME="fm_get_current_list"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_get_length">
!
! <OVERVIEW>
!    A function to return how many elements are contained within the named 
!    list or entry.
! </OVERVIEW>
! <DESCRIPTION>
!    This function returns the list or entry length for the named list or entry.
!    If the named field or entry does not exist, a value of 0 is returned.
! </DESCRIPTION>
!   <TEMPLATE>
!     length = fm_get_length(name)
!   </TEMPLATE>
!
function  fm_get_length(name)                        &
          result (length)  !{
!
!   <OUT NAME="length" TYPE="integer">
!     The number of elements that the field name has.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list or entry that the user wishes to get the length of.
!   </IN>
!
!        Function definition
!
integer                      :: length
!
!        arguments
!
character(len=*), intent(in) :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=13), parameter :: sub_name     = 'fm_get_length'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  length = 0
  return
endif  !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then  !{
!
!        Set the field length
!
  if (temp_field_p%field_type .eq. list_type) then !{
    length = temp_field_p%length
  else !}{
    length = temp_field_p%max_index
  endif !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                            &
         'Could not follow path for ', trim(name)
  endif  !}
  length = 0
endif  !}

end function  fm_get_length  !}
! </FUNCTION> NAME="fm_get_length"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_get_type">
!
! <OVERVIEW>
!    A function to return the type of the named field.
! </OVERVIEW>
! <DESCRIPTION>
!    This function returns the type of the field for name.
!    This indicates whether the named field is a "list" (has children fields),
!    or has values of type "integer", "real", "logical" or "string".
!    If it does not exist it returns a blank string. 
! </DESCRIPTION>
!   <TEMPLATE>
!     name_field_type = fm_get_type(name)
!   </TEMPLATE>
!
function  fm_get_type(name)                        &
          result (name_field_type)  !{
!   <OUT NAME="name_field_type" TYPE="character(len=8)">
!     A string containing the type of the named field.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a field that the user wishes to find the type of.
!   </IN>
!
!        Function definition
!
character(len=8)             :: name_field_type
!
!        arguments
!
character(len=*), intent(in) :: name

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=11), parameter :: sub_name     = 'fm_get_type'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  name_field_type = ' '
  return
endif  !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then  !{
!
!        Set the field type
!
  name_field_type = field_type_name(temp_field_p%field_type)
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                            &
         'Could not follow path for ', trim(name)
  endif  !}
  name_field_type = ' '
endif  !}

end function  fm_get_type  !}
! </FUNCTION> NAME="fm_get_type"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_get_value">
!
! <OVERVIEW>
!    An overloaded function to find and extract a value for a named field.
! </OVERVIEW>
! <DESCRIPTION>
!    Find and extract the value for name. The value may be of type real, 
!    integer, logical or character. If a single value from an array  of values 
!    is required, an optional index can be supplied.
!    Return true for success and false for failure
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_get_value(name, value, index)
!   </TEMPLATE>
!
function  fm_get_value_integer(name, value, index)                 &
          result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a field that the user wishes to get a value for.
!   </IN>
!   <OUT NAME="value" TYPE="integer, real, logical or character">
!     The value associated with the named field.
!   </OUT>
!   <IN NAME="index" TYPE="integer, optional">
!     An optional index to retrieve a single value from an array.
!   </IN>
!
!        Function definition
!
logical                                :: success
!
!        arguments
!
character(len=*), intent(in)           :: name
integer,          intent(out)          :: value
integer,          intent(in), optional :: index

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=20), parameter :: sub_name     = 'fm_get_value_integer'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: index_t
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  value = 0
  success = .false.
  return
endif  !}
!
!        Set index to retrieve
!
if (present(index)) then  !{
  index_t = index
else !}{
  index_t = 1
endif !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then  !{
!
!        check that the field is the correct type
!
  if (temp_field_p%field_type .eq. integer_type) then  !{
    if (index_t .lt. 1) then  !{
!
!        Index is not positive
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Optional index for ', trim(name),                 &
             ' not positive: ', index_t
      endif  !}
      value = 0
      success = .false.
    elseif (index_t .gt. temp_field_p%max_index) then  !}{
!
!        Index is too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                        &
             'Optional index for ', trim(name),                      &
             ' too large: ', index_t, ' > ', temp_field_p%max_index
      endif  !}
      value = 0
      success = .false.
    else  !}{
!
!        extract the value
!
      value = temp_field_p%i_value(index_t)
      success = .true.
    endif !}
  else  !}{
!
!        Field not corrcet type
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Field not type integer ', trim(name)
    endif  !}
    value = 0
    success = .false.
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                            &
         'Could not follow path for ', trim(name)
  endif  !}
  value = 0
  success = .false.
endif  !}

end function  fm_get_value_integer  !}

!#######################################################################
!#######################################################################

function  fm_get_value_logical(name, value, index)                 &
          result (success)  !{
!
!        Function definition
!
logical                                :: success
!
!        arguments
!
character(len=*), intent(in)           :: name
logical,          intent(out)          :: value
integer,          intent(in), optional :: index

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=20), parameter :: sub_name     = 'fm_get_value_logical'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: index_t
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  value = .false.
  success = .false.
  return
endif  !}
!
!        Set index to retrieve
!
if (present(index)) then  !{
  index_t = index
else  !}{
  index_t = 1
endif  !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then  !{
!
!        check that the field is the correct type
!
  if (temp_field_p%field_type .eq. logical_type) then  !{

    if (index_t .lt. 1) then  !{
!
!        Index is not positive
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Optional index for ', trim(name),                 &
             ' not positive: ', index_t
      endif  !}
      value = .false.
      success = .false.

    elseif (index_t .gt. temp_field_p%max_index) then  !}{
!
!        Index is too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                        &
             'Optional index for ', trim(name),                      &
             ' too large: ', index_t, ' > ', temp_field_p%max_index
      endif  !}
      value = .false.
      success = .false.

    else  !}{
!
!        extract the value
!
      value = temp_field_p%l_value(index_t)
      success = .true.
    endif !}
  else  !}{
!
!        Field not correct type
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Field not type logical ', trim(name)
    endif  !}
    value = .false.
    success = .false.
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                            &
         'Could not follow path for ', trim(name)
  endif  !}
  value = .false.
  success = .false.
endif  !}

end function  fm_get_value_logical  !}

!#######################################################################
!#######################################################################

function  fm_get_value_real(name, value, index)                 &
          result (success)  !{
!
!        Function definition
!
logical                                :: success
!
!        arguments
!
character(len=*), intent(in)           :: name
real,             intent(out)          :: value
integer,          intent(in), optional :: index

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=17), parameter :: sub_name     = 'fm_get_value_real'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: index_t
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  value = 0.0
  success = .false.
  return
endif  !}
!
!        Set index to retrieve
!
if (present(index)) then  !{
  index_t = index
else  !}{
  index_t = 1
endif  !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then  !{
!
!        check that the field is the correct type
!
  if (temp_field_p%field_type .eq. real_type) then  !{

    if (index_t .lt. 1) then  !{

!
!        Index is not positive
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                        &
             'Optional index for ', trim(name),                      &
             ' not positive: ', index_t
      endif  !}
      value = 0.0
      success = .false.

    elseif (index_t .gt. temp_field_p%max_index) then  !}{

!
!        Index is too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                        &
             'Optional index for ', trim(name),                      &
             ' too large: ', index_t, ' > ', temp_field_p%max_index
      endif  !}
      value = 0.0
      success = .false.

    else  !}{

!
!        extract the value
!
      value = temp_field_p%r_value(index_t)
      success = .true.
    endif !}
  else  !}{
!
!        Field not correct type
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Field not type real ', trim(name)
    endif  !}
    value = 0.0
    success = .false.
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                            &
         'Could not follow path for ', trim(name)
  endif  !}
  value = 0.0
  success = .false.
endif  !}

end function  fm_get_value_real  !}

!#######################################################################
!#######################################################################

function  fm_get_value_string(name, value, index)                 &
          result (success)  !{
!
!        Function definition
!
logical                                :: success
!
!        arguments
!
character(len=*), intent(in)           :: name
character(len=*), intent(out)          :: value
integer,          intent(in), optional :: index

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=19), parameter :: sub_name     = 'fm_get_value_string'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: index_t
type (field_def), pointer, save :: temp_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  value = ''
  success = .false.
  return
endif  !}
!
!        Set index to retrieve
!
if (present(index)) then  !{
  index_t = index
else  !}{
  index_t = 1
endif  !}
!
!        Get a pointer to the field
!
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then  !{
!
!        check that the field is the correct type
!
  if (temp_field_p%field_type .eq. string_type) then  !{
    if (index_t .lt. 1) then  !{
!
!        Index is not positive
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                        &
             'Optional index for ', trim(name),                      &
             ' not positive: ', index_t
      endif  !}
      value = ''
      success = .false.

    elseif (index_t .gt. temp_field_p%max_index) then  !}{
!
!        Index is too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                        &
             'Optional index for ', trim(name),                      &
             ' too large: ', index_t, ' > ', temp_field_p%max_index
      endif  !}
      value = ''
      success = .false.
    else  !}{
!
!        extract the value
!
      value = temp_field_p%s_value(index_t)
      !if (trim(value) == '') then
        !success = .false.
      !else
        success = .true.
      !endif
    endif !}
  else  !}{
!
!        Field not correct type
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Field not type string ', trim(name)
    endif  !}
    value = ''
    success = .false.
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                            &
         'Could not follow path for ', trim(name)
  endif  !}
  value = ''
  success = .false.
endif  !}

end function  fm_get_value_string  !}
! </FUNCTION> NAME="fm_get_value"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_intersection">
!
! <OVERVIEW>
!    A function to find the common names of the sub-fields in a list 
!    of fields.
! </OVERVIEW>
! <DESCRIPTION>
!    Return a pointer to an fm_array_list of the intersection
!    of an array of lists, ignoring the contents of the values,
!    but just returning the names.
!    Return false on the end of the intersection.
! </DESCRIPTION>
!   <TEMPLATE>
!     return_p => fm_intersection(lists,dim)
!   </TEMPLATE>
!
function fm_intersection(lists, dim)                        &
        result (return_p)  !{
!   <OUT NAME="return_p" TYPE="type (fm_array_list_def), pointer">
!     A pointer to a list of names that are common to the fields provided in 
!     lists.
!   </OUT>
!   <IN NAME="dim" TYPE="dim">
!     The dimension of lists.
!   </IN>
!   <IN NAME="lists" TYPE="character(len=*)" DIM="(dim)">
!     A list of fields that the user wishes to find the common fields of.
!   </IN>
!
!        Function definition
!
type (fm_array_list_def), pointer  :: return_p
!
!        arguments
!
integer,          intent(in)       :: dim
character(len=*), intent(in)       :: lists(dim)

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=15), parameter :: sub_name     = 'fm_intersection'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character (len=fm_field_name_len)  :: name
character (len=fm_field_name_len),                          &
        dimension(:), allocatable  :: names
character (len=fm_type_name_len)   :: field_type
integer                            :: count
integer                            :: error
integer                            :: index
integer                            :: n, ier
integer                            :: shortest
logical                            :: found
type (field_def), pointer, save    :: temp_p 
integer                            :: out_unit

out_unit = stdout()

nullify(return_p)
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        return error if dimension if bad
!
if (dim .le. 0) then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Non-positive dimension: ', dim
  endif  !}
  nullify(return_p)
  return
endif  !}
!
!        make sure that the lists exist, and find the shortest list
!        and its length
!
count = -1
shortest = 0
do n = 1, dim  !{
  temp_p => find_list(lists(n), current_list_p, .false.)
  if (associated(temp_p)) then  !{
    if (count .eq. -1) then  !{
      count = temp_p%length
      shortest = n
    else  !}{
      if (count .gt. temp_p%length) then  !{
        count = temp_p%length
        shortest = n
      endif  !}
    endif  !}
  else  !}{
    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
                         'List does not exist: "', trim(lists(n)), '"'
    endif  !}
    nullify(return_p)
    return
  endif  !}
enddo  !} n
!
!        allocate return pointer
!
allocate( return_p, stat = error)
if (error .ne. 0) then !{
  write (out_unit,*) trim(error_header), 'Error ', error          &
                 , ' allocating memory for return_p '
  nullify(return_p)
  return
endif  !}
if ( associated(return_p%names)) deallocate(return_p%names)
!
!        return if any list is empty
!
if (count .eq. 0) then  !{
  return_p%length = 0
  return
endif  !}
!
!        If there is only one list, then return its names
!
if (dim .eq. 1) then  !{
!
!        allocate space for names in return pointer
!
  allocate( return_p%names(count), stat = error)
  if (error .ne. 0) then !{
    write (out_unit,*) trim(error_header), 'Error ', error        &
                   , ' allocating memory for names in return_p '
    nullify(return_p)
    return
  endif  !}
  count = 0
  do while (fm_loop_over_list(lists(1), name, field_type, index))  !{
    count = count + 1
    return_p%names(count) = name
  enddo  !}
  return
endif  !}
!
!        allocate space for names
!
allocate( names(count), stat = error)
if (error .ne. 0) then !{
  write (out_unit,*) trim(error_header), 'Error ', error          &
                 , ' allocating memory for names '
  nullify(return_p)
  return
endif  !}
!
!        Loop over the shortest list, checking whether its names
!        occur in all of the other lists. If so, then save the name
!
count = 0
do while (fm_loop_over_list(lists(shortest), name, field_type, index))  !{
  found = .true.
  do n = 1, dim  !{
    if (n .ne. shortest) then   !{
      temp_p => find_list(trim(lists(n)) // list_sep // name,        &
                          current_list_p, .false.)
      if (.not. associated(temp_p)) then  !{
        found = .false.
        exit
      endif  !}
    endif  !}
  enddo  !}
  if (found) then  !{
    count = count + 1
    names(count) = name
  endif  !}
enddo  !}
!
!        allocate space for names in return pointer
!
allocate( return_p%names(count), stat = error)
if (error .ne. 0) then !{
  write (out_unit,*) trim(error_header), 'Error ', error  &
                 , ' allocating memory for names in return_p '
  deallocate(names)
  nullify(return_p)
  return
endif  !}
!
!        copy the names to the return pointer and clean up
!
do n = 1, count  !{
  return_p%names(n) = names(n)
enddo  !} n
return_p%length = count
deallocate(names)

end function fm_intersection  !}
! </FUNCTION> NAME="fm_intersection"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_loop_over_list">
!
! <OVERVIEW>
!    A function for looping over a list.
! </OVERVIEW>
! <DESCRIPTION>
!    Loop over the list, setting the name, type and index
!    of the next field. Return false at the end of the loop.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_loop_over_list(list, name, field_type, index)
!   </TEMPLATE>
!
function  fm_loop_over_list(list, name, field_type, index)        &
          result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="list" TYPE="character(len=*)">
!     The name of a list to loop over.
!   </IN>
!   <OUT NAME="name" TYPE="character(len=*)">
!     The name of a field from list.
!   </OUT>
!   <OUT NAME="field_type" TYPE="character(len=fm_type_name_len)">
!     The type of a list entry.
!   </OUT>
!   <OUT NAME="index" TYPE="integer">
!     The index of tje field within the list.
!   </OUT>
!
!        Function definition
!
logical                                      :: success
!
!        arguments
!
character(len=*),                intent(in)  :: list
character(len=*),                intent(out) :: name
character(len=fm_type_name_len), intent(out) :: field_type
integer,                         intent(out) :: index

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=17), parameter :: sub_name     = 'fm_loop_over_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
type (field_def), pointer, save :: temp_list_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}

if (list .eq. loop_list .and. associated(loop_list_p)) then  !{
!
!        We've already started this loop, so continue on
!
  loop_list_p => loop_list_p%next
  success = set_list_stuff()
elseif (list .eq. ' ') then  !{
!
!        If list is empty, then loop over the current list
!
  loop_list = ' '
  loop_list_p => current_list_p%first_field
  success = set_list_stuff()
else  !}{
!
!        Get a pointer to the list
!
  loop_list = list
  loop_list_p => find_list(loop_list, current_list_p, .false.)
  if (associated(loop_list_p)) then  !{
    loop_list_p => loop_list_p%first_field
    success = set_list_stuff()
  else  !}{
!
!        Error following the path
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                        &
           'Could not follow path for ', trim(list)
    endif  !}
    success = .false.
  endif  !}
endif  !}

return

contains

!#######################################################################
!#######################################################################

! <FUNCTION NAME="set_list_stuff">
!
! <DESCRIPTION>
! If the the pointer matches to the right list,
! extract the field information.  Used in fm_loop_over_list
! </DESCRIPTION>
function  set_list_stuff()                                                &
          result (success)  !{
!
!        Function definition
!
  logical        :: success
!
!        arguments
!
  if (associated(loop_list_p)) then  !{
    name = loop_list_p%name
    field_type = field_type_name(loop_list_p%field_type)
    index = loop_list_p%index
    success = .true.
  else  !}{
    name = ' '
    field_type = ' '
    index = 0
    success = .false.
    loop_list = ' '
  endif  !}

end function  set_list_stuff  !}
! </FUNCTION> NAME="set_list_stuff"

end function  fm_loop_over_list  !}
! </FUNCTION> NAME="fm_loop_over_list"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_new_list">
!
! <OVERVIEW>
!    A function to create a new list.
! </OVERVIEW>
! <DESCRIPTION>
!    Allocate and initialize a new list and return the index of the list. 
!    If an error occurs return the parameter NO_FIELD.
! </DESCRIPTION>
!   <TEMPLATE>
!     index = fm_new_list(name, create, keep)
!   </TEMPLATE>
!
function  fm_new_list(name, create, keep)                        &
          result (index)  !{
!   <OUT NAME="index" TYPE="integer">
!     The index of the newly created list.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list that the user wishes to create.
!   </IN>
!   <IN NAME="create" TYPE="logical, optional">
!     If present and .true., create the list if it does not exist.
!   </IN>
!   <IN NAME="keep" TYPE="logical, optional">
!     If present and .true., make this list the current list.
!   </IN>
!
!        Function definition
!
integer                                :: index
!
!        arguments
!
character(len=*), intent(in)           :: name
logical,          intent(in), optional :: create
logical,          intent(in), optional :: keep

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=11), parameter :: sub_name     = 'fm_new_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
logical                          :: create_t
logical                          :: keep_t
character(len=fm_path_name_len)  :: path
character(len=fm_field_name_len) :: base
type (field_def), pointer, save  :: temp_list_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field list name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a list name'
  endif  !}
  index = NO_FIELD
  return
endif  !}
!
!        Check for optional arguments
!
if (present(create)) then  !{
  create_t = create
else  !}{
  create_t = .false.
endif  !}

if (present(keep)) then  !{
  keep_t = keep
else  !}{
  keep_t = .false.
endif  !}
!
!        Get a pointer to the parent list
!
call find_base(name, path, base)

temp_list_p => find_list(path, current_list_p, create_t)

if (associated(temp_list_p)) then  !{
!
!        Create the list
!
  temp_list_p => make_list(temp_list_p, base)
  if (associated(temp_list_p)) then  !{
!
!        Make this list the current list, if requested
!
    if (keep_t) then  !{
      current_list_p => temp_list_p
    endif  !}
    index = temp_list_p%index
  else  !}{
!
!        Error in making the list
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                        &
           'Could not create list ', trim(name)
    endif  !}
    index = NO_FIELD

  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                  &
         'Could not follow path for ', trim(name)
  endif  !}
  index = NO_FIELD

endif  !}

end function  fm_new_list  !}
! </FUNCTION> NAME="fm_new_list"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_new_value">
!
! <OVERVIEW>
!    An overloaded function to assign a value to a field.
! </OVERVIEW>
! <DESCRIPTION>
!    Allocate and initialize a new value and return the index.
!    If an error condition occurs the parameter NO_FIELD is returned.
!
!    If the type of the field is changing (e.g. real values being transformed to
!    integers), then any previous values for the field are removed and replaced 
!    by the value passed in the present call to this function.
!     
!    If append is present and .true., then index cannot be greater than 0 if 
!    it is present.
! </DESCRIPTION>
!   <TEMPLATE>
!     field_index = fm_new_value(name, value, [create], [index], [append])
!   </TEMPLATE>
!
function  fm_new_value_integer(name, value, create, index, append)     &
          result (field_index)  !{
!   <OUT NAME="field_index" TYPE="integer">
!     An index for the named field.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a field that the user wishes to create a value for.
!   </IN>
!   <IN NAME="value" TYPE="integer, real, logical, or character(len=*)">
!     The value that the user wishes to apply to the named field.
!   </IN>
!   <IN NAME="create" TYPE="logical, optional">
!     If present and .true., then a value for this field will be created.
!   </IN>
!   <IN NAME="index" TYPE="integer, optional">
!     The index to an array of values that the user wishes to apply a new value.
!   </IN>
!   <IN NAME="append" TYPE="logical, optional">
!     If present and .true., then append the value to an array of the present 
!     values. If present and .true., then index cannot be greater than 0.
!   </IN>
!
!        Function definition
!
integer                                :: field_index
!
!        arguments
!
character(len=*), intent(in)           :: name
integer,          intent(in)           :: value
logical,          intent(in), optional :: create
integer,          intent(in), optional :: index
logical,          intent(in), optional :: append

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=20), parameter :: sub_name     = 'fm_new_value_integer'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
logical                          :: create_t
integer                          :: i, ier
integer                          :: index_t
integer, pointer, dimension(:)   :: temp_i_value
character(len=fm_path_name_len)  :: path
character(len=fm_field_name_len) :: base
type (field_def), pointer, save  :: temp_list_p 
type (field_def), pointer, save  :: temp_field_p 
integer                          :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  field_index = NO_FIELD
  return
endif  !}
!
!        Check for optional arguments
!
if (present(create)) then  !{
  create_t = create
else  !}{
  create_t = .false.
endif  !}
!
!        Check that append is not true and index non-positive
!

if (present(index) .and. present(append)) then  !{
  if (append .and. index .gt. 0) then  !{
    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Index and Append both set for ', trim(name)
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
endif  !}
!
!        Set index to define
!
if (present(index)) then  !{
  index_t = index
  if (index_t .lt. 0) then  !{
!
!        Index is negative
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Optional index for ', trim(name),                   &
           ' negative: ', index_t
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
else  !}{
  index_t = 1
endif  !}
!
!        Get a pointer to the parent list
!
call find_base(name, path, base)
temp_list_p => find_list(path, current_list_p, create_t)

if (associated(temp_list_p)) then  !{
  temp_field_p => find_field(base, temp_list_p)
  if (.not. associated(temp_field_p)) then  !{
!
!        Create the field if it doesn't exist
!
    temp_field_p => create_field(temp_list_p, base)
  endif  !}
  if (associated(temp_field_p)) then  !{
!
!        Check if the field_type is the same as previously
!        If not then reset max_index to 0
!
    if (temp_field_p%field_type /= integer_type ) then
        temp_field_p%max_index = 0
      if (temp_field_p%field_type /= null_type ) then  !{
        if (verb .gt. verb_level_warn) then  !{
          write (out_unit,*) trim(warn_header),                   &
               'Changing type of ', trim(name), ' from ',         &
               trim(field_type_name(temp_field_p%field_type)),    &
               ' to ', trim(field_type_name(integer_type))
        endif  !}
      endif  !}
    endif
!
!        Assign the type
!
    temp_field_p%field_type = integer_type
!
!        Set the index if appending
!

    if (present(append)) then  !{
      if (append) then  !{
        index_t = temp_field_p%max_index + 1
      endif  !}
    endif  !}

    if (index_t .gt. temp_field_p%max_index + 1) then  !{

!
!        Index too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Index too large for ', trim(name), ': ', index_t
      endif  !}
      field_index = NO_FIELD
      return

    elseif (index_t .eq. 0 .and.                                &
            temp_field_p%max_index .gt. 0) then  !}{
!
!        Can't set non-null field to null
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Trying to nullify a non-null field: ',            &
             trim(name)
      endif  !}
      field_index = NO_FIELD
      return

    elseif (.not. associated(temp_field_p%i_value) .and.        &
            index_t .gt. 0) then  !}{
!
!        Array undefined, so allocate the array
!
      allocate(temp_field_p%i_value(1))
      temp_field_p%max_index = 1
      temp_field_p%array_dim = 1
    elseif (index_t .gt. temp_field_p%array_dim) then  !}{
!
!        Array is too small, so allocate new array and copy over
!        old values
!
      temp_field_p%array_dim = temp_field_p%array_dim + array_increment
      allocate (temp_i_value(temp_field_p%array_dim))
      do i = 1, temp_field_p%max_index  !{
        temp_i_value(i) = temp_field_p%i_value(i)
      enddo  !} i
      if (associated (temp_field_p%i_value)) deallocate(temp_field_p%i_value)
      temp_field_p%i_value => temp_i_value
      temp_field_p%max_index = index_t
    endif  !}
!
!        Assign the value and set the field_index for return
!        for non-null fields (index_t > 0)
!
    if (index_t .gt. 0) then  !{
      temp_field_p%i_value(index_t) = value
      if (index_t .gt. temp_field_p%max_index) then  !{
        temp_field_p%max_index = index_t
      endif  !}
    endif  !}
    field_index = temp_field_p%index

  else  !}{
!
!        Error in making the field
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Could not create integer value field ',             &
           trim(name)
    endif  !}
    field_index = NO_FIELD
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                       &
         'Could not follow path for ',                          &
         trim(name)
  endif  !}
  field_index = NO_FIELD
endif  !}

end function  fm_new_value_integer  !}

!#######################################################################
!#######################################################################

function  fm_new_value_logical(name, value, create, index, append) &
          result (field_index)  !{
!
!        Function definition
!
integer                                :: field_index
!
!        arguments
!
character(len=*), intent(in)           :: name
logical,          intent(in)           :: value
logical,          intent(in), optional :: create
integer,          intent(in), optional :: index 
logical,          intent(in), optional :: append

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=20), parameter :: sub_name     = 'fm_new_value_logical'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len)      :: path
character(len=fm_field_name_len)     :: base
integer                              :: i, ier
integer                              :: index_t
logical                              :: create_t
logical, dimension(:), pointer       :: temp_l_value 
type (field_def),      pointer, save :: temp_list_p 
type (field_def),      pointer, save :: temp_field_p 
integer                              :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  field_index = NO_FIELD
  return
endif  !}
!
!        Check for optional arguments
!
if (present(create)) then  !{
  create_t = create
else  !}{
  create_t = .false.
endif  !}
!
!        Check that append is not true and index greater than 0
!
if (present(index) .and. present(append)) then  !{
  if (append .and. index .gt. 0) then  !{
    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Index and Append both set for ', trim(name)
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
endif  !}
!
!        Set index to define
!

if (present(index)) then  !{
  index_t = index
  if (index_t .lt. 0) then  !{
!
!        Index is negative
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Optional index for ', trim(name),                   &
           ' negative: ', index_t
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
else  !}{
  index_t = 1
endif !}
!
!        Get a pointer to the parent list
!
call find_base(name, path, base)
temp_list_p => find_list(path, current_list_p, create_t)

if (associated(temp_list_p)) then  !{
  temp_field_p => find_field(base, temp_list_p)
  if (.not. associated(temp_field_p)) then  !{
!
!        Create the field if it doesn't exist
!
    temp_field_p => create_field(temp_list_p, base)
  endif  !}
  if (associated(temp_field_p)) then  !{
!
!        Check if the field_type is the same as previously
!        If not then reset max_index to 0
!
    if (temp_field_p%field_type /= logical_type ) then
        temp_field_p%max_index = 0
      if (temp_field_p%field_type /= null_type ) then  !{
        if (verb .gt. verb_level_warn) then  !{
          write (out_unit,*) trim(warn_header),                   &
               'Changing type of ', trim(name), ' from ',         &
               trim(field_type_name(temp_field_p%field_type)),    &
               ' to ', trim(field_type_name(logical_type))
        endif  !}
      endif  !}
    endif
!
!        Assign the type
!
    temp_field_p%field_type = logical_type
!
!        Set the index if appending
!

    if (present(append)) then  !{
      if (append) then  !{
        index_t = temp_field_p%max_index + 1
      endif  !}
    endif  !}

    if (index_t .gt. temp_field_p%max_index + 1) then  !{

!
!        Index too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Index too large for ', trim(name), ': ', index_t
      endif  !}
      field_index = NO_FIELD
      return

    elseif (index_t .eq. 0 .and.                                &
            temp_field_p%max_index .gt. 0) then  !}{

!
!        Can't set non-null field to null
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Trying to nullify a non-null field: ', trim(name)
      endif  !}
      field_index = NO_FIELD
      return

    elseif (.not. associated(temp_field_p%l_value) .and.        &
            index_t .gt. 0) then  !}{

!
!        Array undefined, so allocate the array
!

      allocate(temp_field_p%l_value(1))
      temp_field_p%max_index = 1
      temp_field_p%array_dim = 1

    elseif (index_t .gt. temp_field_p%array_dim) then  !}{

!
!        Array is too small, so allocate new array and copy over
!        old values
!
      temp_field_p%array_dim = temp_field_p%array_dim + array_increment
      allocate (temp_l_value(temp_field_p%array_dim))
      do i = 1, temp_field_p%max_index  !{
        temp_l_value(i) = temp_field_p%l_value(i)
      enddo  !} i
      if (associated(temp_field_p%l_value)) deallocate(temp_field_p%l_value)
      temp_field_p%l_value => temp_l_value
      temp_field_p%max_index = index_t

    endif  !}

!
!        Assign the value and set the field_index for return
!        for non-null fields (index_t > 0)
!

    if (index_t .gt. 0) then  !{
      temp_field_p%l_value(index_t) = value
      if (index_t .gt. temp_field_p%max_index) then  !{
        temp_field_p%max_index = index_t
      endif  !}
    endif  !}
    field_index = temp_field_p%index
  else  !}{
!
!        Error in making the field
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Could not create logical value field ',             &
           trim(name)
    endif  !}
    field_index = NO_FIELD
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                       &
         'Could not follow path for ',                          &
         trim(name)
  endif  !}
  field_index = NO_FIELD
endif  !}

end function  fm_new_value_logical  !}

!#######################################################################
!#######################################################################

function  fm_new_value_real(name, value, create, index, append) &
          result (field_index)  !{
!
!        Function definition
!
integer                                :: field_index
!
!        arguments
!
character(len=*), intent(in)           :: name
real,             intent(in)           :: value
logical,          intent(in), optional :: create
integer,          intent(in), optional :: index
logical,          intent(in), optional :: append
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=17), parameter :: sub_name     = 'fm_new_value_real'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

logical                          :: create_t
integer                          :: i, ier
integer                          :: index_t
real, pointer, dimension(:)      :: temp_r_value
character(len=fm_path_name_len)  :: path
character(len=fm_field_name_len) :: base
type (field_def), pointer, save  :: temp_list_p 
type (field_def), pointer, save  :: temp_field_p 
integer                          :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  field_index = NO_FIELD
  return
endif  !}
!
!        Check for optional arguments
!
if (present(create)) then  !{
  create_t = create
else  !}{
  create_t = .false.
endif  !}
!
!        Check that append is not true and index greater than 0
!
if (present(index) .and. present(append)) then  !{
  if (append .and. index .gt. 0) then  !{
    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                          &
           'Index and Append both set for ', trim(name)
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
endif  !}
!
!        Set index to define
!

if (present(index)) then  !{
  index_t = index
  if (index_t .lt. 0) then  !{
!
!        Index is negative
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Optional index for ', trim(name),                   &
           ' negative: ', index_t
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
else  !}{
  index_t = 1
endif !}

!
!        Get a pointer to the parent list
!
call find_base(name, path, base)
temp_list_p => find_list(path, current_list_p, create_t)

if (associated(temp_list_p)) then  !{
  temp_field_p => find_field(base, temp_list_p)
  if (.not. associated(temp_field_p)) then  !{
!
!        Create the field if it doesn't exist
!
    temp_field_p => create_field(temp_list_p, base)
  endif  !}
  if (associated(temp_field_p)) then  !{
!
!        Check if the field_type is the same as previously
!        If not then reset max_index to 0
!
    if (temp_field_p%field_type /= real_type ) then
        temp_field_p%max_index = 0
      if (temp_field_p%field_type /= null_type ) then  !{
        if (verb .gt. verb_level_warn) then  !{
          write (out_unit,*) trim(warn_header),                   &
               'Changing type of ', trim(name), ' from ',         &
               trim(field_type_name(temp_field_p%field_type)),    &
               ' to ', trim(field_type_name(real_type))
        endif  !}
      endif  !}
    endif
!
!        Assign the type
!
    temp_field_p%field_type = real_type
!
!        Set the index if appending
!
    if (present(append)) then  !{
      if (append) then  !{
        index_t = temp_field_p%max_index + 1
      endif  !}
    endif  !}
    if (index_t .gt. temp_field_p%max_index + 1) then  !{
!
!        Index too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Index too large for ', trim(name), ': ', index_t
      endif  !}
      field_index = NO_FIELD
      return
    elseif (index_t .eq. 0 .and.                                &
            temp_field_p%max_index .gt. 0) then  !}{
!
!        Can't set non-null field to null
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Trying to nullify a non-null field: ',            &
             trim(name)
      endif  !}
      field_index = NO_FIELD
      return
    elseif (.not. associated(temp_field_p%r_value) .and.        &
            index_t .gt. 0) then  !}{
!
!        Array undefined, so allocate the array
!
      allocate(temp_field_p%r_value(1))
      temp_field_p%max_index = 1
      temp_field_p%array_dim = 1
    elseif (index_t .gt. temp_field_p%array_dim) then  !}{
!
!        Array is too small, so allocate new array and copy over
!        old values
!
      temp_field_p%array_dim = temp_field_p%array_dim + array_increment
      allocate (temp_r_value(temp_field_p%array_dim))
      do i = 1, temp_field_p%max_index  !{
        temp_r_value(i) = temp_field_p%r_value(i)
      enddo  !} i
      if (associated(temp_field_p%r_value)) deallocate(temp_field_p%r_value)
      temp_field_p%r_value => temp_r_value
      temp_field_p%max_index = index_t
    endif  !}
!
!        Assign the value and set the field_index for return
!        for non-null fields (index_t > 0)
!
    if (index_t .gt. 0) then  !{
      temp_field_p%r_value(index_t) = value
      if (index_t .gt. temp_field_p%max_index) then  !{
        temp_field_p%max_index = index_t
      endif  !}
    endif  !}
    field_index = temp_field_p%index
  else  !}{
!
!        Error in making the field
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                        &
           'Could not create real value field ', trim(name)
    endif  !}
    field_index = NO_FIELD
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                          &
         'Could not follow path for ', trim(name)
  endif  !}
  field_index = NO_FIELD
endif  !}

end function  fm_new_value_real  !}

!#######################################################################
!#######################################################################

function  fm_new_value_string(name, value, create, index, append) &
          result (field_index)  !{
!
!        Function definition
!
integer                                :: field_index
!
!        arguments
!
character(len=*), intent(in)           :: name
character(len=*), intent(in)           :: value
logical,          intent(in), optional :: create
integer,          intent(in), optional :: index
logical,          intent(in), optional :: append
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=19), parameter :: sub_name     = 'fm_new_value_string'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=fm_string_len), dimension(:), pointer :: temp_s_value
character(len=fm_path_name_len)                     :: path
character(len=fm_field_name_len)                    :: base
integer                                             :: i, ier
integer                                             :: index_t
logical                                             :: create_t
type (field_def),                     save, pointer :: temp_list_p
type (field_def),                     save, pointer :: temp_field_p
integer                         :: out_unit

out_unit = stdout()
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Must supply a field name
!
if (name .eq. ' ') then  !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'Must supply a field name'
  endif  !}
  field_index = NO_FIELD
  return
endif  !}
!
!        Check for optional arguments
!
if (present(create)) then  !{
  create_t = create
else  !}{
  create_t = .false.
endif  !}
!
!        Check that append is not true and index greater than 0
!

if (present(index) .and. present(append)) then  !{
  if (append .and. index .gt. 0) then  !{
    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Index and Append both set for ', trim(name)
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
endif  !}
!
!        Set index to define
!
if (present(index)) then  !{
  index_t = index
  if (index_t .lt. 0) then  !{
!
!        Index is negative
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Optional index for ', trim(name),                   &
           ' negative: ', index_t
    endif  !}
    field_index = NO_FIELD
    return
  endif  !}
else  !}{
  index_t = 1
endif  !}

!
!        Get a pointer to the parent list
!
call find_base(name, path, base)
temp_list_p => find_list(path, current_list_p, create_t)

if (associated(temp_list_p)) then  !{
  temp_field_p => find_field(base, temp_list_p)
  if (.not. associated(temp_field_p)) then  !{
!
!        Create the field if it doesn't exist
!
    temp_field_p => create_field(temp_list_p, base)
  endif  !}
  if (associated(temp_field_p)) then  !{
!
!        Check if the field_type is the same as previously
!        If not then reset max_index to 0
!
    if (temp_field_p%field_type /= string_type ) then
        temp_field_p%max_index = 0
      if (temp_field_p%field_type /= null_type ) then  !{
        if (verb .gt. verb_level_warn) then  !{
          write (out_unit,*) trim(warn_header),                   &
               'Changing type of ', trim(name), ' from ',         &
               trim(field_type_name(temp_field_p%field_type)),    &
               ' to ', trim(field_type_name(string_type))
        endif  !}
      endif  !}
    endif
!
!        Assign the type
!
    temp_field_p%field_type = string_type
!
!        Set the index if appending
!

    if (present(append)) then  !{
      if (append) then  !{
        index_t = temp_field_p%max_index + 1
      endif  !}
    endif  !}

    if (index_t .gt. temp_field_p%max_index + 1) then  !{

!
!        Index too large
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Index too large for ', trim(name), ': ', index_t
      endif  !}
      field_index = NO_FIELD
      return

    elseif (index_t .eq. 0 .and.                                &
            temp_field_p%max_index .gt. 0) then  !}{

!
!        Can't set non-null field to null
!

      if (verb .gt. verb_level_warn) then  !{
        write (out_unit,*) trim(warn_header),                   &
             'Trying to nullify a non-null field: ',            &
             trim(name)
      endif  !}
      field_index = NO_FIELD
      return

    elseif (.not. associated(temp_field_p%s_value) .and.        &
            index_t .gt. 0) then  !}{

!
!        Array undefined, so allocate the array
!

      allocate(temp_field_p%s_value(1))
      temp_field_p%max_index = 1
      temp_field_p%array_dim = 1

    elseif (index_t .gt. temp_field_p%array_dim) then  !}{

!
!        Array is too small, so allocate new array and copy over
!        old values
!
      temp_field_p%array_dim = temp_field_p%array_dim + array_increment
      allocate (temp_s_value(temp_field_p%array_dim))
      do i = 1, temp_field_p%max_index  !{
        temp_s_value(i) = temp_field_p%s_value(i)
      enddo  !} i
      if (associated(temp_field_p%s_value)) deallocate(temp_field_p%s_value)
      temp_field_p%s_value => temp_s_value
      temp_field_p%max_index = index_t

    endif  !}

!
!        Assign the value and set the field_index for return
!        for non-null fields (index_t > 0)
!

    if (index_t .gt. 0) then  !{
      temp_field_p%s_value(index_t) = value
      if (index_t .gt. temp_field_p%max_index) then  !{
        temp_field_p%max_index = index_t
      endif  !}
    endif  !}
    field_index = temp_field_p%index
  else  !}{
!
!        Error in making the field
!

    if (verb .gt. verb_level_warn) then  !{
      write (out_unit,*) trim(warn_header),                     &
           'Could not create string value field ',              &
           trim(name)
    endif  !}
    field_index = NO_FIELD
  endif  !}
else  !}{
!
!        Error following the path
!

  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                       &
         'Could not follow path for ', trim(name)
  endif  !}
  field_index = NO_FIELD
endif  !}

end function  fm_new_value_string  !}
! </FUNCTION> NAME="fm_new_value"


!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_reset_loop">
!
! <OVERVIEW>
!    Resets the loop variable. For use in conjunction with fm_loop_over_list.
! </OVERVIEW>
! <DESCRIPTION>
!    Resets the loop variable. For use in conjunction with fm_loop_over_list.
! </DESCRIPTION>
!   <TEMPLATE>
!     call fm_reset_loop 
!   </TEMPLATE>
!
subroutine  fm_reset_loop
!
!        arguments
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=13), parameter :: sub_name     = 'fm_reset_loop'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        Reset the variables
!
loop_list = ' '
nullify(loop_list_p)

end subroutine  fm_reset_loop  !}
! </FUNCTION> NAME="fm_reset_loop"

!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_return_root">
!
! <OVERVIEW>
!    Return the root list to the value at initialization
! </OVERVIEW>
! <DESCRIPTION>
!    Return the root list to the value at initialization. 
!    For use in conjunction with fm_change_root. 
!
!    Users should use this routine before leaving their routine if they 
!    previously used fm_change_root.
! </DESCRIPTION>
!   <TEMPLATE>
!     call fm_return_root
!   </TEMPLATE>
!
subroutine  fm_return_root  !{
!
!        arguments
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=14), parameter :: sub_name     = 'fm_return_root'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
!
!        restore the saved values to the current root
!
root_p%name = save_root_name
root_p%parent => save_root_parent_p
!
!        set the pointer to the original root field
!
root_p => root
!
!        reset the save root name and parent variables
!
save_root_name = ' '
nullify(save_root_parent_p)

end subroutine  fm_return_root  !}
! </FUNCTION> NAME="fm_return_root"

!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="get_field">
!
! <OVERVIEW>
!    Return a pointer to the field if it exists relative to this_list_p,
!    null otherwise
! </OVERVIEW>
! <DESCRIPTION>
!    Return a pointer to the field if it exists relative to this_list_p,
!    null otherwise
! </DESCRIPTION>
!   <TEMPLATE>
!     list_p => get_field(name, this_list_p)
!   </TEMPLATE>
!
function get_field(name, this_list_p)                                        &
        result (list_p)  !{
!   <OUT NAME="list_p" TYPE="type (field_def)">
!     A pointer to the field name.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list that the user wishes to get information for.
!   </IN>
!   <IN NAME="this_list_p" TYPE="type (field_def)">
!     A pointer to a list that serves as the base point for searching for name.
!   </IN>
!
!        Function definition
!
type (field_def), pointer        :: list_p
!
!        arguments
!
character(len=*), intent(in)     :: name
type (field_def), pointer        :: this_list_p
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=9),  parameter :: sub_name     = 'get_field'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len)  :: path
character(len=fm_field_name_len) :: base
type (field_def), pointer, save  :: temp_p 

nullify(list_p)
!
!        Get the path and base for name
!
call find_base(name, path, base)
!
!        Find the list if path is not empty
!
if (path .ne. ' ') then  !{
  temp_p => find_list(path, this_list_p, .false.)
  if (associated(temp_p)) then  !{
    list_p => find_field(base, temp_p)
  else  !}{
    nullify(list_p)
  endif  !}
else  !}{
  list_p => find_field(base, this_list_p)
endif  !}

end function get_field  !}
! </FUNCTION> NAME="get_field"
!</PRIVATE>


!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_modify_name">
!
! <OVERVIEW>
!    This function allows a user to rename a field without modifying the 
!    contents of the field.
! </OVERVIEW>
! <DESCRIPTION>
!    Function to modify the name of a field. 
!    Should be used with caution.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_modify_name(oldname, newname)
!   </TEMPLATE>
!
function fm_modify_name(oldname, newname)                                        &
        result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="oldname" TYPE="character(len=*)">
!     The name of a field that the user wishes to change the name of.
!   </IN>
!   <IN NAME="newname" TYPE="character(len=*)">
!     The name that the user wishes to change the name of the field to.
!   </IN>
!
!        Function definition
!
logical                          :: success
!
!        arguments
!
character(len=*), intent(in)     :: oldname
character(len=*), intent(in)     :: newname
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=14), parameter :: sub_name     = 'fm_modify_name'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len)  :: path
character(len=fm_field_name_len) :: base
type (field_def), pointer, save  :: list_p 
type (field_def), pointer, save  :: temp_p 
!
!        Get the path and base for name
!
call find_base(oldname, path, base)
!
!        Find the list if path is not empty
!
success = .false.
if (path .ne. ' ') then  !{
  temp_p => find_list(path, current_list_p, .false.)
  if (associated(temp_p)) then  !{
    list_p => find_field(base, temp_p)
    if (associated(list_p)) then !{
      list_p%name = newname
      success = .true.
    endif!}
  else  !}{
    nullify(list_p)
  endif  !}
else  !}{
  list_p => find_field(base, current_list_p)
  if (associated(list_p)) then !{
    list_p%name = newname
    success = .true.
  endif !} 
endif  !}

end function fm_modify_name  !}
! </FUNCTION> NAME="fm_modify_name"


!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="initialize">
!
! <OVERVIEW>
!    A function to initialize the values of the pointers. This will remove
!    all fields and reset the field tree to only the root field.
! </OVERVIEW>
! <DESCRIPTION>
!    A function to initialize the values of the pointers. This will remove
!    all fields and reset the field tree to only the root field.
! </DESCRIPTION>
!   <TEMPLATE>
!     call initialize
!   </TEMPLATE>
!
subroutine initialize  !{
!
!        arguments
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=10), parameter :: sub_name     = 'initialize'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer :: ier
!
!        Initialize the root field
!
if (.not. module_is_initialized) then  !{
  root_p => root

  field_type_name(integer_type) = 'integer'
  field_type_name(list_type) = 'list'
  field_type_name(logical_type) = 'logical'
  field_type_name(real_type) = 'real'
  field_type_name(string_type) = 'string'

  root%name = ' '
  root%index = 1
  root%parent => root_p

  root%field_type = list_type

  root%length = 0
  nullify(root%first_field)
  nullify(root%last_field)
  root%max_index = 0
  root%array_dim = 0
  if (associated(root%i_value)) deallocate(root%i_value)
  if (associated(root%l_value)) deallocate(root%l_value)
  if (associated(root%r_value)) deallocate(root%r_value)
  if (associated(root%s_value)) deallocate(root%s_value)

  nullify(root%next)
  nullify(root%prev)

  current_list_p => root

  nullify(loop_list_p)
  loop_list = ' '

  nullify(save_root_parent_p)
  save_root_name = ' '

  module_is_initialized = .true.

endif  !}

end subroutine initialize  !}
! </FUNCTION> NAME="initialize"
!</PRIVATE>

!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="make_list">
!
! <OVERVIEW>
!    This function creates a new field and returns a pointer to that field.
! </OVERVIEW>
! <DESCRIPTION>
!    Allocate and initialize a new list in this_list_p list.
!    Return a pointer to the list on success, or a null pointer
!    on failure
! </DESCRIPTION>
!   <TEMPLATE>
!     list_p => make_list(this_list_p, name)
!   </TEMPLATE>
!
function  make_list(this_list_p, name)                        &
          result (list_p)  !{
!   <OUT NAME="list_p" TYPE="type (field_def), pointer">
!     A pointer to the list that has been created. 
!   </OUT>
!   <IN NAME="this_list_p" TYPE="type (field_def), pointer">
!     The base of a list that the user wishes to add a list to.
!   </IN>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list that the user wishes to create.
!   </IN>
!
!        Function definition
!
type (field_def), pointer    :: list_p
!
!        arguments
!
type (field_def), pointer    :: this_list_p
character(len=*), intent(in) :: name
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=9),  parameter :: sub_name     = 'make_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer :: ier
type (field_def), pointer, save :: dummy_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Check to see whether there is already a list with
!        this name, and if so, return an error as list names
!        must be unique
!
dummy_p => find_field(name, this_list_p )
if (associated(dummy_p)) then  !{
!
!        This list is already specified, return an error
!
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header), 'List ',                 &
         trim(name), ' already exists'
  endif  !}
!  nullify(list_p)
  list_p => dummy_p
  return
endif  !}
!
!        Create a field for the new list
!
nullify(list_p)
list_p => create_field(this_list_p, name)
if (.not. associated(list_p)) then !{
  if (verb .gt. verb_level_warn) then  !{
    write (out_unit,*) trim(warn_header),                          &
         'Could not create field ', trim(name)
  endif  !}
  nullify(list_p)
  return
endif  !}
!
!        Initialize the new list
!
list_p%length = 0
list_p%field_type = list_type
if (associated(list_p%i_value)) deallocate(list_p%i_value)
if (associated(list_p%l_value)) deallocate(list_p%l_value)
if (associated(list_p%r_value)) deallocate(list_p%r_value)
if (associated(list_p%s_value)) deallocate(list_p%s_value)

end function  make_list  !}
! </FUNCTION> NAME="make_list"
!</PRIVATE>


!#######################################################################
!#######################################################################

! <FUNCTION NAME="fm_query_method">
!
! <OVERVIEW>
!    This is a function that provides the capability to return parameters 
!    associated with a field in a pair of strings.
! </OVERVIEW>
! <DESCRIPTION>
!    Given a name return a list of method names and control strings.
!    This function should return strings similar to those in the field
!    table if a comma delimited format is being used.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_query_method(name, method_name, method_control)
!   </TEMPLATE>
!
function fm_query_method(name, method_name, method_control)                &
          result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list that the user wishes to change to.
!   </IN>
!   <OUT NAME="method_name" TYPE="character(len=*)">
!     The name of a parameter associated with the named field.
!   </OUT>
!   <OUT NAME="method_control" TYPE="character(len=*)">
!     The value of parameters associated with the named field.
!   </OUT>
!
!        Function definition
!
logical                       :: success
!
!        arguments
!
character(len=*), intent(in)  :: name
character(len=*), intent(out) :: method_name
character(len=*), intent(out) :: method_control
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=15), parameter :: sub_name     = 'fm_query_method'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len) :: path
character(len=fm_path_name_len) :: base
character(len=fm_path_name_len) :: name_loc
logical                         :: recursive_t
type (field_def), pointer, save :: temp_list_p 
type (field_def), pointer, save :: temp_value_p 
type (field_def), pointer, save :: this_field_p 
integer                         :: out_unit

  out_unit = stdout()
  success     = .false.
  recursive_t = .true.
  method_name = " "
  method_control = " "
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}
name_loc = lowercase(name)
call find_base(name_loc, path, base)

  temp_list_p => find_list(name_loc, current_list_p, .false.)

if (associated(temp_list_p)) then
! Find the entry values for the list.
  success = query_method(temp_list_p, recursive_t, base, method_name, method_control)
else  !}{
! This is not a list but it may be a parameter with a value
! If so put the parameter value in method_name.

  temp_value_p => find_list(path, current_list_p, .false.)
  if (associated(temp_value_p)) then  !{
! Find the entry values for this item.
  this_field_p => temp_value_p%first_field

  do while (associated(this_field_p))  !{
    if ( this_field_p%name == base ) then !{
      method_name = this_field_p%s_value(1)
      method_control = ""
      success = .true.
      exit
    else !}{
      success = .false.
    endif !}
    this_field_p => this_field_p%next
  enddo

  else  !}{
!
!        Error following the path
!
    if (verb .gt. verb_level_warn) then
      write (out_unit,*) trim(warn_header), 'Could not follow path for ', trim(path)
    endif
    success = .false.
  endif  !}
endif  !}

end function  fm_query_method  !}
! </FUNCTION> NAME="fm_query_method"

!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME="query_method">
!
! <OVERVIEW>
!    A private function that can recursively recover values for parameters 
!    associated with a field.
! </OVERVIEW>
! <DESCRIPTION>
!    A private function that can recursively recover values for parameters 
!    associated with a field.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = query_method(list_p, recursive, name, method_name, method_control)
!   </TEMPLATE>
!
recursive function query_method(list_p, recursive, name, method_name, method_control) &
          result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="list_p" TYPE="type (field_def), pointer">
!     A pointer to the field that is of interest.
!   </IN>
!   <IN NAME="name" TYPE="character(len=*)">
!     The name of a list that the user wishes to change to.
!   </IN>
!   <OUT NAME="method_name" TYPE="character(len=*)">
!     The name of a parameter associated with the named field.
!   </OUT>
!   <OUT NAME="method_control" TYPE="character(len=*)">
!     The value of parameters associated with the named field.
!   </OUT>
!
!        Function definition
!
logical                       :: success
!
!        arguments
!
type (field_def), pointer     :: list_p
logical,          intent(in)  :: recursive
character(len=*), intent(in)  :: name
character(len=*), intent(out) :: method_name, method_control
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=12), parameter :: sub_name     = 'query_method'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
integer,                  parameter :: max_depth = 64
character(len=max_depth), parameter :: blank = '    '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: depthp1
integer                         :: first
integer                         :: i
integer                         :: last
character(len=64)               :: scratch
type (field_def), pointer :: this_field_p 
integer                         :: out_unit

out_unit = stdout()

!
!        Check for a valid list
!
if (.not. associated(list_p)) then  !{
  if (verb .gt. verb_level_warn) then
    write (out_unit,*) trim(warn_header), 'Invalid list pointer'
  endif
  success = .false.
elseif (list_p%field_type .ne. list_type) then  !}{
  if (verb .gt. verb_level_warn) then
    write (out_unit,*) trim(warn_header), trim(list_p%name)//' is not a list'
  endif
  success = .false.
else  !}{
!
!        set the default return value
!
  success = .true.

  this_field_p => list_p%first_field

  do while (associated(this_field_p))  !{
    select case(this_field_p%field_type)
    case(list_type)
!
!        If this is a list, then this is the method name
!
      if (recursive) then  !{
        if (.not. query_method(this_field_p, .true., this_field_p%name, method_name, method_control)) then  !{
          success = .false.
          exit
        else  !}{
          !write (method_name,'(a,a)') method_name(1:LEN_TRIM(method_name)), &
          i = LEN_TRIM(method_name)
          if ( i .gt. 0 ) then
            write (method_name,'(a,a)') method_name(1:i), &
                    trim(this_field_p%name)
          else
            write (method_name,'(a)') trim(this_field_p%name)
          endif        
        endif  !}
      endif  !}

    case(integer_type)
        write (scratch,*) this_field_p%i_value
        call strip_front_blanks(scratch)
        write (method_control,'(a,a,a,a,a)') trim(method_control),comma, &
                trim(this_field_p%name), ' = ', trim(scratch)


    case(logical_type)

        write (method_control,'(a,a,a,a,l1)') trim(method_control),comma, &
                trim(this_field_p%name), ' = ', this_field_p%l_value

    case(real_type)

        write (scratch,*) this_field_p%r_value
        call strip_front_blanks(scratch)
        write (method_control,'(a,a,a,a,a)') trim(method_control),comma, &
                trim(this_field_p%name), ' = ', trim(scratch)


    case(string_type)
        write (method_control,'(a,a,a,a,a,$)') trim(method_control),comma, &
                trim(this_field_p%name), ' = ',trim(this_field_p%s_value(1))
        do i = 2, this_field_p%max_index
          write (method_control,'(a,a,$)') comma//trim(this_field_p%s_value(i))
        enddo


    case default
        if (verb .gt. verb_level_warn) then
          write (out_unit,*) trim(warn_header), 'Undefined type for ', trim(this_field_p%name)
        endif
        success = .false.
        exit

    end select 
    this_field_p => this_field_p%next
  enddo  !}
endif  !}

end function query_method  !}
! </FUNCTION> NAME="query_method"
!</PRIVATE>

!#######################################################################
!#######################################################################

! <FUNCTION NAME = "fm_copy_list" >
! <OVERVIEW>
!    A function that allows the user to copy a field and add a suffix to 
!    the name of the new field.
! </OVERVIEW>
! <DESCRIPTION>
!    Given the name of a pre-existing field and a suffix, this function
!    will create a new field. The name of the new field will be that of 
!    the old field with a suffix supplied by the user.
! </DESCRIPTION>
!   <TEMPLATE>
!     index = fm_copy_list(list_name, suffix, create)
!   </TEMPLATE>
!
function fm_copy_list(list_name, suffix, create ) &
         result(index)   !{
!   <OUT NAME="index" TYPE="integer">
!     The index of the field that has been created by the copy.
!   </OUT>
!   <IN NAME="list_name" TYPE="character(len=*)">
!     The name of a field that the user wishes to copy..
!   </IN>
!   <IN NAME="suffix" TYPE="character(len=*)">
!     The suffix that will be added to list_name when the field is copied.
!   </IN>
!
!        Function definition
!
integer        :: index
!
!        arguments
!
character(len=*), intent(in)           :: list_name
character(len=*), intent(in)           :: suffix
logical,          intent(in), optional :: create

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=12), parameter :: sub_name     = 'fm_copy_list'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_string_len), dimension(MAX_FIELD_METHODS) :: control
character(len=fm_string_len), dimension(MAX_FIELD_METHODS) :: method
character(len=fm_string_len)                               :: head
character(len=fm_string_len)                               :: list_name_new
character(len=fm_string_len)                               :: tail
character(len=fm_string_len)                               :: val_str
integer                                                    :: n
integer                                                    :: num_meth
integer                                                    :: val_int
logical                                                    :: found_methods
logical                                                    :: got_value
logical                                                    :: recursive_t
logical                                                    :: success
logical                                                    :: val_logical
real                                                       :: val_real
type (field_def), pointer, save                            :: temp_field_p 
type (field_def), pointer, save                            :: temp_list_p 
integer                                                    :: out_unit

out_unit = stdout()


num_meth= 1
list_name_new = trim(list_name)//trim(suffix)
!
  recursive_t = .true.
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}

if (list_name .eq. ' ') then  !{
!
!        If list is empty, then dump the current list
!
  temp_list_p => current_list_p
  success = .true.
else  !}{
!
!        Get a pointer to the list
!
  temp_list_p => find_list(list_name, current_list_p, .false.)
  if (associated(temp_list_p)) then  !{
    success = .true.
  else  !}{
!
!        Error following the path
!
    if (verb .gt. verb_level_warn) then
      write (out_unit,*) trim(warn_header), 'Could not follow path for ', trim(list_name)
    endif
    success = .false.
  endif  !}
endif  !}

!
!        Find the list
!
if (success) then  !{
  method(:) = ' '
  control(:) = ' '
  found_methods = fm_find_methods(trim(list_name), method, control)
  do n = 1, MAX_FIELD_METHODS
    if (LEN_TRIM(method(n)) > 0 ) then
      index = fm_new_list(trim(list_name_new)//list_sep//method(n), create = create)
      call find_base(method(n), head, tail)
      temp_field_p => find_list(trim(list_name)//list_sep//head,temp_list_p, .false.)
      temp_field_p => find_field(tail,temp_field_p)
      select case (temp_field_p%field_type)
        case (integer_type)
          got_value = fm_get_value( trim(list_name)//list_sep//method(n), val_int)
          if ( fm_new_value( trim(list_name_new)//list_sep//method(n), val_int, &
                             create = create, append = .true.) < 0 ) &
            call mpp_error(FATAL, trim(error_header)//'Could not set the '//trim(method(n))//&
                                  ' for '//trim(list_name)//trim(suffix))
  
        case (logical_type)
          got_value = fm_get_value( trim(list_name)//list_sep//method(n), val_logical)
          if ( fm_new_value( trim(list_name_new)//list_sep//method(n), val_logical, &
                             create = create, append = .true.) < 0 ) &
            call mpp_error(FATAL, trim(error_header)//'Could not set the '//trim(method(n))//&
                                  ' for '//trim(list_name)//trim(suffix))
  
        case (real_type)
          got_value = fm_get_value( trim(list_name)//list_sep//method(n), val_real)
          if ( fm_new_value( trim(list_name_new)//list_sep//method(n), val_real, &
                             create = create, append = .true.) < 0 ) &
            call mpp_error(FATAL, trim(error_header)//'Could not set the '//trim(method(n))//&
                                  ' for '//trim(list_name)//trim(suffix))
  
        case (string_type)
          got_value = fm_get_value( trim(list_name)//list_sep//method(n), val_str)
          if ( fm_new_value( trim(list_name_new)//list_sep//method(n), val_str, &
                             create = create, append = .true.) < 0 ) &
            call mpp_error(FATAL, trim(error_header)//'Could not set the '//trim(method(n))//&
                                  ' for '//trim(list_name)//trim(suffix))
        case default
      end select
  
    endif
  enddo
endif  !}

end function fm_copy_list !}         
! </FUNCTION > NAME = "fm_copy_list"

!#######################################################################
!#######################################################################

! <FUNCTION NAME = "fm_find_methods" >
! <OVERVIEW>
!    This function retrieves all the methods associated with a field.
! </OVERVIEW>
! <DESCRIPTION>
!    This function retrieves all the methods associated with a field.
!    This is different from fm_query_method in that this function gets all
!    the methods associated as opposed to 1 method.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = fm_find_methods(list_name, methods, control )
!   </TEMPLATE>
!
function fm_find_methods(list_name, methods, control ) &
         result(success)   !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="list_name" TYPE="character(len=*)">
!     The name of a list that the user wishes to find methods for.
!   </IN>
!   <OUT NAME="methods" TYPE="character(len=*)">
!     An array of the methods associated with list_name.
!   </OUT>
!   <OUT NAME="control" TYPE="character(len=*)">
!     An array of the parameters associated with methods.
!   </OUT>
!
!        Function definition
!
logical                                     :: success
!
!        arguments
!
character(len=*), intent(in)                :: list_name
character(len=*), intent(out), dimension(:) :: methods
character(len=*), intent(out), dimension(:) :: control

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=15), parameter :: sub_name     = 'fm_find_methods'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: num_meth
logical                         :: recursive_t
type (field_def), pointer, save :: temp_list_p 
integer                         :: out_unit

out_unit = stdout()
num_meth= 1
!
!        Check whether to do things recursively
!
  recursive_t = .true.
!  recursive_t = .false.
!
!        Initialize the field manager if needed
!
if (.not. module_is_initialized) then  !{
  call initialize
endif  !}

if (list_name .eq. ' ') then  !{
!
!        If list is empty, then dump the current list
!
  temp_list_p => current_list_p
  success = .true.
else  !}{
!
!        Get a pointer to the list
!
  temp_list_p => find_list(list_name, current_list_p, .false.)
  if (associated(temp_list_p)) then  !{
    success = .true.
  else  !}{
!
!        Error following the path
!
    if (verb .gt. verb_level_warn) then
      write (out_unit,*) trim(warn_header), 'Could not follow path for ', trim(list_name)
    endif
    success = .false.
  endif  !}
endif  !}

!
!        Find the list
!
if (success) then  !{
  success = find_method(temp_list_p, recursive_t, num_meth, methods, control)
endif  !}

end function fm_find_methods !}         
! </FUNCTION > NAME = "fm_find_methods"

!#######################################################################
!#######################################################################

! <PRIVATE><FUNCTION NAME = "find_method">
!
! <OVERVIEW>
!    Given a field list pointer this function retrieves methods and 
!    associated parameters for the field list.
! </OVERVIEW>
! <DESCRIPTION>
!    Given a field list pointer this function retrieves methods and 
!    associated parameters for the field list.
! </DESCRIPTION>
!   <TEMPLATE>
!     success = find_method(list_p, recursive, num_meth, method, control)
!   </TEMPLATE>
!
recursive function find_method(list_p, recursive, num_meth, method, control)   &
          result (success)  !{
!   <OUT NAME="success" TYPE="logical">
!     A flag to indicate whether the function operated with (FALSE) or 
!     without (TRUE) errors.
!   </OUT>
!   <IN NAME="list_p" TYPE="type (field_def), pointer">
!     A pointer to the field of interest
!   </IN>
!   <IN NAME="recursive" TYPE="logical">
!     If true, then recursively search for methods.
!   </IN>
!   <INOUT NAME="num_meth" TYPE="integer">
!     The number of methods found.
!   </INOUT>
!   <OUT NAME="method" TYPE="character(len=*)" DIM="(:)">
!     The methods associated with the field pointed to by list_p
!   </OUT>
!   <OUT NAME="control" TYPE="character(len=*)" DIM="(:)">
!     The control parameters for the methods found.
!   </OUT>
!
!        Function definition
!
logical                                     :: success
!
!        arguments
!
type (field_def), pointer                   :: list_p
logical,          intent(in)                :: recursive
integer,          intent(inout)             :: num_meth
character(len=*), intent(out), dimension(:) :: method
character(len=*), intent(out), dimension(:) :: control
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=11), parameter :: sub_name     = 'find_method'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer, parameter                          :: max_depth = 64
character(len=max_depth), parameter         :: blank = '    '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
character(len=fm_path_name_len) :: scratch
integer                         :: depthp1
integer                         :: first
integer                         :: i
integer                         :: last
integer                         :: n
type (field_def), pointer, save :: this_field_p 
integer                         :: out_unit

out_unit = stdout()
!
!        Check for a valid list
!
if (.not. associated(list_p)) then  !{
  if (verb .gt. verb_level_warn) then
    write (out_unit,*) trim(warn_header), 'Invalid list pointer'
  endif
  success = .false.
elseif (list_p%field_type .ne. list_type) then  !}{
  if (verb .gt. verb_level_warn) then
    write (out_unit,*) trim(warn_header), trim(list_p%name), ' is not a list'
  endif
  success = .false.
else  !}{
!
!        set the default return value
!
  success = .true.

  this_field_p => list_p%first_field

  do while (associated(this_field_p))  !{
    select case(this_field_p%field_type)
    case(list_type)
!
!        If this is a list, then this is the method name
!
        if ( this_field_p%length > 1) then
           do n = num_meth+1, num_meth + this_field_p%length - 1
              write (method(n),'(a,a,a,$)') trim(method(num_meth)), &
                                            trim(this_field_p%name), list_sep
           enddo
           write (method(num_meth),'(a,a,a,$)') trim(method(num_meth)), &
                                                trim(this_field_p%name), list_sep
        else
           write (method(num_meth),'(a,a,a,$)') trim(method(num_meth)), &
                                                trim(this_field_p%name), list_sep
        endif
        success = find_method(this_field_p, .true., num_meth, method, control)

    case(integer_type)
        write (scratch,*) this_field_p%i_value
        call strip_front_blanks(scratch)
        write (method(num_meth),'(a,a)') trim(method(num_meth)), &
                trim(this_field_p%name)
        write (control(num_meth),'(a)') &
                trim(scratch)
        num_meth = num_meth + 1


    case(logical_type)

        write (method(num_meth),'(a,a)') trim(method(num_meth)), &
                trim(this_field_p%name)
        write (control(num_meth),'(l1)') &
                this_field_p%l_value
        num_meth = num_meth + 1

    case(real_type)

        write (scratch,*) this_field_p%r_value
        call strip_front_blanks(scratch)
        write (method(num_meth),'(a,a)') trim(method(num_meth)), &
                trim(this_field_p%name)
        write (control(num_meth),'(a)') &
                trim(scratch)
        num_meth = num_meth + 1


    case(string_type)
        write (method(num_meth),'(a,a)') trim(method(num_meth)), &
                trim(this_field_p%name)
        write (control(num_meth),'(a)') &
                 trim(this_field_p%s_value(1))
        do i = 2, this_field_p%max_index
          write (control(num_meth),'(a,a,$)') comma//trim(this_field_p%s_value(i))
        enddo
        num_meth = num_meth + 1


    case default
        if (verb .gt. verb_level_warn) then
          write (out_unit,*) trim(warn_header), 'Undefined type for ', trim(this_field_p%name)
        endif
        success = .false.
        exit

    end select 

    this_field_p => this_field_p%next
  enddo  !}
endif  !}

end function find_method !}
! </FUNCTION > NAME = "find_method"
!</PRIVATE>

!#######################################################################
! <SUBROUTINE NAME="fm_set_verbosity">
!
! <OVERVIEW>
!   A subroutine to set the verbosity of the field manager output.
! </OVERVIEW>
! <DESCRIPTION>
!   This subroutine will set the level of verbosity in the module.
!   Currently, verbosity is either on (1) or off (0). However,
!   in the future, "on" may have more granularity. If no argument
!   is given, then, if verbosity is on it will be turned off, and
!   is off, will be turned to the default on level.
!   If verbosity is negative then it is turned off.
!   Values greater than the maximum will be set to the maximum.
! </DESCRIPTION>
!   <TEMPLATE>
!     call fm_set_verbosity(verbosity)
!   </TEMPLATE>
!
subroutine  fm_set_verbosity(verbosity)  !{
!   <IN NAME="verbosity" TYPE="integer, optional">
!     The level of verbosity required by the user.
!   </IN>
!
!       arguments
!

integer, intent(in), optional :: verbosity

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!       local parameters
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

character(len=16), parameter :: sub_name     = 'fm_set_verbosity'
character(len=64), parameter :: error_header = '==>Error from ' // trim(module_name)   //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: warn_header  = '==>Warning from ' // trim(module_name) //  &
                                               '(' // trim(sub_name) // '): '
character(len=64), parameter :: note_header  = '==>Note from ' // trim(module_name)    //  &
                                               '(' // trim(sub_name) // '): '
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!        local variables
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
integer                         :: out_unit

out_unit = stdout()

!
!       Check whether an argument has been given
!

if (present(verbosity)) then  !{

  if (verbosity .le. 0) then  !{
    verb = 0
  elseif (verbosity .ge. max_verbosity) then  !}{
    verb = max_verbosity
  else  !}{
    verb = verbosity
  endif  !}

else  !}{

  if (verb .eq. 0) then  !{
    verb = default_verbosity
  else  !}{
    verb = 0
  endif  !}

endif  !}

write (out_unit,*) 
write (out_unit,*) trim(note_header),                          &
     'Verbosity now at level ', verb
write (out_unit,*) 

end subroutine  fm_set_verbosity  !}
! </SUBROUTINE> NAME="fm_set_verbosity"

end module field_manager_mod

#ifdef test_field_manager

program test

use field_manager_mod
use mpp_mod, only : mpp_exit, mpp_pe, mpp_root_pe, mpp_error, NOTE

implicit none
!#include "mpif.h"


integer :: i, j, nfields, num_methods, model
character(len=fm_string_len) :: field_type, field_name, str, name_field_type, path
character(len=512) :: method_name, method_control
real :: param
integer :: flag, index
logical :: success
type(method_type), dimension(20) :: methods

call field_manager_init(nfields)

! Dump the list of fields produced from reading the field_table

! Here are the lists that propagate off the root "/"
! By calling fm_dump_list with a single argument you only get the 
! lists branching off this argument in the list.
write(*,*) "Here's a baseline listing"
success = fm_dump_list("/")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! By adding the optional .true. argument you get a recursive listing of the fields.
write(*,*) "Here's a recursive listing"
success = fm_dump_list("/", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Using fm_dump_list with a blank first argument returns the last field accessed by field manager.
write(*,*) 'Dumping last field changed to by field_manager using fm_change_list'
success = fm_dump_list("", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Change list to look at the land model fields
write(*,*) 'Changing list to land_mod'
success = fm_change_list("/land_mod")
write(*,*) 'Dumping last list changed to by field_manager using fm_change_list i.e list of land model fields'
success = fm_dump_list("", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Now let's modify some of the field entries.
! 
!In this example we add a field ( convection = 'off' ) to the radon list
write(*,*) "ADDING convection = off TO RADON LIST"
!if ( fm_change_list('/atmos_mod/tracer/radon')) then
if ( fm_exists('/atmos_mod/tracer/radon')) then
   write(*,*) "'/atmos_mod/tracer/radon' exists "
   success = fm_change_list('/atmos_mod/tracer/radon')
! The next line creates a new field branching off radon.
   index = fm_new_value('convection','off')
endif

success = fm_query_method('radon',method_name,method_control)
if (success ) then
call mpp_error(NOTE, "Method names for radon is/are "//trim(method_name))
call mpp_error(NOTE, "Method controls for radon is/are "//trim(method_control))
else
call mpp_error(NOTE, "There is no atmos model radon field defined in the field_table")
endif
! Dump the listing of the modified tracer
success = fm_dump_list("/atmos_mod/tracer/radon", .true.)
if (.not. success ) call mpp_error(NOTE, "There is no atmos model radon field defined in the field_table")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


! Find out what the current path is. Should be '/atmos_mod/tracer/radon' as set in fm_change_list above.
path = fm_get_current_list()
write(*,*) 'Current path is ',trim(path)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Now let's modify the value of the field we just added.
write(*,*) "MODIFYING RADON FIELD CONVECTION ATTRIBUTE TO convection = RAS_off "
index = fm_new_value('convection','RAS_off')

! Dump the listing of the modified tracer
success = fm_dump_list("/atmos_mod/tracer/radon", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'





write(*,*) "ORIGINAL OCEAN MODEL TRACER FIELDS"

! Dump the listing of the original ocean model tracers
success = fm_dump_list("/ocean_mod/tracer", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


index = fm_get_length("/ocean_mod/tracer") 
write(*,*) "The length of the current list '/ocean_mod/tracer' is ",index," i.e."
success = fm_dump_list("/ocean_mod/tracer")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Find out what type of field this is. Possibilities are real, integer, string, logical, and list
name_field_type = fm_get_type('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope')
write(*,*) 'The type for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is ',name_field_type

success = fm_get_value('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope',str)
write(*,*) 'The value for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is (character) ',str


write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) "MODIFYING BIOTIC1 FIELD slope ATTRIBUTE TO slope = 0.95 "
if ( fm_change_list('/ocean_mod/tracer/biotic1/diff_horiz/linear')) &
   index = fm_new_value('slope',0.95, index = 1)

! Dump the listing of the modified ocean model tracer attribute
success = fm_dump_list("/ocean_mod/tracer/biotic1", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

name_field_type = fm_get_type('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope')
write(*,*) 'Now the type for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is ',name_field_type
success =  fm_get_value('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope',param)
write(*,*) 'The value for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is (real) ',param
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) 'Changing the name of biotic1 to biotic_control'
success = fm_modify_name('/ocean_mod/tracer/biotic1', 'biotic_control')

! Dump the listing of the modified tracer
success = fm_dump_list("/ocean_mod/tracer/biotic_control", .true.)

! Double check to show that the tracer has been renamed and the original doesn't exist anymore. 
success = fm_dump_list("/ocean_mod/tracer/biotic1", .true.)
if (.not. success ) call mpp_error(NOTE, "Ocean model tracer biotic1 does not exist anymore.")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


if ( fm_change_list("/ocean_mod/tracer/age_ctl") ) then
success = fm_dump_list("", .true.)
write(*,*) "Now we'll add a new list to this list"
index = fm_new_list("units",create = .true.)

success = fm_dump_list("", .true.)

write(*,*) "Now we'll give it a value"
if (success) index = fm_new_value('units','days')

success = fm_dump_list("", .true.)


write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
endif
!errorcode = 121
!CALL MPI_ERROR_STRING(errorcode, string, resultlen, ierror)
!write(*,*) string
call field_manager_end

call mpp_exit

end program test

#endif
