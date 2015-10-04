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

module tracer_manager_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Matt Harrison
! </REVIEWER>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Bruce Wyman
! </REVIEWER>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Peter Phillipps
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Code to manage the simple addition of tracers to the FMS code.
!     This code keeps track of the numbers and names of tracers included
!     in a tracer table.
! </OVERVIEW>

! <DESCRIPTION>
!     This code is a grouping of calls which will allow the simple
!     introduction of tracers into the FMS framework. It is designed to
!     allow users of a variety of component models interact easily with
!     the dynamical core of the model. 
!     
!     In calling the tracer manager routines the user must provide a
!     parameter identifying the model that the user is working with. This
!     parameter is defined within field_manager as MODEL_X 
!     where X is one of [ATMOS, OCEAN, LAND, ICE].
!
!     In many of these calls the argument list includes model and tracer_index. These 
!     are the parameter corresponding to the component model and the tracer_index N is 
!     the Nth tracer within the component model. Therefore a call with MODEL_ATMOS and 5 
!     is different from a call with MODEL_OCEAN and 5.
!
! </DESCRIPTION>


!----------------------------------------------------------------------

use           mpp_mod, only : mpp_error,          &
                              mpp_pe,             &
                              mpp_root_pe,        &
                              FATAL,              &
                              WARNING,            &
                              NOTE,               &
                              stdlog
use        mpp_io_mod, only : mpp_open,           &
                              mpp_close,          &
                              MPP_ASCII,          &
                              MPP_APPEND,         &
                              MPP_RDONLY
use           fms_mod, only : lowercase,          &
                              write_version_number

use field_manager_mod, only : field_manager_init, &
                              get_field_info,     &
                              get_field_methods,  &
                              MODEL_ATMOS,        &
                              MODEL_LAND,         &
                              MODEL_OCEAN,        &
                              MODEL_ICE,          &
                              MODEL_COUPLER,      &
                              NUM_MODELS,         &
                              method_type,        &
                              default_method,     &
                              parse,              &
                              fm_copy_list,       &
                              fm_change_list,     &
                              fm_modify_name,     &
                              fm_query_method,    &
                              fm_new_value,       &
                              fm_exists,          &
                              MODEL_NAMES

implicit none
private

!-----------------------------------------------------------------------

public  tracer_manager_init, &
        tracer_manager_end,  &
        check_if_prognostic, &
        get_tracer_indices,  &
        get_tracer_index,    &
        get_tracer_names,    &
        get_tracer_name,     &
        query_method,        &
        set_tracer_atts,     &
        set_tracer_profile,  &
        register_tracers,    &
        get_number_tracers,  &
        NO_TRACER,           &
        MAX_TRACER_FIELDS

!-----------------------------------------------------------------------
interface get_tracer_index
  module procedure get_tracer_index_integer, get_tracer_index_logical
end interface
!-----------------------------------------------------------------------

integer            :: num_tracer_fields = 0
integer, parameter :: MAX_TRACER_FIELDS = 120
integer, parameter :: MAX_TRACER_METHOD = 20
integer, parameter :: NO_TRACER         = 1-HUGE(1)
integer, parameter :: NOTRACER          = -HUGE(1)

integer :: total_tracers(NUM_MODELS), prog_tracers(NUM_MODELS), diag_tracers(NUM_MODELS)
logical :: model_registered(NUM_MODELS) = .FALSE.

type, private ::  tracer_type
   character(len=32)        :: tracer_name, tracer_units
   character(len=128)       :: tracer_longname
   integer                  :: num_methods, model, instances
   logical                  :: is_prognostic, instances_set
   logical                  :: needs_init
end type tracer_type

type, private ::  tracer_name_type
   character(len=32)  :: model_name, tracer_name, tracer_units
   character(len=128) :: tracer_longname
end type tracer_name_type


type, private :: inst_type
   character(len=128) :: name
   integer            :: instances
end type inst_type

type(tracer_type), save  :: tracers(MAX_TRACER_FIELDS)
type(inst_type)  , save  :: instantiations(MAX_TRACER_FIELDS)

character(len=128) :: version = '$Id: tracer_manager.F90,v 16.0 2008/07/30 22:48:11 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'
logical            :: module_is_initialized = .false.

logical            :: verbose_local
integer            :: TRACER_ARRAY(NUM_MODELS,MAX_TRACER_FIELDS)

contains

!
!#######################################################################
!
! <SUBROUTINE NAME="tracer_manager_init">
!   <OVERVIEW>
!      It is not necessary to call this routine.
!      It is included only for backward compatability.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine writes the version and tagname to the logfile and 
!     sets the module initialization flag.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_manager_init
!   </TEMPLATE>
subroutine tracer_manager_init
integer :: model, num_tracers, num_prog, num_diag

  if(module_is_initialized) return
  module_is_initialized = .TRUE.

  call write_version_number (version, tagname)
  call field_manager_init()
  TRACER_ARRAY = NOTRACER
  do model=1,NUM_MODELS 
    call get_tracer_meta_data(model, num_tracers, num_prog, num_diag)
  enddo

end subroutine tracer_manager_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_tracer_meta_data">
!   <OVERVIEW>
! read tracer table and store tracer information associated with "model"
! in "tracers" array. 
!   </OVERVIEW>
subroutine get_tracer_meta_data(model, num_tracers,num_prog,num_diag)

integer,  intent(in) :: model ! model being used
integer, intent(out) :: num_tracers, num_prog, num_diag
character(len=256)    :: warnmesg

character(len=32)  :: name_type, type, name
integer :: n, m, mod, num_tracer_methods, nfields, swop
integer :: j, log_unit, num_methods
logical :: flag_type
type(method_type), dimension(MAX_TRACER_METHOD) :: methods
integer :: instances, siz_inst,i
character(len = 32) :: digit,suffnam

character(len=128) :: list_name , control
integer            :: index_list_name
logical :: fm_success

!   <ERROR MSG="invalid model type" STATUS="FATAL">
!     The index for the model type is invalid.
!   </ERROR>
if (model .ne. MODEL_ATMOS .and. model .ne. MODEL_LAND .and. &
    model .ne. MODEL_OCEAN .and. model .ne. MODEL_ICE  .and. &
    model .ne. MODEL_COUPLER) call mpp_error(FATAL,'tracer_manager_init : invalid model type')

! One should only call get_tracer_meta_data once for each model type
! Therefore need to set up an array to stop the subroutine being 
! unnecssarily called multiple times.

if ( model_registered(model) ) then
! This routine has already been called for the component model.
! Fill in the values from the previous registration and return.
  num_tracers = total_tracers(model)
  num_prog    = prog_tracers(model)
  num_diag    = diag_tracers(model) 
  return
endif

! Initialize the number of tracers to zero.
num_tracers = 0; num_prog = 0; num_diag = 0

call field_manager_init(nfields=nfields)

!   <ERROR MSG="No tracers are available to be registered." STATUS="NOTE">
!      No tracers are available to be registered. This means that the field
!      table does not exist or is empty.
!   </ERROR>
if (nfields == 0 ) then
if (mpp_pe() == mpp_root_pe()) &
  call mpp_error(NOTE,'tracer_manager_init : No tracers are available to be registered.')
  return
endif

! search through field entries for model tracers
total_tracers(model) = 0

do n=1,nfields
   call get_field_info(n,type,name,mod,num_methods)

   if (mod == model .and. type == 'tracer') then
         num_tracer_fields = num_tracer_fields + 1
         total_tracers(model) = total_tracers(model) + 1
         TRACER_ARRAY(model,total_tracers(model))  = num_tracer_fields
!   <ERROR MSG="MAX_TRACER_FIELDS exceeded" STATUS="FATAL">
!     The maximum number of tracer fields has been exceeded.
!   </ERROR>
         if(num_tracer_fields > MAX_TRACER_FIELDS) call mpp_error(FATAL,'tracer_manager_init: MAX_TRACER_FIELDS exceeded')
         tracers(num_tracer_fields)%model          = model
         tracers(num_tracer_fields)%tracer_name    = name
         tracers(num_tracer_fields)%tracer_units   = 'none'
         tracers(num_tracer_fields)%tracer_longname = tracers(num_tracer_fields)%tracer_name
         tracers(num_tracer_fields)%instances_set   = .FALSE.
         num_tracer_methods     = 0
         methods = default_method ! initialize methods array
         call get_field_methods(n,methods)
         do j=1,num_methods
            select case (methods(j)%method_type) 
            case ('units')
               tracers(num_tracer_fields)%tracer_units   = methods(j)%method_name
            case ('longname')
               tracers(num_tracer_fields)%tracer_longname = methods(j)%method_name
            case ('instances')
!               tracers(num_tracer_fields)%instances = methods(j)%method_name
               siz_inst = parse(methods(j)%method_name,"",instances)
               tracers(num_tracer_fields)%instances = instances
               tracers(num_tracer_fields)%instances_set   = .TRUE.
            case default
               num_tracer_methods = num_tracer_methods+1
!               tracers(num_tracer_fields)%methods(num_tracer_methods) = methods(j)
            end select
         enddo
         tracers(num_tracer_fields)%num_methods = num_tracer_methods
         tracers(num_tracer_fields)%needs_init = .false.
         flag_type = query_method ('tracer_type',model,total_tracers(model),name_type)
         if (flag_type .and. name_type == 'diagnostic') then
            tracers(num_tracer_fields)%is_prognostic = .false.
         else   
            tracers(num_tracer_fields)%is_prognostic = .true.
         endif   
         if (tracers(num_tracer_fields)%is_prognostic) then
            num_prog = num_prog+1
         else
            num_diag = num_diag+1
         endif
   endif
enddo

! Now cycle through the tracers and add additional instances of the tracers.

do n = 1, num_tracer_fields !{
!   call get_field_info(n,type,name,mod,num_methods)

  if ( model == tracers(n)%model .and. tracers(n)%instances_set ) then !{ We have multiple instances of this tracer

    if ( num_tracer_fields + tracers(n)%instances > MAX_TRACER_FIELDS ) then
      write(warnmesg, '("tracer_manager_init: Number of tracers will exceed MAX_TRACER_FIELDS with &
                       &multiple (",I3," instances) setup of tracer ",A)') tracers(n)%instances,tracers(n)%tracer_name
      call mpp_error(FATAL, warnmesg)
    endif                        

    do i = 2, tracers(n)%instances !{
      num_tracer_fields = num_tracer_fields + 1
      total_tracers(model) = total_tracers(model) + 1
      TRACER_ARRAY(model,total_tracers(model))  = num_tracer_fields
      ! Copy the original tracer type to the multiple instances.
      tracers(num_tracer_fields) = tracers(n)
      if ( query_method ('instances', model,model_tracer_number(model,n),name, control)) then !{
          
        if (i .lt. 10) then  !{
           write (suffnam,'(''suffix'',i1)') i
           siz_inst = parse(control, suffnam,digit)
           if (siz_inst == 0 ) then
             write (digit,'(''_'',i1)') i
           else
             digit = "_"//trim(digit)
           endif  
        elseif (i .lt. 100) then  !}{
           write (suffnam,'(''suffix'',i2)') i
           siz_inst = parse(control, suffnam,digit)
           if (siz_inst == 0 ) then
             write (digit,'(''_'',i2)') i
           else
             digit = "_"//trim(digit)
           endif
        else  !}{
          call mpp_error(FATAL, 'tracer_manager_init: MULTIPLE_TRACER_SET_UP exceeds 100 for '//tracers(n)%tracer_name )
        endif  !}

        select case(model)
          case (MODEL_COUPLER)
            list_name = "/coupler_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_ATMOS)
            list_name = "/atmos_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_OCEAN)
            list_name = "/ocean_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_ICE  )
            list_name = "/ice_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case (MODEL_LAND )
            list_name = "/land_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          case default
            list_name = "/default/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
        end select

        if (mpp_pe() == mpp_root_pe() ) write (*,*) "Creating list name = ",trim(list_name)//trim(digit)

        index_list_name = fm_copy_list(trim(list_name),digit, create = .true.)
        tracers(num_tracer_fields)%tracer_name = trim(tracers(num_tracer_fields)%tracer_name)//trim(digit)
      endif !}
         
      if (tracers(num_tracer_fields)%is_prognostic) then !{
         num_prog = num_prog+1
      else !}{
         num_diag = num_diag+1
      endif !}
    enddo !}
    ! Multiple instances of tracers were found so need to rename the original tracer.
    digit = "_1" 
    siz_inst = parse(control, "suffix1",digit)
    if (siz_inst > 0 ) then !{
      digit = "_"//trim(digit)
    endif !}
    fm_success = fm_modify_name(trim(list_name), trim(tracers(n)%tracer_name)//trim(digit))
    tracers(n)%tracer_name = trim(tracers(n)%tracer_name)//trim(digit)
  endif !}
enddo !}

! Find any field entries with the instances keyword.
do n=1,nfields
   call get_field_info(n,type,name,mod,num_methods)

   if ( mod == model .and. type == 'instances' ) then
      call get_field_methods(n,methods)
      do j=1,num_methods

         if (.not.get_tracer_index(mod,methods(j)%method_type,m)) then 
           call mpp_error(FATAL,'tracer_manager_init: The instances keyword was found for undefined tracer '&
           //trim(methods(j)%method_type))
         else
           if ( tracers(m)%instances_set ) &
              call mpp_error(FATAL,'tracer_manager_init: The instances keyword was found for '&
              //trim(methods(j)%method_type)//' but has previously been defined in the tracer entry')
           siz_inst = parse(methods(j)%method_name,"",instances)
           tracers(m)%instances = instances
           call mpp_error(NOTE,'tracer_manager_init: '//trim(instantiations(j)%name)// &
                               ' will have '//trim(methods(j)%method_name)//' instances')
         endif
         if ( num_tracer_fields + instances > MAX_TRACER_FIELDS ) then
           write(warnmesg, '("tracer_manager_init: Number of tracers will exceed MAX_TRACER_FIELDS with &
                       &multiple (",I3," instances) setup of tracer ",A)') tracers(m)%instances,tracers(m)%tracer_name
           call mpp_error(FATAL, warnmesg)
         endif                        
! We have found a valid tracer that has more than one instantiation.
! We need to modify that tracer name to tracer_1 and add extra tracers for the extra instantiations.
         if (instances .eq. 1) then
           siz_inst = parse(methods(j)%method_control, 'suffix1',digit)
           if (siz_inst == 0 ) then
             digit = '_1'
           else
             digit = "_"//trim(digit)
           endif  
         endif
         do i = 2, instances
           num_tracer_fields = num_tracer_fields + 1
           total_tracers(model) = total_tracers(model) + 1
           TRACER_ARRAY(model,total_tracers(model))  = num_tracer_fields
           tracers(num_tracer_fields)                =  tracers(m)
           
           if (i .lt. 10) then  !{
             write (suffnam,'(''suffix'',i1)') i
             siz_inst = parse(methods(j)%method_control, suffnam,digit)
             if (siz_inst == 0 ) then
               write (digit,'(''_'',i1)') i
             else
               digit = "_"//trim(digit)
             endif  
          elseif (i .lt. 100) then  !}{
             write (suffnam,'(''suffix'',i2)') i
             siz_inst = parse(methods(j)%method_control, suffnam,digit)
             if (siz_inst == 0 ) then
               write (digit,'(''_'',i2)') i
             else
               digit = "_"//trim(digit)
             endif
          else  !}{
            call mpp_error(FATAL, 'tracer_manager_init: MULTIPLE_TRACER_SET_UP exceeds 100 for '&
                                  //tracers(num_tracer_fields)%tracer_name )
          endif  !}

          select case(model)
            case (MODEL_COUPLER)
              list_name = "/coupler_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_ATMOS)
              list_name = "/atmos_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_OCEAN)
              list_name = "/ocean_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_ICE  )
              list_name = "/ice_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case (MODEL_LAND )
              list_name = "/land_mod/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
            case default
              list_name = "/default/tracer/"//trim(tracers(num_tracer_fields)%tracer_name)
          end select

          if (mpp_pe() == mpp_root_pe() ) write (*,*) "Creating list name = ",trim(list_name)

          index_list_name = fm_copy_list(trim(list_name),digit, create = .true.)

          tracers(num_tracer_fields)%tracer_name    =  trim(tracers(num_tracer_fields)%tracer_name)//digit
          if (tracers(num_tracer_fields)%is_prognostic) then
            num_prog = num_prog+1
          else
            num_diag = num_diag+1
          endif
        enddo
!Now rename the original tracer to tracer_1 (or if suffix1 present to tracer_'value_of_suffix1')
        siz_inst = parse(methods(j)%method_control, 'suffix1',digit)
        if (siz_inst == 0 ) then
          digit = '_1'
        else
          digit = "_"//trim(digit)
        endif  
        fm_success = fm_modify_name(trim(list_name), trim(tracers(m)%tracer_name)//trim(digit))
        tracers(m)%tracer_name    =  trim(tracers(m)%tracer_name)//trim(digit)
      enddo
   endif
enddo

num_tracers = num_prog + num_diag
! Make the number of tracers available publicly.
total_tracers(model)    = num_tracers
prog_tracers(model)     = num_prog
diag_tracers(model)     = num_diag
model_registered(model) = .TRUE.

! Now sort through the tracer fields and sort them so that the 
! prognostic tracers are first.

do n=1, num_tracers
  if (.not.check_if_prognostic(model,n) .and. n.le.num_prog) then 
  ! This is a diagnostic tracer so find a prognostic tracer to swop with
    do m = n, num_tracers
       if (check_if_prognostic(model,m) .and. .not.check_if_prognostic(model,n)) then
           swop = TRACER_ARRAY(model,n)
           TRACER_ARRAY(model,n) = TRACER_ARRAY(model,m)
           TRACER_ARRAY(model,m) = swop
           cycle
       endif
    enddo
  endif
enddo

do n=1, num_tracer_fields
  call print_tracer_info(model,n)
enddo

log_unit = stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (log_unit,15) trim(MODEL_NAMES(model)),total_tracers(model)
endif

15 format ('Number of tracers in field table for ',A,' model = ',i4)

end subroutine get_tracer_meta_data
!</SUBROUTINE>


function model_tracer_number(model,n)
integer, intent(in) :: model, n
integer model_tracer_number

integer :: i

model_tracer_number = NO_TRACER

do i = 1, MAX_TRACER_FIELDS
  if ( TRACER_ARRAY(model,i) == n ) then
    model_tracer_number = i
    return
  endif
enddo

end function model_tracer_number

!#######################################################################
!
! <SUBROUTINE NAME="register_tracers">

!   <OVERVIEW>
!      It is not necessary to call this routine.
!      It is included only for backward compatability.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine returns the total number of valid tracers,
!     the number of prognostic and diagnostic tracers.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call register_tracers(model, num_tracers,num_prog,num_diag)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="num_tracers" TYPE="integer">
!    The total number of valid tracers within the component model.
!   </OUT>
!   <OUT NAME="num_prog" TYPE="integer">
!     The number of prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_diag" TYPE="integer">
!     The number of diagnostic tracers within the component model.
!   </OUT>
subroutine register_tracers(model, num_tracers, num_prog, num_diag, num_family)
integer, intent(in) :: model
integer, intent(out) :: num_tracers, num_prog, num_diag
integer, intent(out), optional :: num_family

if(.not.module_is_initialized) call tracer_manager_init

call get_number_tracers(model, num_tracers, num_prog, num_diag, num_family)

end subroutine register_tracers
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="get_number_tracers">
!   <OVERVIEW>
!      A routine to return the number of tracers included in a component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine returns the total number of valid tracers,
!     the number of prognostic and diagnostic tracers
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_number_tracers(model, num_tracers,num_prog,num_diag)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="num_tracers" TYPE="integer, optional">
!    The total number of valid tracers within the component model.
!   </OUT>
!   <OUT NAME="num_prog" TYPE="integer, optional">
!     The number of prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_diag" TYPE="integer, optional">
!     The number of diagnostic tracers within the component model.
!   </OUT>
subroutine get_number_tracers(model, num_tracers, num_prog, num_diag, num_family)

integer,  intent(in) :: model
integer, intent(out), optional :: num_tracers, num_prog, num_diag, num_family

if(.not.module_is_initialized) call tracer_manager_init

!   <ERROR MSG="Model number is invalid." STATUS="FATAL">
!     The index of the component model is invalid.
!   </ERROR>
if (model .ne. MODEL_ATMOS .and. model .ne. MODEL_LAND .and. &
    model .ne. MODEL_OCEAN .and. model .ne. MODEL_ICE  .and. &
    model .ne. MODEL_COUPLER)  &
    call mpp_error(FATAL,"get_number_tracers : Model number is invalid.")

if (present(num_tracers)) num_tracers = total_tracers(model)
if (present(num_prog))    num_prog    = prog_tracers(model)
if (present(num_diag))    num_diag    = diag_tracers(model)
if (present(num_family))  num_family  = 0 ! Needed only for backward compatability with lima

end subroutine get_number_tracers
!</SUBROUTINE>


! <SUBROUTINE NAME="get_tracer_indices">

!   <OVERVIEW>
!     Routine to return the component model tracer indices as defined within
!     the tracer manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If several models are being used or redundant tracers have been written to
! the tracer_table, then the indices in the component model and the tracer
! manager may not have a one to one correspondence. Therefore the component
! model needs to know what index to pass to calls to tracer_manager routines in
! order that the correct tracer information be accessed.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_tracer_indices(model, ind, prog_ind, diag_ind)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             all the tracers within the component model.
!   </OUT>
!   <OUT NAME="prog_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="diag_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the diagnostic tracers within the component model.
!   </OUT>
subroutine get_tracer_indices(model, ind, prog_ind, diag_ind, fam_ind)

integer, intent(in) :: model
integer, intent(out), dimension(:), optional :: ind, prog_ind, diag_ind, fam_ind

integer :: i, j, np, nd, n

if(.not.module_is_initialized) call tracer_manager_init

nd=0;np=0;n=0

! Initialize arrays with dummy values
if (PRESENT(ind))      ind      = NO_TRACER
if (PRESENT(prog_ind)) prog_ind = NO_TRACER
if (PRESENT(diag_ind)) diag_ind = NO_TRACER
if (PRESENT(fam_ind))  fam_ind  = NO_TRACER

do i = 1, MAX_TRACER_FIELDS
j = TRACER_ARRAY(model,i)
 if ( j /= NOTRACER) then
   if ( model == tracers(j)%model) then
      if (PRESENT(ind)) then
         n=n+1
!   <ERROR MSG="index array size too small in get_tracer_indices" STATUS="FATAL">
!     The global index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (n > size(ind(:))) call mpp_error(FATAL,'get_tracer_indices : index array size too small in get_tracer_indices')
         ind(n) = i
      endif

      if (tracers(j)%is_prognostic.and.PRESENT(prog_ind)) then
         np=np+1
!   <ERROR MSG="prognostic array size too small in get_tracer_indices" STATUS="FATAL">
!     The prognostic index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if ( np > size( prog_ind(:)))call mpp_error(FATAL,&
                                          'get_tracer_indices : prognostic array size too small in get_tracer_indices')
         prog_ind(np) = i
      else if (.not.tracers(j)%is_prognostic .and. PRESENT(diag_ind)) then
         nd = nd+1
!   <ERROR MSG="diagnostic array size too small in get_tracer_indices" STATUS="FATAL">
!     The diagnostic index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (nd > size(diag_ind(:))) call mpp_error(FATAL,&
                                         'get_tracer_indices : diagnostic array size too small in get_tracer_indices')
         diag_ind(nd) = i
      endif
   endif
 endif
enddo

return
end subroutine get_tracer_indices
!</SUBROUTINE>

!<FUNCTION NAME= "get_tracer_index">
!   <OVERVIEW>
!     Function which returns the number assigned to the tracer name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This is a function which returns the index, as implied within the component model.
!     There are two overloaded interfaces: one of type integer, one logical.
!   </DESCRIPTION>
!   <TEMPLATE>
!     integer: index = get_tracer_index(model, name,        indices, verbose)
!     logical:    if ( get_tracer_index(model, name, index, indices, verbose) ) then
!   </TEMPLATE>
!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     The name of the tracer (as assigned in the field table).
!   </IN>
!   <IN NAME="indices" TYPE="integer, optional" DIM="(:)">
!     An array indices.
!     When present, the returned index will limit the search for the tracer
!     to those tracers whos indices are amoung those in array "indices".
!     This would be useful when it is desired to limit the search to a subset
!     of the tracers. Such a subset might be the diagnostic or prognostic tracers.
!     (Note that subroutine get_tracer_indices returns these subsets)
!   </IN>
!   <IN NAME="verbose" TYPE="logical, optional">
!     A flag to allow the message saying that a tracer with this name has not 
!     been found. This should only be used for debugging purposes.
!   </IN>
!   <OUT NAME="get_tracer_index" TYPE="integer">
!     integer function:
!       The index of the tracer named "name". 
!       If no tracer by that name exists then the returned value is NO_TRACER.
!     logical function:
!       If no tracer by that name exists then the returned value is .false.,
!       otherwise the returned value is .true.
!   </OUT>
function get_tracer_index_integer(model, name, indices, verbose)

integer, intent(in)                         :: model
character(len=*), intent(in)                :: name
integer, intent(in), dimension(:), optional :: indices
logical, intent(in), optional               :: verbose
integer :: get_tracer_index_integer

integer :: i

if(.not.module_is_initialized) call tracer_manager_init

get_tracer_index_integer = NO_TRACER

if (PRESENT(indices)) then
    do i = 1, size(indices(:))
       if (model == tracers(indices(i))%model .and. lowercase(trim(name)) == trim(tracers(indices(i))%tracer_name)) then
           get_tracer_index_integer = i
           exit
       endif
    enddo
else
    do i=1, num_tracer_fields
       if(TRACER_ARRAY(model,i) == NOTRACER) cycle
       if (lowercase(trim(name)) == trim(tracers(TRACER_ARRAY(model,i))%tracer_name)) then
           get_tracer_index_integer = i!TRACER_ARRAY(model,i)
           exit
       endif
    enddo
end if

verbose_local=.FALSE.
if (present(verbose)) verbose_local=verbose

if (verbose_local) then
! <ERROR MSG="tracer with this name not found: X" STATUS="NOTE">
  if (get_tracer_index_integer == NO_TRACER ) then
    call mpp_error(NOTE,'get_tracer_index : tracer with this name not found: '//trim(name))
  endif
! </ERROR>
endif
   
return

end function get_tracer_index_integer

!#######################################################################
function get_tracer_index_logical(model, name, index, indices, verbose)

integer, intent(in)                         :: model
character(len=*), intent(in)                :: name
integer, intent(out)                        :: index
integer, intent(in), dimension(:), optional :: indices
logical, intent(in), optional               :: verbose
logical :: get_tracer_index_logical

index = get_tracer_index_integer(model, name, indices, verbose)
if(index == NO_TRACER) then
  get_tracer_index_logical = .false.
else
  get_tracer_index_logical = .true.
endif

end function get_tracer_index_logical
!</FUNCTION>

!#######################################################################
! <SUBROUTINE NAME="tracer_manager_end" >
!   <OVERVIEW>
!     Routine to write to the log file that the tracer manager is ending.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Routine to write to the log file that the tracer manager is ending.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_manager_end
!   </TEMPLATE>
subroutine tracer_manager_end

integer :: log_unit

log_unit = stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (log_unit,'(/,(a))') 'Exiting tracer_manager, have a nice day ...'
endif

module_is_initialized = .FALSE.

end subroutine tracer_manager_end
!</SUBROUTINE>

!#######################################################################
!
subroutine print_tracer_info(model,n)
!
! Routine to print out the components of the tracer.
! This is useful for informational purposes.
! Used in get_tracer_meta_data.
!
! Arguments:
! INTENT IN
!  i            : index of the tracer that is being printed.
!
integer, intent(in) :: model,n
integer :: i,log_unit

if(.not.module_is_initialized) call tracer_manager_init

if(mpp_pe()==mpp_root_pe() .and. TRACER_ARRAY(model,n)> 0 ) then
  i = TRACER_ARRAY(model,n)
  log_unit = stdlog()
  write(log_unit, *)'----------------------------------------------------'
  write(log_unit, *) 'Contents of tracer entry ', i
  write(log_unit, *) 'Model type and field name'
  write(log_unit, *) 'Model                : ', tracers(i)%model
  write(log_unit, *) 'Field name           : ', trim(tracers(i)%tracer_name)
  write(log_unit, *) 'Tracer units         : ', trim(tracers(i)%tracer_units)
  write(log_unit, *) 'Tracer longname      : ', trim(tracers(i)%tracer_longname)
  write(log_unit, *) 'Tracer is_prognostic : ', tracers(i)%is_prognostic
  write(log_unit, *)'----------------------------------------------------'
endif

900 FORMAT(A,2(1x,E12.6))
901 FORMAT(E12.6,1x,E12.6)


end subroutine print_tracer_info

!#######################################################################
!
! <SUBROUTINE NAME="get_tracer_names" >
!   <OVERVIEW>
!     Routine to find the names associated with a tracer number.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine can return the name, long name and units associated
!     with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_tracer_names(model,n,name,longname, units)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="name" TYPE="character" >
!     Field name associated with tracer number.
!   </OUT>
!   <OUT NAME="longname" TYPE="character, optional" >
!     The long name associated with tracer number.
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional" >
!     The units associated with tracer number.
!   </OUT>

subroutine get_tracer_names(model,n,name,longname, units, err_msg)

integer,          intent(in)  :: model, n
character (len=*),intent(out) :: name
character (len=*), intent(out), optional :: longname, units, err_msg
character (len=128) :: err_msg_local
integer :: n1
character(len=11) :: chn

if(.not.module_is_initialized) call tracer_manager_init

 if (n < 1 .or. n > total_tracers(model)) then
   write(chn, '(i11)') n
   err_msg_local = ' Invalid tracer index.  Model name = '//trim(MODEL_NAMES(model))//',  Index='//trim(chn)
   if(error_handler('get_tracer_names', err_msg_local, err_msg)) return
 endif
 n1 = TRACER_ARRAY(model,n)

name = trim(tracers(n1)%tracer_name)
if (PRESENT(longname)) longname = trim(tracers(n1)%tracer_longname)
if (PRESENT(units))    units    = trim(tracers(n1)%tracer_units)

end subroutine get_tracer_names
!</SUBROUTINE>
!
!#######################################################################
!
! <FUNCTION NAME="get_tracer_name" >
!   <OVERVIEW>
!     Routine to find the names associated with a tracer number.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine can return the name, long name and units associated with a tracer.
!     The return value of get_tracer_name is .false. when a FATAL error condition is
!     detected, otherwise the return value is .true.
!   </DESCRIPTION>
!   <TEMPLATE>
!     if(.not.get_tracer_name(model,n,name,longname, units, err_msg)) call mpp_error(.....
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="name" TYPE="character" >
!     Field name associated with tracer number.
!   </OUT>
!   <OUT NAME="longname" TYPE="character, optional" >
!     The long name associated with tracer number.
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional" >
!     The units associated with tracer number.
!   </OUT>
!   <OUT NAME="err_msg" TYPE="character, optional" >
!     When present:
!       If a FATAL error condition is detected then err_msg will contain an error message
!       and the return value of get_tracer_name will be .false.
!       If no FATAL error is detected err_msg will be filled with space characters and
!       and the return value of get_tracer_name will be .true.
!     When not present:
!       A FATAL error will result in termination inside get_tracer_name without returning.
!       If no FATAL error is detected the return value of get_tracer_name will be .true.
!   </OUT>

function get_tracer_name(model,n,name,longname, units, err_msg)

logical :: get_tracer_name
integer,          intent(in)  :: model, n
character (len=*),intent(out) :: name
character (len=*), intent(out), optional :: longname, units, err_msg
character (len=128) :: err_msg_local
integer :: n1
character(len=11) :: chn

if(.not.module_is_initialized) call tracer_manager_init

 if (n < 1 .or. n > total_tracers(model)) then
   write(chn, '(i11)') n
   err_msg_local = ' Invalid tracer index.  Model name = '//trim(MODEL_NAMES(model))//',  Index='//trim(chn)
   if(error_handler('get_tracer_name', err_msg_local, err_msg)) then
     get_tracer_name = .false.
     return
   endif
 else
   get_tracer_name = .true.
 endif
 n1 = TRACER_ARRAY(model,n)

name = trim(tracers(n1)%tracer_name)
if (PRESENT(longname)) longname = trim(tracers(n1)%tracer_longname)
if (PRESENT(units))    units    = trim(tracers(n1)%tracer_units)

end function get_tracer_name
!</FUNCTION>
!
!#######################################################################
!
!<FUNCTION NAME= "check_if_prognostic">
!   <OVERVIEW>
!    Function to see if a tracer is prognostic or diagnostic.
!   </OVERVIEW>
!   <DESCRIPTION>
!    All tracers are assumed to be prognostic when read in from the field_table
!    However a tracer can be changed to a diagnostic tracer by adding the line
!    "tracer_type","diagnostic"
!    to the tracer description in field_table.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =check_if_prognostic(model, n)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number
!   </IN>
!   <OUT NAME="check_if_prognostic" TYPE="logical">
!     A logical flag set TRUE if the tracer is 
!                        prognostic.
!   </OUT>
function check_if_prognostic(model, n, err_msg)

integer, intent(in) :: model, n
logical             :: check_if_prognostic
character(len=*), intent(out), optional :: err_msg
character(len=128) :: err_msg_local
character(len=11) :: chn

if(.not.module_is_initialized) call tracer_manager_init

if (n < 1 .or. n > total_tracers(model)) then
  write(chn, '(i11)') n
  err_msg_local = ' Invalid tracer index.  Model name = '//trim(MODEL_NAMES(model))//',  Index='//trim(chn)
  check_if_prognostic = .true.
  if(error_handler('check_if_prognostic', err_msg_local, err_msg)) return
endif

!Convert local model index to tracer_manager index

check_if_prognostic = tracers(TRACER_ARRAY(model,n))%is_prognostic

end function check_if_prognostic
!</FUNCTION>
!
!#######################################################################
!
! <SUBROUTINE NAME="set_tracer_profile" >
!   <OVERVIEW>
!     Subroutine to set the tracer field to the wanted profile.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If the profile type is 'fixed' then the tracer field values are set 
! equal to the surface value.
! If the profile type is 'profile' then the top/bottom of model and
! surface values are read and an exponential profile is calculated,
! with the profile being dependent on the number of levels in the
! component model. This should be called from the part of the dynamical
! core where tracer restarts are called in the event that a tracer
! restart file does not exist.
!
!  This can be activated by adding a method to the field_table
! e.g.
!  "profile_type","fixed","surface_value = 1e-12"
!  would return values of surf_value = 1e-12 and a multiplier of 1.0
!  One can use these to initialize the entire field with a value of 1e-12.
!
!  "profile_type","profile","surface_value = 1e-12, top_value = 1e-15"
!   In a 15 layer model this would return values of surf_value = 1e-12 and 
!   multiplier = 0.6309573 i.e 1e-15 = 1e-12*(0.6309573^15)
!   In this case the model should be MODEL_ATMOS as you have a "top" value.
!
!   If you wish to initialize the ocean model, one can use bottom_value instead
!   of top_value.

!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_profile(model, n, tracer)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <INOUT NAME="tracer_array" TYPE="real">
!     The initialized tracer array.
!   </INOUT>

subroutine set_tracer_profile(model, n, tracer, err_msg)

integer,  intent(in)  :: model, n
   real, intent(inout), dimension(:,:,:) :: tracer
character(len=*), intent(out), optional :: err_msg

real    :: surf_value, multiplier
integer :: numlevels, k, n1, flag
real    :: top_value, bottom_value
character(len=80) :: scheme, control,profile_type
character(len=128) :: err_msg_local
character(len=11) :: chn

if(.not.module_is_initialized) call tracer_manager_init

if (n < 1 .or. n > total_tracers(model)) then
  write(chn, '(i11)') n
  err_msg_local = ' Invalid tracer index.  Model name = '//trim(MODEL_NAMES(model))//',  Index='//trim(chn)
  if(error_handler('set_tracer_profile', err_msg_local, err_msg)) return
endif
n1 = TRACER_ARRAY(model,n)

!default values
profile_type  = 'Fixed'
surf_value = 0.0E+00
top_value  = surf_value
bottom_value = surf_value
multiplier = 1.0

tracer = surf_value

if ( query_method ( 'profile_type',model,n,scheme,control)) then
!Change the tracer_number to the tracer_manager version

  if(lowercase(trim(scheme(1:5))).eq.'fixed') then
    profile_type                   = 'Fixed'
    flag =parse(control,'surface_value',surf_value)
    multiplier = 1.0
    tracer = surf_value
  endif

  if(lowercase(trim(scheme(1:7))).eq.'profile') then
    profile_type                   = 'Profile'
    flag=parse(control,'surface_value',surf_value)
    if (surf_value .eq. 0.0) &
      call mpp_error(FATAL,'set_tracer_profile : Cannot have a zero surface value for an exponential profile. Tracer '&
                           //tracers(n1)%tracer_name//" "//control//" "//scheme)
    select case (tracers(n1)%model)
      case (MODEL_ATMOS)
        flag=parse(control,'top_value',top_value)
        if(mpp_pe()==mpp_root_pe() .and. flag == 0) &
           call mpp_error(NOTE,'set_tracer_profile : Parameter top_value needs to be defined for the tracer profile.')
      case (MODEL_OCEAN)
        flag =parse(control,'bottom_value',bottom_value)
        if(mpp_pe() == mpp_root_pe() .and. flag == 0) &
           call mpp_error(NOTE,'set_tracer_profile : Parameter bottom_value needs to be defined for the tracer profile.')
      case default
!   Should there be a NOTE or WARNING message here?
    end select

! If profile type is profile then set the surface value to the input
! value and calculate the vertical multiplier.
! 
! Assume an exponential decay/increase from the surface to the top level
!  C = C0 exp ( -multiplier* level_number)
!  => multiplier = exp [ ln(Ctop/Csurf)/number_of_levels]
!
numlevels = size(tracer,3) -1
    select case (tracers(n1)%model)
      case (MODEL_ATMOS)
        multiplier = exp( log (top_value/surf_value) /numlevels)
        tracer(:,:,1) = surf_value
        do k = 2, size(tracer,3)
          tracer(:,:,k) = tracer(:,:,k-1) * multiplier
        enddo
      case (MODEL_OCEAN)
        multiplier = exp( log (bottom_value/surf_value) /numlevels)
        tracer(:,:,size(tracer,3)) = surf_value
        do k = size(tracer,3) - 1, 1, -1
          tracer(:,:,k) = tracer(:,:,k+1) * multiplier
        enddo
      case default
    end select
  endif !scheme.eq.profile

  if (mpp_pe() == mpp_root_pe() ) write(*,700) 'Tracer ',trim(tracers(n1)%tracer_name),    &
                            ' initialized with surface value of ',surf_value, &
                            ' and vertical multiplier of ',multiplier
  700 FORMAT (3A,E12.6,A,F10.6)

endif ! end of query scheme

end subroutine set_tracer_profile
!</SUBROUTINE>

!
!#######################################################################
!
! <FUNCTION NAME="query_method" >
!   <OVERVIEW>
!     A function to query the "methods" associated with each tracer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to query the "methods" associated with each tracer. The
!  "methods" are the parameters of the component model that can be
!  adjusted by user by placing formatted strings, associated with a
!  particular tracer, within the field table.
!  These methods can control the advection, wet deposition, dry
!  deposition or initial profile of the tracer in question. Any
!  parametrization can use this function as long as a routine for parsing
!  the name and control strings are provided by that routine.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =query_method  (method_type, model, n, name, control)
!   </TEMPLATE>

!   <IN NAME="method_type" TYPE="character">
!     The method that is being requested.
!   </IN>
!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number
!   </IN>
!   <OUT NAME="name" TYPE="character">
!     A string containing the modified name to be used with
!     method_type. i.e. "2nd_order" might be the default for 
!     advection. One could use "4th_order" here to modify 
!     that behaviour.
!   </OUT>
!   <OUT NAME="control" TYPE="character, optional">
!     A string containing the modified parameters that are 
!     associated with the method_type and name.
!   </OUT>
!   <OUT NAME="query_method" TYPE="logical">
!      A flag to show whether method_type exists with regard to
!      tracer n. If method_type is not present then one must
!      have default values.
!   </OUT>

!<NOTE>
!  At present the tracer manager module allows the initialization of a tracer
!  profile if a restart does not exist for that tracer. 
!  Options for this routine are as follows
!
!  Tracer profile setup
!  ==================================================================
!  |method_type  |method_name  |method_control                      |
!  ==================================================================
!  |profile_type |fixed        |surface_value = X                   |
!  |profile_type |profile      |surface_value = X, top_value = Y    |(atmosphere)
!  |profile_type |profile      |surface_value = X, bottom_value = Y |(ocean)
!  ==================================================================
!
!</NOTE>
 function query_method  (method_type, model, n, name, control, err_msg)
!
!  A function to query the schemes associated with each tracer. 
!  
!  INTENT IN
!   method_type  : The method that is being requested.
!   model        : The model that you are calling this function from.
!   n            : The tracer number.
!  INTENT OUT
!   name         : A string containing the modified name to be used with
!                  method_type. i.e. "2nd_order" might be the default for 
!                  advection. One could use "4th_order" here to modify 
!                  that behaviour.
!   control      : A string containing the modified parameters that are 
!                  associated with the method_type and name.
!   query_method : A flag to show whether method_type exists with regard 
!                  to tracer n. If method_type is not present then one
!                  must have default values.

 character(len=*), intent(in)            :: method_type
 integer         , intent(in)            :: model, n
 character(len=*), intent(out)           :: name
 character(len=*), intent(out), optional :: control, err_msg
 logical                                 :: query_method

 integer :: n1
 character(len=256) :: list_name, control_tr
 character(len=11)  :: chn
 character(len=128) :: err_msg_local

 if(.not.module_is_initialized) call tracer_manager_init

!Convert the local model tracer number to the tracer_manager version.

 if (n < 1 .or. n > total_tracers(model)) then
   write(chn, '(i11)') n
   err_msg_local = ' Invalid tracer index.  Model name = '//trim(MODEL_NAMES(model))//',  Index='//trim(chn)
   if(error_handler('query_method', err_msg_local, err_msg)) return
 endif

 n1 = TRACER_ARRAY(model,n)

 select case(model)
  case (MODEL_COUPLER)
   list_name = "/coupler_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
  case (MODEL_ATMOS)
   list_name = "/atmos_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
  case (MODEL_OCEAN)
   list_name = "/ocean_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
  case (MODEL_ICE  )
   list_name = "/ice_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
  case (MODEL_LAND )
   list_name = "/land_mod/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
  case default
   list_name = "/default/tracer/"//trim(tracers(n1)%tracer_name)//"/"//trim(method_type)
 end select

 name = ''
 control_tr = ''
 query_method = fm_query_method(list_name, name, control_tr)

 if ( present(control)) control = trim(control_tr)

 end function query_method
!</FUNCTION>

!<SUBROUTINE NAME="set_tracer_atts">
!   <OVERVIEW>
!     A subroutine to allow the user set the tracer longname and units from the 
!     tracer initialization routine.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to allow the user set the tracer longname and units from the 
!     tracer initialization routine. It seems sensible that the user who is 
!     coding the tracer code will know what units they are working in and it 
!     is probably safer to set the value in the tracer code rather than in 
!     the field table.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_atts(model, name, longname, units)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     Tracer name.
!   </IN>
!   <OUT NAME="longname" TYPE="character, optional">
!     A string describing the longname of the tracer for output to NetCDF files
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional">
!     A string describing the units of the tracer for output to NetCDF files
!   </OUT>
subroutine set_tracer_atts(model, name, longname, units)

integer, intent(in)                    :: model
character(len=*), intent(in)           :: name
character(len=*), intent(in), optional :: longname, units

integer :: n, index
logical :: success
character(len=128) :: list_name

if ( get_tracer_index(model,name,n) ) then
    tracers(TRACER_ARRAY(model,n))%tracer_units   = units
    tracers(TRACER_ARRAY(model,n))%tracer_longname = longname
  select case(model)
    case(MODEL_COUPLER) 
      list_name = "/coupler_mod/tracer/"//trim(name)
    case(MODEL_ATMOS) 
      list_name = "/atmos_mod/tracer/"//trim(name)
    case(MODEL_OCEAN) 
      list_name = "/ocean_mod/tracer/"//trim(name)
    case(MODEL_LAND) 
      list_name = "/land_mod/tracer/"//trim(name)
    case(MODEL_ICE) 
      list_name = "/ice_mod/tracer/"//trim(name)
    case DEFAULT 
      list_name = "/"//trim(name)
  end select      

! Method_type is a list, method_name is a name of a parameter and method_control has the value.
!    list_name = trim(list_name)//"/longname"
  if ( fm_exists(list_name)) then
    success = fm_change_list(list_name)
    if ( present(longname) ) then
      if ( longname .ne. "" ) index = fm_new_value('longname',longname)
    endif
    if ( present(units) ) then
      if (units .ne. "" ) index = fm_new_value('units',units)
    endif
  endif  
    
else
    call mpp_error(NOTE,'set_tracer_atts : Trying to set longname and/or units for non-existent tracer : '//trim(name))
endif

end subroutine set_tracer_atts
!</SUBROUTINE>

!<SUBROUTINE NAME="set_tracer_method">
!   <OVERVIEW> 
!      A subroutine to allow the user to set some tracer specific methods.
!   </OVERVIEW>
!   <DESCRIPTION>
!      A subroutine to allow the user to set methods for a specific tracer. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_method(model, name, method_type, method_name, method_control)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     Tracer name.
!   </IN>
!   <IN NAME="method_type" TYPE="character">
!     The type of the method to be set.
!   </IN>
!   <IN NAME="method_name" TYPE="character">
!     The name of the method to be set.
!   </IN>
!   <IN NAME="method_control" TYPE="character">
!     The control parameters of the method to be set.
!   </IN>
     
subroutine set_tracer_method(model, name, method_type, method_name, method_control)

integer, intent(in)                    :: model
character(len=*), intent(in)           :: name
character(len=*), intent(in)           :: method_type
character(len=*), intent(in)           :: method_name
character(len=*), intent(in)           :: method_control

integer :: n, num_method, index
logical :: success
character(len=128) :: list_name

if ( get_tracer_index(model,name,n) ) then
  tracers(n)%num_methods = tracers(n)%num_methods + 1
  num_method = tracers(n)%num_methods

  select case(model)
    case(MODEL_COUPLER)
      list_name = "/coupler_mod/tracer/"//trim(name)
    case(MODEL_ATMOS)
      list_name = "/atmos_mod/tracer/"//trim(name)
    case(MODEL_OCEAN)
      list_name = "/ocean_mod/tracer/"//trim(name)
    case(MODEL_LAND)
      list_name = "/land_mod/tracer/"//trim(name)
    case(MODEL_ICE)
      list_name = "/ice_mod/tracer/"//trim(name)
    case DEFAULT
      list_name = "/"//trim(name)
  end select      

  if ( method_control .ne. "" ) then
! Method_type is a list, method_name is a name of a parameter and method_control has the value.
    list_name = trim(list_name)//"/"//trim(method_type)
    if ( fm_exists(list_name)) then
      success = fm_change_list(list_name)
      index = fm_new_value(method_type,method_control)
    endif
  else
    call mpp_error(NOTE,'set_tracer_method : Trying to set a method for non-existent tracer : '//trim(name))
  endif
endif

end subroutine set_tracer_method
!</SUBROUTINE>

function error_handler(routine, err_msg_local, err_msg)
logical :: error_handler
character(len=*), intent(in) :: routine, err_msg_local
character(len=*), intent(out), optional :: err_msg

if(present(err_msg)) then
  err_msg = err_msg_local
  error_handler = .true.    
else
  call mpp_error(FATAL,trim(routine)//': '//trim(err_msg_local))
endif

end function error_handler

end module tracer_manager_mod
