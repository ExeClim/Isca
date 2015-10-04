#!/bin/csh -f
unalias *

set argv = (`getopt Hbehmsf:t: $*`)

#-----------------------------------------------------------------------

 set sep =  
 set format = standard

 set rec = tail
 set help = 0

 set hours   = 0
 set minutes = 0
 set seconds = 0

#-----------------------------------------------------------------------

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -H:
            set help = 1; breaksw
        case -b:
            set rec = 2; breaksw
        case -e:
            set rec = 1; breaksw
        case -h:
            set hours = 1; breaksw
        case -m:
            set hours = 1; set minutes = 1; breaksw
        case -s:
            set hours = 1; set minutes = 1; set seconds = 1; breaksw
        case -f:
            set format = $argv[2]; shift argv; breaksw
        case -t:
            set sep = $argv[2]; shift argv; breaksw
    endsw
    shift argv
end
shift argv

#  --- help output ---

 if ( $help ) then
   cat << END
   time_stamp.csh [ -behms -f format -t separator ]

      -H help (no execution)
      -b beginning date
      -e ending date (default)
      -h hours
      -m hours,minutes
      -s hours,minutes,seconds
      -f format=standard(default),european,digital
      -t separator (default=blank)

END
   exit
 endif

#  --- check format ---

 if ( $format != "standard" &&  \
      $format != "european" && $format != "digital" ) then
    echo ERROR invalid format
    exit (4)
 endif
 
 set hsep = $sep
 if ( $format == "standard" || $format == "european" ) set hsep = h

#-----------------------------------------------------------------------

 if ( -e time_stamp.out ) then
     set time_stamp = `tail -$rec time_stamp.out`

     set month_name = `echo $time_stamp[7] | tr "[A-Z]" "[a-z]"`
     set  month_num = `printf "%.2d" $time_stamp[2]`

#    ---- day can have more than 2 digits ----

     set day_num = `printf "%.2d" $time_stamp[3]`
     if ( $month_name == "day"  && $format == "standard" ) \
        set day_num = `printf "%.4d" $time_stamp[3]`

#    ---- hours,min,sec can have only 2 digits ----

     set  hour_num = `printf "%.2d" $time_stamp[4]`
     set  min_num = `printf "%.2d" $time_stamp[5]`
     set  sec_num = `printf "%.2d" $time_stamp[6]`

#    ---- pad ISO years to 4 digits ----
     if ( $format == "digital" ) then
         set year = `printf %.4d $time_stamp[1]` # will work even if year>9999
     endif

#    ---- create date label ----

     set date_name

     if ( $format == "standard" ) then
        if ( $month_name != "day" ) set date_name = $time_stamp[1]
        set date_name = $date_name$month_name$day_num
     else if ( $format == "european" ) then
        set date_name = $day_num$month_name$time_stamp[1]
     else if ( $format == "digital" ) then
        set date_name = $year$sep$month_num$sep$day_num
     endif

        if ( $hours   ) set date_name = $date_name$hsep$hour_num
        if ( $minutes ) set date_name = $date_name$sep$min_num
        if ( $seconds ) set date_name = $date_name$sep$sec_num

 else
#    --- dummy values ---
     set month_name = "xxx"
     set date_name  = "no_time_stamp"
 endif

     echo $date_name

