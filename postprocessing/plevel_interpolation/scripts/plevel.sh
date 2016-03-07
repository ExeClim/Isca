#!/bin/sh
#
#  Script for converting to pressure levels
#-----------------------------------------------------------------------

fields=
more_output=1
ofile=plevel.nc
ifiles=
do_all_3d_fields=.false.
do_all_fields=.false.
allow_zero=.false.
default_missval=.false.
mask_extrap=.true.
tlist=

# default is ncep reanalysis levels
plevs="100000 92500 85000 70000 60000 50000 40000 30000 25000 \
        20000 15000 10000  7000  5000  3000  2000  1000"

#-----------------------------------------------------------------------

while getopts 03amxi:o:d:p:t: arg
do
    case $arg in
       0) allow_zero=.true.;;
       3) do_all_3d_fields=.true.;;
       a) do_all_fields=.true.;;
       m) default_missval=.true.;;
       x) mask_extrap=.false.;;
       i) ifiles=$OPTARG;;
       o) ofile=$OPTARG;;
       d) more_output=$OPTARG;;
       p) plevs=$OPTARG;;
       t) tlist=$OPTARG;;
       *) exit 1;;
    esac
done
shift `expr $OPTIND - 1`
fields=$@

#-----------------------------------------------------------------------

if [ "${ifiles:-NULL}" = "NULL" ]; then
name=`basename $0`
cat << EOF

Interpolates data from model levels to pressure levels.
The input model grid is a hybrid sigma-pressure coordinate
and the output pressure levels may be specified.
The minimum required input fields are "bk", "pk", and "ps".

Usage:  $name [-a] [-3] [-0] [-f] [-d #] -i file [-o ofile] [-m] [fields.....]

        -a        = Output all fields converting 3d fields to pressure levels.
        -3        = Output and convert all 3d fields to pressure levels.
        -0        = When fields sphum or zsurf do not exist use zero, otherwise fail.
        -m        = Default missing value is used for all fields (the _FillValue).
        -x        = DO NOT set data extrapolated beneath the surface to missing values.
        -i file   = Input netcdf file, the file must contain the required variables 
                    (pk,bk,ps,...).
                    If the -i files option is omitted a usage message is printed.
        -o ofile  = The output file name. (Default: plevel.nc)
        -p plevs  = A list of output pressure levels in pascals (with no decimal point).
                    The default is the 17 NCEP reananalysis levels (bottom to top).
                    The list must be in quotes and values must be separated by a space.
        -d value    The verbosity level, use an integer number where value >= 0.
                    (Default: value=1)
        -t #,#,#    The starting, ending, and increment index for time axis processing
                    where # is a positive number.
                    (The default is to process all time indices.)

        fields    = A list of (additional) output fields. If this list is not supplied,
                    then the "-a" or "-3" option must be specified.  Possible list entries
                    included any fields in the input files, and additional fields: slp, hght.
                    Additional input fields may be required for these output fields.

Example:  $name -a -i atmos.nc slp hght

EOF
exit 1
fi

# location of executable

source /opt/modules/default/init/sh
source /home/fms/local/ia64/v12/moduleloads.sh
executable=$FMSLOCALPATH/PLEV.exe
if [ ! -x "$executable" ]; then
   echo "ERROR: executable does not exist"
   echo "       executable=$executable"
   exit 1
fi

#  make sure input files are present

list=
for file in $ifiles;
do
   if [ ! -e $file ]; then
       list="$list $file"
   fi
done
if [ ${#list} -gt 0 ]; then
    echo ERROR: the following input files do not exist
    for file in $list;
    do
      echo $file
    done
    exit 1
fi

# process time loop limits (create array of length 3)

tlist=`echo ${tlist} | sed -e "s/,,,/,,0,/"`
tlist=`echo ${tlist} | sed -e "s/^,/0,/"`
tlist=`echo ${tlist} | sed -e "s/,,/,0,/"`
tloop=`echo ${tlist} | sed -e "s/,/ /g"`

#-----------------------------------------------------------------------
#   ---- namelist for pressure interp program ----
    
    namelist="plev.input.nml"

    echo " &input" > $namelist

# input file names
    echo "     in_file_name =  '$file' ," >> $namelist
    echo "     out_file_name =   '$ofile' ,"  >> $namelist

# input field names
i=0
for field in $fields;
do
    let "i=$i+1"
    echo "    field_names($i) =  '$field' ," >> $namelist
done

# pressure level values
i=0
for prs in $plevs;
do
    let "i=$i+1"
    echo "    pout($i) =  $prs.," >> $namelist
done

# more namelist values
cat >> $namelist << EOF
    do_all_3d_fields = $do_all_3d_fields,
    do_all_fields = $do_all_fields,
    allow_zero_sphum = $allow_zero,
    allow_zero_topog = $allow_zero,
    mask_extrap = $mask_extrap,
    use_default_missing_value = $default_missval,
    verbose = $more_output,
EOF

#--- time loop limits ---
i=0
for index in $tloop;
do
    let "i=$i+1"
    if [ $index -gt 0 ]; then 
       case $i in
          1) echo "    time_beg = $index" >> $namelist;;
          2) echo "    time_end = $index" >> $namelist;;
          3) echo "    time_inc = $index" >> $namelist;;
       esac    
    fi
done

   echo " /" >> $namelist

$executable

#rm -f $namelist

