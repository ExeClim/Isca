"""
Recursively scans the argument directory for f90 files with `namelist`
declarations.

The default values of those are then found from the surrounding file and
all compiled and written to output `defaults.nml`.

e.g. Running from this directory on the src tree:
````
$ python get_namelist_defaults.py ../../../../src
defaults.nml written
$ head -10 defaults.nml
&dry_convection_nml
    tau = 'UNDEFINED'
/

&fms_io_nml
    max_files_r = 40
    read_data_bug = .false.
    format = 'netcdf'
    show_open_namelist_file_warning = .false.
    time_stamp_restart = .true.
```

"""

import fnmatch
import mmap
import os
import re
import sys

import f90nml

parsers = (int, f90nml.parser.pybool, f90nml.parser.pyfloat, f90nml.parser.pycomplex, f90nml.parser.pystr)

base_dir = sys.argv[1]

namelist_defaults = {}
for root, dirnames, filenames in os.walk(base_dir, followlinks=True):
        for filename in fnmatch.filter(filenames, '*.[Ff]90') :
            #print('Searching file %r' % filename)
            with open(os.path.join(root, filename), 'r+') as f:
                data = mmap.mmap(f.fileno(), 0)
                mo = re.search(br'namelist\s*/\s*(\w+)\s*/\s*([\w\s,\&!]+?[^&])\n', data,  re.MULTILINE| re.IGNORECASE)
                if mo:
                    sector = mo.group(1).decode()
                    nml = namelist_defaults.setdefault(sector, {})
                    #nml['_filename'] = filename

                    # process the list of namelist parameters
                    params = mo.group(2).replace(b'&', b' ').replace(b'\n', b' ') #.replace(b',', b' ')
                    params = params.split(b',')
                    # remove comments and additional whitespace
                    params = [p.partition(b'!')[0].strip() for p in params if p.strip()]
                    for p in params:
                        # search for first reference to parameter in the file
                        # and get the default value

                        #r = br'(\w+).*::\s*' + p + br'\s*=\s*([\w\d\.]+)'
                        if p:
                            r = br'(.*)' + p + br'.*=\s*((?:\(\/[\d\s,+\-.]+\/\))|(?:\'.*\')|(?:\".*\")|[\-+\w\d.]+)'
                            mp = re.search(r, data, re.IGNORECASE)
                            if mp:
                                if '!' in mp.group(1).decode():
                                    continue
                                value =  mp.group(2).decode()
                                for parse in parsers:
                                    try:
                                        value = parse(value)
                                        break
                                    except:
                                        continue
                                nml[p.decode()] = value
                            else:
                                nml[p.decode()] = 'UNDEFINED'


#print(namelist_defaults)
n = f90nml.Namelist()
n.update(namelist_defaults)
n.write('defaults.nml')
print('defaults.nml written')