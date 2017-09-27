from contextlib import contextmanager
import copy
import logging
import os
import re

import f90nml
from jinja2 import Environment, FileSystemLoader
import sh
import pdb

import create_alert
import getpass

P = os.path.join
_module_directory = os.path.dirname(os.path.realpath(__file__))


mkdir = sh.mkdir.bake('-p')

log = logging.getLogger('gfdl')
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(ch)

def clean_log_info(s):
    if s.strip():
       log.info(s.strip())

def clean_log_error(s):
    if s.strip():
       log.error(s.strip())

def clean_log_debug(s):
    if s.strip():
       log.debug(s.strip())

try:
    GFDL_BASE        = os.environ['GFDL_BASE']
    GFDL_WORK        = os.environ['GFDL_WORK']
    GFDL_DATA        = os.environ['GFDL_DATA']


except Exception, e:
    print('Environment variables GFDL_BASE, GFDL_WORK, GFDL_DATA must be set')
    exit(0)


try:
    _screen_window = os.environ['WINDOW']
    _using_screen = True
except:
    _using_screen = False


def set_screen_title(title):
    if _using_screen:
        try:
            sh.screen('-p', _screen_window, '-X', 'title', title)
            sh.screen('-X', 'redisplay')
        except Exception as e:
            log.warning('Screen title could not be changed.')
            log.debug(e.message)

class CompilationError(Exception):
    pass

class Experiment(object):
    """A basic GFDL experiment"""

    RESOLUTIONS = {
        'T170': {
            'num_lon': 512,
            'num_lat': 256,
            'num_fourier': 170,
            'num_spherical': 171
        },

        'T85': {
            'num_lon': 256,
            'num_lat': 128,
            'num_fourier': 85,
            'num_spherical': 86
        },

        'T42': {
            'num_lon': 128,
            'num_lat': 64,
            'num_fourier': 42,
            'num_spherical': 43,
        },
    }

    def __init__(self, name, commit=None, repo=None, overwrite_data=False,run_idb=False):
        super(Experiment, self).__init__()
        self.name = name


        # set the default locations of working directory,
        # executable directory, restart file storage, and
        # output data directory.
        # These can be overridden e.g. if an experiment is to use
        # the same executable as another
        self.workdir = P(GFDL_WORK, self.name)

        # where executable will be compiled to / fectched from
        if run_idb:
            self.execdir = P(self.workdir, 'exec_debug')  #compiled with debug flags
        else:
            self.execdir = P(self.workdir, 'exec')

        self.restartdir = P(self.workdir, 'restarts') # where restarts will be stored
        self.rundir = P(self.workdir, 'run')          # temporary area an individual run will be performed
        self.datadir = P(GFDL_DATA, self.name)        # where run data will be moved to upon completion

        self.log = log

        if os.path.isdir(self.workdir):
            log.warning('Working directory for exp %r already exists' % self.name)
        else:
            log.debug('Making directory %r' % self.workdir)
            mkdir(self.workdir)

        # Running an experiment from a git repo
        self.repo = repo
        self.commit = commit
        if repo and commit:
            # do a checkout of the specific code and work from that source tree
            self.srcdir = P(self.workdir, 'source')
            self.clone_and_checkout()
        else:
            self.srcdir = GFDL_BASE

        try:
            git_dir = self.srcdir+'/.git'
            commit_id = sh.git("--git-dir="+git_dir, "log", "--pretty=format:'%H'", "-n 1")
            commit_id = str(commit_id).split("'")[1]
        except:
            commit_id = None

        if commit:
            commit_consistency_check = commit_id[0:len(commit)]==commit
            if not commit_consistency_check:
                raise ValueError('commit id specified and commit id actually used are not the same:' +commit+commit_id[0:len(commit)])

        self.commit_id = commit_id

        try:
            git_dir = GFDL_BASE+'/.git'
            commit_id_base = sh.git("--git-dir="+git_dir, "log", "--pretty=format:'%H'", "-n 1")
            commit_id_base = str(commit_id_base).split("'")[1]
            git_diff_output = sh.git("--no-pager", "--git-dir="+git_dir, "--work-tree="+self.srcdir, "diff", "--no-color", self.srcdir)
            git_diff_output = str(git_diff_output).split("\n")            
        except:
            commit_id_base = None
            git_diff_output  = None
            
        self.git_diff_output   = git_diff_output            

        if commit_id==commit_id_base:
            self.commit_id_base = self.commit_id
        else:
            self.commit_id_base = commit_id_base

        if not (repo and commit):
            try:
                git_dir = self.srcdir+'/.git'
                git_status_output = sh.git("--git-dir="+git_dir, "--work-tree="+self.srcdir, "status", "-b", "--porcelain")
                git_status_output = str(git_status_output)
                git_status_output = git_status_output.split("\n")
                git_status_final = ["Running from GFDL_BASE repo, so adding git status output.\n", git_status_output[0]]
                for suffix_to_search_for in ['.f90', '.inc']:
                    git_status_relevant = [str_entry for str_entry in git_status_output if suffix_to_search_for in str_entry.lower()]
                    git_status_final.extend(git_status_relevant)
            except:
                git_status_final = None
        else:
            git_status_final = ['run using specific commit, as specified above, so no git status output for f90 or inc files']
        
        self.git_status_output = git_status_final      

        self.template_dir = P(_module_directory, 'templates')

        self.templates = Environment(loader=FileSystemLoader(self.template_dir))

        self.diag_table = DiagTable()
        self.field_table_file = P(self.template_dir, 'field_table')

        if run_idb:
            self.run_idb = P('true')
        else:
            self.run_idb = P('false')


        # Setup the path_names for compilation
        # 1. read the default path_names file
        # 2. self.path_names can be changed before compilation
        # 3. when self.compile() is called, if more than
        #    `missing_file_tol` paths are not found, abort the compilation.
        #
        # * set `missing_file_tol = -1` to ignore all missing files.
        self.path_names = self._read_pathnames(P(self.template_dir, 'path_names'))
        self.missing_file_tol = 3


        self.namelist_files = [P(self.template_dir, 'core.nml'), P(self.template_dir, 'phys.nml')]
        self.namelist = self.rebuild_namelist()

        self.inputfiles = []

        self.compile_flags=[]

        self.overwrite_data = overwrite_data

        # a short prefix to go before the month and day in
        # the `screen` title. Default rm is for "runmonth".
        self.screen_runmonth_prefix = 'rm'

    def use_template_file(self, filename, values):
        """Use a template file for compilation of the code"""
        new_pathname = P(self.rundir, filename)
        self.templates.get_template(filename).stream(values).dump(new_pathname)
        self.path_names.insert(0, new_pathname)

    def rm_workdir(self):
        try:
            sh.rm(['-r', self.workdir])
        except sh.ErrorReturnCode:
            log.warning('Tried to remove working directory but it doesnt exist')

    def rm_datadir(self):
        try:
            sh.rm(['-r', self.datadir])
        except sh.ErrorReturnCode:
            log.warning('Tried to remove data directory but it doesnt exist')

    def clear_workdir(self):
        self.rm_workdir()
        mkdir(self.workdir)
        log.info('Emptied working directory %r' % self.workdir)

    def keep_only_certain_restart_files(self, max_num_files, interval=12):
        try:
#	    sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

            #First defines a list of ALL the restart file numbers
	    files_to_remove=range(0,max_num_files)

            #Then defines a list of the ones we want to KEEP
	    files_to_keep  =range(0,max_num_files,interval) 

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
	    for x in files_to_keep:
               files_to_remove.remove(x) 

            #Then we remove them.
	    for entry in files_to_remove:
                sh.rm(P(self.workdir,'restarts','res_'+str(entry)+'.cpio'))

        except sh.ErrorReturnCode_1:
            log.warning('Tried to remove some restart files, but the last one doesnt exist')

    def clear_rundir(self):
        sh.cd(self.workdir)
        try:
            sh.rm(['-r', self.rundir])
        except sh.ErrorReturnCode:
            log.warning('Tried to remove run directory but it doesnt exist')
        mkdir(self.rundir)
        log.info('Emptied run directory %r' % self.rundir)

    def _read_pathnames(self, path_names_file):
        with open(path_names_file) as pn:
            return [l.strip() for l in pn]

    def rebuild_namelist(self):
        namelist = f90nml.read(self.namelist_files[0])
        for nl in self.namelist_files[1:]:
            namelist.update(f90nml.read(nl))
        return namelist

    def clone_and_checkout(self, commit=None):
        c = commit or self.commit
        try:
            sh.cd(self.srcdir)
            sh.git.status()
        except Exception as e:
            log.info('Repository not found at %r. Cloning.' % self.srcdir)
            log.debug(e.message)
            try:
                sh.git.clone(self.repo, self.srcdir)
            except Exception as e:
                log.error('Unable to clone repository %r' % self.repo)
                raise e
        try:
            log.info('Checking out commit %r' % c)
            sh.cd(self.srcdir)
            sh.git.checkout(c)
        except Exception as e:
            log.error('Unable to checkout commit %r' % c)
            raise e

    def get_restart_file(self, month):
        return P(self.restartdir, 'res_%d.cpio' % (month))

    def write_path_names(self, outdir):
        log.info('Writing path_names to %r' % P(outdir, 'path_names'))
        with open(P(outdir, 'path_names'), 'w') as pn:
            pn.writelines('\n'.join(self.path_names))

    def write_namelist(self, outdir):
        log.info('Writing namelist to %r' % P(outdir, 'input.nml'))
        self.namelist.write(P(outdir, 'input.nml'))

    def set_resolution(self, res):
        delta = self.RESOLUTIONS[res]
        self.update_namelist({'spectral_dynamics_nml': delta})

    def write_diag_table(self, outdir):
        outfile = P(outdir, 'diag_table')
        log.info('Writing diag_table to %r' % outfile)
        if len(self.diag_table.files):
            template = self.templates.get_template('diag_table')
            calendar = not self.namelist['main_nml']['calendar'].lower().startswith('no_calendar')
            vars = {'calendar': calendar, 'outputfiles': self.diag_table.files.values()}
            template.stream(**vars).dump(outfile)
        else:
            log.error("No output files defined in the DiagTable. Stopping.")
            raise ValueError()

    def write_field_table(self, outdir):
        log.info('Writing field_table to %r' % P(outdir, 'field_table'))
        sh.cp(self.field_table_file, P(outdir, 'field_table'))

    def write_details_file(self, outdir):
        info = self.get_info()
        with open(P(outdir, 'details.txt'), 'w') as f:
            f.writelines(info)

    def get_info(self):
        git_commit = sh.git(['rev-parse', 'HEAD'])
        git_branch = sh.git(['symbolic-ref', 'HEAD'])
        git_desc   = sh.git(['describe', '--always'])
        git_show   = sh.git.show('--pretty')
        return {'git_desc': git_desc,
                'git_show': git_show}

    def disable_rrtm(self):
        # remove all rrtm paths
        self.path_names = [p for p in self.path_names if not 'rrtm' in p]

        # add no compile flag
        self.compile_flags.append('-DRRTM_NO_COMPILE')

        # set the namelist to use gray radiation scheme
        # self.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
        # self.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False

        log.info('RRTM compilation disabled.  Namelist set to gray radiation.')

    def check_path_names(self):
        new_path_names = []
        missing_files = []
        for f in self.path_names:
            if os.path.isfile(os.path.join(self.srcdir, 'src', f)):
                new_path_names.append(f)
            else:
                log.warn('"%r" in path_names but not found in checked out src' % f)
                missing_files.append(f)
        if len(missing_files) >= self.missing_file_tol >= 0:
            log.error('Too many files in "path_names" not found.\nTo prevent this error set exp.missing_file_tol = -1')
            for f in missing_files:
                log.error('File "%r" not found.' % f)
            raise CompilationError('Missing %d compilation files (tolerance %d).' % (len(missing_files), self.missing_file_tol))
        else:
            self.path_names = new_path_names

    def compile(self):
        mkdir(self.execdir)
        vars = {
            'execdir': self.execdir,
            'template_dir': self.template_dir,
            'srcdir': self.srcdir,
            'workdir': self.workdir,
            'compile_flags': ' '.join(self.compile_flags),
            'run_idb': self.run_idb
        }

        self.check_path_names()
        self.write_path_names(self.workdir)

        self.templates.get_template('compile.sh').stream(**vars).dump(P(self.workdir, 'compile.sh'))
        log.info('Running compiler')
        try:
            set_screen_title('compiling')
            sh.bash(P(self.workdir, 'compile.sh'), _out=clean_log_debug, _err=clean_log_debug)
        except sh.ErrorReturnCode as e:
            log.critical('Compilation failed.')
            raise e
        log.debug('Compilation complete.')


    def use_diag_table(self, diag_table):
        """Use a DiagTable object for the diagnostic output specification."""
        self.diag_table = diag_table.copy()


    _day_re = re.compile('Integration completed through\s+([0-9]+) days')
    _month_re = re.compile('Integration completed through\s+([\w\s]+)\s([0-9]+:)')
    def _log_runmonth(self, outputstring):
        s = outputstring.strip()
        m = self._day_re.search(s)
        try:
            days = int(m.group(1))
            set_screen_title('%s:%d:%d' %
                (self.screen_runmonth_prefix, self._cur_month, days))
        except:
            pass
        m = self._month_re.search(s)
        try:
            date = m.group(1)
            set_screen_title('%s:%s' %
                (self.screen_runmonth_prefix, date.strip()))
        except:
            pass
        return clean_log_debug(outputstring)




    def run(self, month, restart_file=None, use_restart=True, num_cores=8, overwrite_data=False, light=False, run_idb=False, experiment_restart=None, email_alerts=True, email_address_for_alerts=None, disk_space_limit=20, disk_space_cutoff_limit=5):

        indir = P(self.rundir, 'INPUT')
        outdir = P(self.datadir, 'run%03d' % month)

        if os.path.isdir(outdir):
            if self.overwrite_data or overwrite_data:
                log.warning('Data for month %d already exists and overwrite_data is True. Overwriting.' % month)
                sh.rm('-r', outdir)
            else:
                log.error('Data for month %d already exists but overwrite_data is False. Stopping.' % month)
                # exit(4)
                return False

        # make the output run folder and copy over the input files
        mkdir([indir, P(self.rundir, 'RESTART'), self.restartdir])

        self.write_namelist(self.rundir)
        self.write_field_table(self.rundir)
        self.write_diag_table(self.rundir)

        for filename in self.inputfiles:
            sh.cp([filename, P(indir, os.path.split(filename)[1])])


        if use_restart:
            if not restart_file:
                # get the restart from previous month
                restart_file = self.get_restart_file(month - 1)
            if not os.path.isfile(restart_file):
                log.error('Restart file not found, expecting file %r' % restart_file)
                exit(2)
            else:
                log.info('Using restart file %r' % restart_file)
        else:
            log.info('Running month %r without restart file' % month)
            restart_file = None

        vars = {
            'month': month,
            'datadir': outdir,
            'rundir': self.rundir,
            'execdir': self.execdir,
            'srcdir': self.srcdir,
            'restart_file': restart_file,
            'num_cores': num_cores,
            'run_idb': run_idb,
            'experiment_restart': experiment_restart
        }

        runmonth = self.templates.get_template('runmonth.sh')

        # employ the template to create a runscript
        t = runmonth.stream(**vars).dump(P(self.rundir, 'runmonth.sh'))

    # Check scratch space has enough disk space
        if email_alerts:
            if email_address_for_alerts is None:
                email_address_for_alerts = getpass.getuser()+'@exeter.ac.uk'
            create_alert.run_alerts(self.execdir, GFDL_BASE, self.name, month, email_address_for_alerts, disk_space_limit, disk_space_cutoff_limit)

        log.info("Running GFDL for month %r" % month)
        self._cur_month = month
        try:
            set_screen_title('month:%d' % month)
            proc = sh.bash(P(self.rundir, 'runmonth.sh'), _out=self._log_runmonth, _err=self._log_runmonth, _bg=True)
            proc.wait()  # allow the background thread to complete
            log.info("Run for month %r complete" % month)
        except KeyboardInterrupt:
            log.error("Manual interrupt, killing process.")
            proc.process.kill()
            log.info("Cleaning run directory.")
            self.clear_rundir()
            return False
        except sh.ErrorReturnCode as e:
            log.error("Run failed for month %r" % month)
            proc.process.kill()
            raise e

        mkdir(outdir)

        #restart_file = P(self.restartdir, 'res_%d.tar.gz' % month)
        #sh.tar('zcvf', restart_file, 'RESTART')
        restart_file = P(self.restartdir, 'res_%d.cpio' % month)
        sh.cd(P(self.rundir, 'RESTART'))
        state_files = sh.glob('*.res*')
        sh.cpio('-ov', _in='\n'.join(state_files), _out=restart_file)
        log.info("Saved restart file %r" % restart_file)
        sh.rm('-r', P(self.rundir, 'RESTART'))

        try:
            git_hash_file = open(P(outdir, 'git_hash_used.txt'), "w")
            if self.commit_id!=self.commit_id_base:
                git_hash_file.write("*---hash of specified commit used for fortran code in workdir---*:\n")
                git_hash_file.write(self.commit_id)
                git_hash_file.write("\n")                 
                git_hash_file.write("\n*---hash of commit used for code in GFDL_BASE, including this python script---*:\n")
                git_hash_file.write(self.commit_id_base)
                git_hash_file.write("\n")            
            else:
                git_hash_file.write("*---hash of commit used for code in GFDL_BASE, including this python script---*:\n")
                git_hash_file.write(self.commit_id)
                git_hash_file.write("\n")
            if self.git_status_output is not None:
                git_hash_file.write("\n"+"*---git status output (only f90 and inc files)---*:\n")
                git_hash_file.writelines( line_out+"\n" for line_out in self.git_status_output)
            if self.git_diff_output is not None:
                git_hash_file.write("\n"+"*---git diff output run in GFDL_BASE (everything)---*:\n")            
                git_hash_file.writelines( line_out+"\n" for line_out in self.git_diff_output)                                
            git_hash_file.close()
        except:
            log.info("Could not output git commit hash")        
        
        if light:
            os.system("cp -a "+self.rundir+"/*.nc "+outdir)
            sh.cp(['-a', P(self.restartdir, 'res_%d.cpio' % (month)), outdir])
            if month > 1:
                try:
                    sh.rm( P(self.restartdir, 'res_%d.cpio' % (month-1)))
                except sh.ErrorReturnCode:
                    log.warning('Previous months restart already removed')

                try:
                    sh.rm( P(self.datadir, 'run03%d' % (month-1) , 'res_%d.cpio' % (month-1)))
                except sh.ErrorReturnCode:
                    log.warning('Previous months restart already removed')

        else:
            sh.cp(['-a', self.rundir+'/.', outdir])
        self.clear_rundir()
        sh.cd(self.rundir)
        return True

    runmonth = run

    def derive(self, new_experiment_name):
        """Derive a new experiment based on this one."""
        new_exp = Experiment(new_experiment_name)
        new_exp.execdir = self.execdir
        new_exp.namelist = self.namelist.copy()
        new_exp.use_diag_table(self.diag_table.copy())
        new_exp.inputfiles = self.inputfiles
        new_exp.commit_id = self.commit_id
        new_exp.commit_id_base = self.commit_id_base        
        new_exp.git_status_output = self.git_status_output
        new_exp.git_diff_output   = self.git_diff_output        
        return new_exp

    def run_parameter_sweep(self, parameter_values, runs=10, num_cores=16):
        # parameter_values should be a namelist fragment, with multiple values
        # for each study e.g. to vary obliquity:
        # exp.run_parameter_sweep({'astronomy_nml': {'obliq': [0.0, 5.0, 10.0, 15.0]}})
        # will run 4 independent studies and create data e.g.
        # <exp_name>/astronomy_nml_obliq_<0.0, 5.0 ...>/run[1-10]/daily.nc
        params = [(sec, name, values) for sec, parameters in parameter_values.items()
                        for name, values in parameters.items()]
        # make a list of lists of namelist section, parameter names and values
        params = [[(a,b,val) for val in values] for a,b,values in params]
        parameter_space = itertools.product(params)
        for combo in parameter_space:
            title = '_'.join(['%s_%s_%r' % (sec[:3], name[:5], val) for sec, name, val in combo])
            exp = self.derive(self.name + '_' + title)
            for sec, name, val in combo:
                exp.namelist[sec][name] = val
            exp.clear_rundir()
            exp.run(1, use_restart=False, num_cores=num_cores)
            for i in range(runs-1):
                exp.run(i+2)

    def update_namelist(self, new_vals):
        """Update the namelist sections, overwriting existing values."""
        for sec in new_vals:
            if sec not in self.namelist:
                self.namelist[sec] = {}
            nml = self.namelist[sec]
            nml.update(new_vals[sec])

    def runinterp(self, month, infile, outfile, var_names = '-a', p_levs = "EVEN", rm_input=False):
        """Interpolate data from sigma to pressure levels. Includes option to remove original file."""
        import subprocess
        pprocess = P(GFDL_BASE,'postprocessing/plevel_interpolation/scripts')
        interper = 'source '+pprocess+'/plevel.sh -i '
        inputfile = P(self.datadir, 'run%03d' % month, infile)
        outputfile = P(self.datadir, 'run%03d' % month, outfile)
        
        # Select from pre-chosen pressure levels, or input new ones in hPa in the format below.
        if p_levs == "MODEL":
            plev = ' -p "2 9 18 38 71 125 206 319 471 665 904 1193 1532 1925 2375 2886 3464 4115 4850 5679 6615 7675 8877 10244 11801 13577 15607 17928 20585 23630 27119 31121 35711 40976 47016 53946 61898 71022 81491 93503" '
        elif p_levs == "EVEN":
            plev = ' -p "100000 95000 90000 85000 80000 75000 70000 65000 60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 10000 5000" '
        else:
            plev = p_levs
        command = interper + inputfile + ' -o ' + outputfile + plev + var_names
        subprocess.call([command], shell=True)
        if rm_input:
            sh.rm( inputfile)


class DiagTable(object):
    def __init__(self):
        super(DiagTable, self).__init__()
        self.files = {}

    def add_file(self, name, freq, units="hours", time_units=None):
        if time_units is None:
            time_units = units
        self.files[name] = {
            'name': name,
            'freq': freq,
            'units': units,
            'time_units': time_units,
            'fields': []
        }

    def add_field(self, module, name, time_avg=False, files=None):
        if files is None:
            files = self.files.keys()

        for fname in files:
            self.files[fname]['fields'].append({
                'module': module,
                'name': name,
                'time_avg': time_avg
                })

    def copy(self):
        d = DiagTable()
        d.files = copy.deepcopy(self.files)
        return d
