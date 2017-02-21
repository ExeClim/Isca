import sh
import os

P = os.path.join

class temporary_exp_object(object):
    def __init__(self, basedir, workdir, datadir, exp_name):
        self.basedir = basedir
        self.workdir = workdir
        self.datadir = datadir
        self.expname = exp_name


def create_exp_object(exp_name):

    basedir = os.environ['GFDL_BASE']
    workdir = os.environ['GFDL_WORK']
    datadir = os.environ['GFDL_DATA']
    expname = '/'+exp_name+'/'

    exp_object = temporary_exp_object(basedir, workdir, datadir, exp_name)
    
    return exp_object


def keep_only_certain_restart_files(exp_object, max_num_files, interval=12):

    #        sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

            #First defines a list of ALL the restart file numbers
        files_to_remove=range(0,max_num_files)

            #Then defines a list of the ones we want to KEEP
        files_to_keep  =range(0,max_num_files,interval) 

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
        for x in files_to_keep:
               files_to_remove.remove(x) 

            #Then we remove them.
        for entry in files_to_remove:
            try:
                sh.rm(P(exp_object.workdir,exp_object.expname,'restarts','res_'+str(entry)+'.cpio'))
            except sh.ErrorReturnCode_1:
                print 'Tried to remove some restart files, but number '+str(entry)+' does not exist'
                
def keep_only_certain_restart_files_data_dir(exp_object, max_num_files, interval=12):

    #        sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

            #First defines a list of ALL the restart file numbers
        files_to_remove=range(0,max_num_files)

            #Then defines a list of the ones we want to KEEP
        files_to_keep  =range(0,max_num_files,interval) 

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
        for x in files_to_keep:
               files_to_remove.remove(x) 

            #Then we remove them.
        for entry in files_to_remove:
            try:
                sh.rm(P(exp_object.datadir,exp_object.expname,'run'+str(entry),'INPUT','res'))
#                 print 'would be removing ' + P(exp_object.datadir,exp_object.expname,'run'+str(entry),'INPUT','res')
            except sh.ErrorReturnCode_1:
                print 'Tried to remove some restart files, but number '+str(entry)+' does not exist'                

            
if __name__=="__main__":

    max_num_files_input = 1200
    
#     exp_name_list=['simple_continents_post_princeton_qflux_anoms_'+str(x) for x in range(1,21)]

#     exp_name_list=['annual_mean_ice_post_princeton_qflux_anoms_'+str(x) for x in range(1,21)]

#     exp_name_list = ['simple_continents_post_princeton_qflux_control_1','simple_continents_post_princeton_fixed_sst_1', 'simple_continents_post_princeton_qflux_control_nod_1', 'simple_continents_post_princeton_qflux_control_scf_1']
#     
#     exp_name_list.extend(['annual_mean_ice_princeton_qflux_control_matrix_qflux_2017_code_1', 'annual_mean_ice_post_princeton_fixed_sst_1', 'annual_mean_ice_post_princeton_qflux_control_1'])
# 
#     exp_name_list.extend(['annual_mean_ice_post_princeton_fixed_sst_el_nino_1'])

#     exp_name_list = ['giant_drag_exp_chai_values_1_bar_damping_without_dc_bug_latest_1', 'giant_drag_exp_chai_values_1_bar_damping_without_dc_bug_latest_2']

    exp_name_list = ['annual_mean_ice_princeton_fixed_sst_1', 'annual_mean_ice_princeton_qflux_control_1']


    for exp_name_input in exp_name_list:    
        temp_obj = create_exp_object(exp_name_input)
#         keep_only_certain_restart_files(temp_obj, max_num_files_input)
        keep_only_certain_restart_files_data_dir(temp_obj, max_num_files_input)

    