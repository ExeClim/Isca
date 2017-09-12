import sh
import os
import pdb

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
        files_to_remove=list(range(0,max_num_files))

            #Then defines a list of the ones we want to KEEP
        files_to_keep  =list(range(0,max_num_files,interval)) 

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
        for x in files_to_keep:
               files_to_remove.remove(x) 

            #Then we remove them.
        for entry in files_to_remove:
            try:
                sh.rm(P(exp_object.workdir,exp_object.expname,'restarts','res_'+str(entry)+'.cpio'))
#                 print P(exp_object.workdir,exp_object.expname,'restarts','res_'+str(entry)+'.cpio')

            except sh.ErrorReturnCode_1:
                pass
#                 print 'Tried to remove some restart files, but number '+str(entry)+' does not exist'
                
def keep_only_certain_restart_files_data_dir(exp_object, max_num_files, interval=12):

    #        sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

            #First defines a list of ALL the restart file numbers
        files_to_remove=list(range(0,max_num_files))

            #Then defines a list of the ones we want to KEEP
        files_to_keep  =list(range(0,max_num_files,interval)) 

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
        for x in files_to_keep:
               files_to_remove.remove(x) 

            #Then we remove them.
        for entry in files_to_remove:
            try:
                sh.rm(P(exp_object.datadir,exp_object.expname,'run%03d' % entry,'INPUT','res'))
#                 print 'would be removing ' + P(exp_object.datadir,exp_object.expname,'run'+str(entry),'INPUT','res')
            except sh.ErrorReturnCode_1:
                pass
#                 print 'Tried to remove some restart files, but number '+str(entry)+' does not exist'                

def keep_only_certain_daily_data_uninterp(exp_object, max_num_files, interval=None, file_name = 'atmos_daily.nc'):

    #        sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

            #First defines a list of ALL the restart file numbers
        files_to_remove=list(range(0,max_num_files))

        if interval is not None:
            #Then defines a list of the ones we want to KEEP
            files_to_keep  =list(range(0,max_num_files,interval)) 
            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
            for x in files_to_keep:
                   files_to_remove.remove(x) 

            #Then we remove them.
        for entry in files_to_remove:           
            try:
                sh.rm(P(exp_object.datadir,exp_object.expname,'run%03d' % entry,file_name))
                print(('Removed '+P(exp_object.datadir,exp_object.expname,'run%03d' % entry,file_name)))        
            except sh.ErrorReturnCode_1:
                pass
#                 print 'Tried to remove some atmos_daily files, but number '+str(entry)+' does not exist'


            
if __name__=="__main__":

    max_num_files_input = 325
    
#     exp_name_list=['simple_continents_post_princeton_qflux_anoms_'+str(x) for x in range(31,32)]

#     exp_name_list=['aquaplanet_qflux_anoms_'+str(x) for x in [12,18,23,32,8]]

#     exp_name_list = ['simple_continents_post_princeton_qflux_control_1','simple_continents_post_princeton_fixed_sst_1', 'simple_continents_post_princeton_qflux_control_nod_1', 'simple_continents_post_princeton_qflux_control_scf_1']
#     
#     exp_name_list = ['annual_mean_ice_princeton_qflux_control_matrix_qflux_2017_code_1', 'annual_mean_ice_post_princeton_fixed_sst_1', 'annual_mean_ice_princeton_fixed_sst_1'] 
# 
#     exp_name_list.extend(['annual_mean_ice_post_princeton_fixed_sst_el_nino_1'])

#     exp_name_list = ['simple_continents_post_princeton_qflux_control_1']

#     exp_name_list = ['annual_mean_ice_princeton_qflux_control_1']#, 'annual_mean_ice_post_princeton_qflux_control_1']

#     exp_name_list = ['annual_mean_ice_post_princeton_fixed_sst_TEST_1', 'annual_mean_ice_princeton_qflux_control_matrix_qflux_1']

#     exp_name_list.extend(['simple_continents_post_princeton_fixed_sst_1'])

#     exp_name_list = ['giant_drag_exp_chai_values_without_dc_bug_latest_1']
#     exp_name_list = ['aquaplanet_qflux_control_1']
    
    exp_name_list = ['giant_drag_exp_chai_values_with_dc_bug_latest_start_to_finish_1', 'giant_drag_exp_chai_values_without_dc_bug_latest_start_to_finish_1']
    
    for exp_name_input in exp_name_list:    
        temp_obj = create_exp_object(exp_name_input)
        keep_only_certain_restart_files(temp_obj, max_num_files_input)
        keep_only_certain_restart_files_data_dir(temp_obj, max_num_files_input)
        keep_only_certain_daily_data_uninterp(temp_obj, max_num_files_input, file_name = 'fms_moist.x')
#         keep_only_certain_daily_data_uninterp(temp_obj, max_num_files_input)

    