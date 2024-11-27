import sh
import os
import pdb
from glob import glob

P = os.path.join

class temporary_exp_object(object):
    def __init__(self, basedir, workdir, datadir, exp_name):
        self.basedir = basedir
        self.workdir = workdir
        self.datadir = datadir
        self.expname = exp_name


def create_exp_object(exp_name, data_directory=None):

    if data_directory is None:
        datadir = os.environ['GFDL_DATA']
    else:
        datadir = data_directory

    workdir = os.environ['GFDL_WORK']
    basedir = os.environ['GFDL_BASE']
    expname = '/'+exp_name+'/'

    exp_object = temporary_exp_object(basedir, workdir, datadir, exp_name)
    
    return exp_object


# def keep_only_certain_restart_files(exp_object, max_num_files, interval=12):

#     #        sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

#             #First defines a list of ALL the restart file numbers
#         files_to_remove=list(range(0,max_num_files))

#             #Then defines a list of the ones we want to KEEP
#         files_to_keep  =list(range(0,max_num_files,interval)) 

#             #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
#         for x in files_to_keep:
#                files_to_remove.remove(x) 

#             #Then we remove them.
#         for entry in files_to_remove:
#             try:
#                 sh.rm(P(exp_object.workdir,exp_object.expname,'restarts','res_'+str(entry)+'.cpio'))
# #                 print P(exp_object.workdir,exp_object.expname,'restarts','res_'+str(entry)+'.cpio')

#             except sh.ErrorReturnCode_1:
#                 pass
# #                 print 'Tried to remove some restart files, but number '+str(entry)+' does not exist'
                
def keep_only_certain_restart_files_data_dir(exp_object, max_num_files=None, interval=12):

    #        sh.ls(sh.glob(P(self.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

        if max_num_files is None:
            month_list = glob(P(exp_object.datadir,exp_object.expname, 'restarts')+'/res*.tar.gz')
            if len(month_list)==0:
                return
            else:
                final_month = month_list[-1].split('/res')
            max_num_files = int(final_month[-1].split('.tar.gz')[0])

            #First defines a list of ALL the restart file numbers
        files_to_remove=list(range(0,max_num_files))

            #Then defines a list of the ones we want to KEEP
        files_to_keep  =list(range(0,max_num_files,interval)) 

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
        for x in files_to_keep:
               files_to_remove.remove(x) 

        first_to_be_removed = True
        number_removed = 0
        number_not_removed = 0
            #Then we remove them.
        for entry in files_to_remove[1:-1]:
            try:
                file_to_remove = P(exp_object.datadir,exp_object.expname, 'restarts', 'res%04d.tar.gz' % entry)
                if os.path.isfile(file_to_remove) and first_to_be_removed:
                    first_to_be_removed=False
                    number_not_removed+=1
                    # print('would have removed '+file_to_remove+' but wanted to make sure not to delete the first restart')
                else:
                    sh.rm(file_to_remove)
                    number_removed+=1
                    # print('have removed ' + file_to_remove)
            except sh.ErrorReturnCode_1:
                number_not_removed+=1                
                # print('could not remove ' + file_to_remove)
                pass
        print(P(exp_object.datadir,exp_object.expname), 'number removed '+str(number_removed), 'number not removed '+str(number_not_removed))
                
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
                sh.rm(P(exp_object.datadir,exp_object.expname,'run%04d' % entry,file_name))
                print(('Removed '+P(exp_object.datadir,exp_object.expname,'run%04d' % entry,file_name)))        
            except sh.ErrorReturnCode_1:
                pass
#                 print 'Tried to remove some atmos_daily files, but number '+str(entry)+' does not exist'


            
if __name__=="__main__":

    max_num_files_input = None
    
    # exp_name_list = ['']
    exp_name_list = glob('/disca/share/sit204/data_isca_from_gv5/frierson_post_soc_fix_*/')
    
    
    for exp_name_input in exp_name_list:    
        print('Percentage progress through list:'+str(exp_name_list.index(exp_name_input)/len(exp_name_list)))
        temp_obj = create_exp_object(exp_name_input, data_directory='/disca/share/sit204/data_from_isca_cpu/')
        # keep_only_certain_restart_files(temp_obj, max_num_files_input)
        keep_only_certain_restart_files_data_dir(temp_obj, max_num_files_input)
        # keep_only_certain_daily_data_uninterp(temp_obj, max_num_files_input, file_name = 'fms_moist.x')
#         keep_only_certain_daily_data_uninterp(temp_obj, max_num_files_input)

    