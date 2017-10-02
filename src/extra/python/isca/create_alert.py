import send_email as send
from check_disk_space import disk_usage
import os
import pdb

def run_alerts(execdir, basedir, exp_name, month, recipient_email_address, disk_space_limit, disk_space_cutoff_limit):

    disk_space_alert(execdir, basedir, exp_name, month, recipient_email_address, disk_space_limit, disk_space_cutoff_limit)


def disk_space_alert(dir, basedir, exp_name, month, recipient_email_address,limit, cutoff_limit):

    remaining_space_in_dir=disk_usage(dir)
    free_space_in_gb=remaining_space_in_dir.free/1e9
    
    if free_space_in_gb < limit and free_space_in_gb > cutoff_limit:
        alert_message="Disk space less than "+str(limit)+"Gb before running month number "+str(month)+" in experiment "+exp_name
        print(alert_message+", sending email")
        send.send_email_fn(recipient_email_address, alert_message, basedir)
    elif free_space_in_gb < limit and free_space_in_gb < cutoff_limit:
        alert_message="Disk space less than cutoff limit of "+str(cutoff_limit)+"Gb before running month number "+str(month)+" in experiment "+exp_name+", therefore run will be killed."
        send.send_email_fn(recipient_email_address, alert_message, basedir)    
        raise IOError(alert_message)
    else:
        print('Disk space more than ' + str(limit) + 'Gb - not sending alert email.')
        
        
if __name__ == '__main__':

    dir=os.getcwd()
    recipient_email_address="sit204@exeter.ac.uk"
    
    try:
        basedir = os.environ['GFDL_BASE']
    except:
        basedir = './../../../../'
    
    run_alerts(dir, basedir, 'test', 1, recipient_email_address, 2000, 5)
