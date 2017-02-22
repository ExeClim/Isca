import send_email as send
from check_disk_space import disk_usage
import os

def run_alerts(dir,recipient_email_address, disk_space_limit):

    disk_space_alert(dir,recipient_email_address, disk_space_limit)


def disk_space_alert(dir,recipient_email_address,limit):

    remaining_space_in_dir=disk_usage(dir)
    free_space_in_gb=remaining_space_in_dir.free/1e9
    
    if free_space_in_gb < limit:
        alert_message="Disk space less than "+str(limit)+"Gb on machine " + os.uname()[1]
        print alert_message+", sending email"
        send.send_email_fn(recipient_email_address, alert_message)
    else:
        print 'Disk space more than ' + str(limit) + 'Gb - not sending alert email.'
        
        
if __name__ == '__main__':

    dir=os.getcwd()
    recipient_email_address="sit204@exeter.ac.uk"
    
    run_alerts(dir,recipient_email_address)
