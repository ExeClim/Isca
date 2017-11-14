import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import socket
import datetime
import pdb

from isca import GFDL_BASE

def get_paz():

    F = open(GFDL_BASE+'/src/extra/python/gfdl/'+'mima_pz.txt','r')
    code = F.read()
    code = code.translate(None, '\n')

    return code

def send_email_fn(to_email,alert_message):

    machine_name=socket.gethostname()
    current_time = datetime.datetime.now().isoformat()

    from_email="mima.python.alerts@gmail.com"

    msg = MIMEMultipart()
    msg['From'] = from_email
    msg['To'] = to_email
    msg['Subject'] = "[Mima-alert] "+alert_message+" on "+machine_name+" at time " + current_time

    body = "This is an automated message. \n"+alert_message
    msg.attach(MIMEText(body, 'plain'))

    try:
        server = smtplib.SMTP('localhost')
    except:
        server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
        server.ehlo()

        try:
            code = get_paz()
        except IOError as error_msg:
            print('Password file is missing - email will not send. Error message: '+ error_msg.strerror+': '+error_msg.filename)
            raise
        server.login(from_email, code)

    try:
        text = msg.as_string()
        server.sendmail(from_email, to_email, text)
        server.quit()
    except:
        print('something went wrong with sending the email')
