import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText 
import socket
import datetime

 
def send_email_fn(to_email,alert_message):

	machine_name=socket.gethostname()
	current_time = datetime.datetime.now().isoformat()

	from_email="mima.python.alerts@gmail.com"

	msg = MIMEMultipart()
	msg['From'] = from_email
	msg['To'] = to_email
	msg['Subject'] = "[Mima-alert] "+alert_message+" on "+machine_name+" at time " + current_time
 
	body = "This is an automated message."
	msg.attach(MIMEText(body, 'plain'))

	server = smtplib.SMTP('smtp.gmail.com', 587)
	server.starttls()
	server.login(from_email, "mima101_win")

	text = msg.as_string() 
	server.sendmail(from_email, to_email, text)
	server.quit()