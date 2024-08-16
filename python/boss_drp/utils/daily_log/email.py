from boss_drp.utils.daily_log.html import daily_log_html
from boss_drp.utils.daily_log.file import daily_log_to_file

import numpy as np
import os.path as ptt
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from os import getenv

def daily_log_email(subject, attachment, logger, obs, mjd,
                    email_file = None, topdir=None, epoch=False,
                    run2d=None, run1d=None, content=None, custom=None,
                    from_domain="chpc.utah.edu",  redux = None):
  
    
    body, _ = daily_log_html(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
                         redux=redux, email=True, epoch=epoch, custom=custom)
    
    daily_log_to_file(obs, mjd, topdir=topdir, run2d=run2d, run1d=run1d,
                     redux=redux, epoch=epoch, custom=custom)#, html_log=body)
    
    try:
        emails = open(email_file).read().splitlines()
    except:
        emails = []
        logger.info(email_file+' does not exist')

    msg = MIMEMultipart("alternative")
    msg['Subject'] = subject
    msg['From'] = f"BOSS Pipeline <{getenv('USER')}@{from_domain}>"
    msg['BCC'] = ', '.join(emails)



    part1 = MIMEText("An html compatible email view is required to view the full email","plain")
    part2 = MIMEText("<br>\n".join(body), "html")
    msg.attach(part1)
    msg.attach(part2)

    if attachment is not None:
        attachment = np.atleast_1d(attachment)
        
        
        for f in attachment:
            with open(f, "rb") as fil:
                part = MIMEApplication(
                    fil.read(),
                    Name=ptt.basename(f)
                )
            # After the file is closed
            part['Content-Disposition'] = 'attachment; filename="%s"' % ptt.basename(f)
            msg.attach(part)
        
    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()
    return(None)
