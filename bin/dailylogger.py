#!/usr/bin/env python3
import logging
import collections
from os import popen, getenv
import os.path as ptt
import smtplib
from email.message import EmailMessage


class Formatter(logging.Formatter):
    def __init__(self):
        super().__init__(fmt="%(levelno)d: %(msg)s", datefmt=None, style='%')
    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        elif record.levelno == logging.DEBUG:
            self._style._fmt = '%(funcName)s: %(message)s'
        else:
            self._style._fmt = "%(levelname)s: %(message)s"
        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)
    
        # Restore the original format configured by the user
        self._style._fmt = format_orig
        return result

class emailLogHandler(logging.Handler):
    def __init__(self, log_queue):
        logging.Handler.__init__(self)
        self.log_queue = log_queue
    def emit(self, record):
        self.log_queue.append(self.format(record))


#def send_email(subject, email_file, attachment, logger):
#
#    try:
#        emails = open(email_file).read().splitlines()
#    except:
#        emails = []
#        logger.info(email_file+' does not exist')
#    for email_add in emails:
#        if len(email_add) == 0:
#            continue
#        try:
#            #cmd = 'echo "'+log+'"| mail -s "'+subject+'" '+email_add+' -A '+attachment
#            cmd = 'echo "" |  mail -v -s \"'+subject+'\" -A '+attachment+' '+email_add
#            stream = popen(cmd)
#            output = stream.read()
#            logger.info(output)
#        except:
#            log = ['ERROR Building Email Log']
#            logger.info(log[0])
#            cmd = 'echo "'+log+'"| mail -v -s "'+subject+'" '+email_add#+' -A '+attachment
#            #cmd = 'echo mail -s "'+subject+'" '+email_add+' -A '+attachment
#            stream = popen(cmd)
#            output = stream.read()
#            logger.info(output)
#    return(None)


def send_email(subject, email_file, attachment, logger, content=None, from_domain="chpc.utah.edu"):


    try:
        emails = open(email_file).read().splitlines()
    except:
        emails = []
        logger.info(email_file+' does not exist')
        
    emails = ' '.join(emails).split()
    msg = EmailMessage()
    if content is None:
        content = subject
    msg.set_content(content)
    msg['Subject'] = subject
    msg['From'] = f"BOSS Pipeline <{getenv('USER')}@{from_domain}>"
    msg['BCC'] = ', '.join(emails)
    if attachment is not None:
        msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
        with open(attachment, 'rb') as fp:
            logdata = fp.read()
            msg.add_attachment(logdata, maintype='text', subtype='plain', filename=ptt.basename(attachment))
    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()
    return(None)


class emailLogger(object):
    def __init__(self, maxlen=None):
        self._log_queue = collections.deque(maxlen=maxlen)
        self._log_handler = emailLogHandler(self._log_queue)
    def contents(self):
        return '\n'.join(self._log_queue)
    @property
    def log_handler(self):
        return self._log_handler

    def send(self, subject, email_file,log):
        try:
            emails = open(email_file).read().splitlines()
        except:
            log.error(email_file+' does not exist')
            emails = []
            
        try:
            send_email(subject, email_file, None, log, content=self.contents(), from_domain="chpc.utah.edu")

#            emails = ' '.join(emails).split()
#            msg = EmailMessage()
#            msg.set_content(self.contents())
#            msg['Subject'] = subject
#            msg['From'] = f"BOSS Pipeline <{getenv('USER')}@{from_domain}>"
#            msg['To'] = ', '.join([emails])
#            s = smtplib.SMTP('localhost')
#            s.send_message(msg)
#            s.quit()
        except:
            outputs = []
            for line in self.contents:
                if 'slurm.session.Client:' in line:
                    continue
                if 'slurm.session.Client: task #' in line:
                    continue
                outputs.append(line)
            self.contents = outputs
            try:
                send_email(subject, email_file, None, log, content=self.contents(), from_domain="chpc.utah.edu")
            except:
                self.contents = ['ERROR Building Email Log']
                send_email(subject, email_file, None, log, content=self.contents(), from_domain="chpc.utah.edu")
#
#        for email_add in emails:
#            if len(email_add) == 0:
#                continue
#            try:
#                cmd = 'echo "'+self.contents()+'"| mail -v -s "'+subject+'" '+email_add
#                stream = popen(cmd)
#                output = stream.read()
#                output
#            except:
#                outputs = []
#                for line in self.contents:
#                    if 'slurm.session.Client:' in line:
#                        continue
#                    if 'slurm.session.Client: task #' in line:
#                        continue
#                    outputs.append(line)
#                self.contents = outputs
#                cmd = 'echo "'+self.contents()+'"| mail -v -s "'+subject+'" '+email_add
#                try:
#                    stream = popen(cmd)
#                    output = stream.read()
#                    output
#                except:
#                    self.contents = ['ERROR Building Email Log']
#                    cmd = 'echo "'+self.contents()+'"| mail -v -s "'+subject+'" '+email_add
#                    stream = popen(cmd)
#                    output = stream.read()
#                    output
        return
