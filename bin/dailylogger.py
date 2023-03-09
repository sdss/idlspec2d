#!/usr/bin/env python3
import logging
import collections
from os import popen

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
        for email_add in emails:
            if len(email_add) == 0:
                continue
            cmd = 'echo "'+self.contents()+'"| mail -s "'+subject+'" '+email_add
            stream = popen(cmd)
            output = stream.read()
            output
                
        return
