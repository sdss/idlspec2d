import logging
import io
import sys
import os.path as ptt
import inspect
from os import rename
import sys


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


    

def backup(logfile):
    if not ptt.exists(logfile): return
    else:
        num = 0
        backname = logfile
        while ptt.exists(backname):
            num = num + 1
            backname = logfile + '.' + str(num)
        rename(logfile, backname)
    return


class Splog:
    def __init__(self):

        if inspect.stack()[1].function == '<module>': name = ptt.splitext(ptt.basename(inspect.stack()[1].filename))[0]
        else: name = inspect.stack()[1].function

        self._log = logging.getLogger(name)
        self._log.setLevel(logging.DEBUG)
        
        
        self.debug = self._log.debug
        self.info = self._log.info
        self.log  = self._log.info
        self.warning = self._log.warning
        self.error = self._log.error
        self.debug = self._log.critical
        
        self._formatter = logging.Formatter('%(funcName)s: %(message)s')
        # create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(self._formatter)
        self._log.addHandler(ch)        
        

    def exception(self, exctype, value, tb):
        filename = tb.tb_frame.f_code.co_filename
        name     = tb.tb_frame.f_code.co_name
        line_no  = tb.tb_lineno
        error = ('--------------------------------------------\n'+
                '    Traceback'+'\n'+
                '    Type:      '+ str(exctype.__name__)+'\n'+
                '    Value:     '+ str(value)+'\n'+
                '    Traceback: '+f"File {filename} line {line_no}, in {name}")
        while True:
            try:
                tb = tb.tb_next
                filename = tb.tb_frame.f_code.co_filename
                name     = tb.tb_frame.f_code.co_name
                line_no  = tb.tb_lineno
                error = error+'\n               '+f"File {filename} line {line_no}, in {name}"
            except:
                break
        self._log.error(error)
        
        
    def open(self, logfile=None, logprint=False):
                
        if logfile is not None:
        
            backup(logfile)
        
            # create file handler 
            fh = logging.FileHandler(logfile)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(self._formatter)
            self._log.addHandler(fh)

            
        # Install exception handler
        self._bkexecpthook = sys.excepthook
        sys.excepthook = self.exception
#         # create console handler
#         ch = logging.StreamHandler()
#         ch.setLevel(logging.DEBUG)
#         ch.setFormatter(formatter)
#         self._log.addHandler(ch)
	
        if logprint is True:
            sys.stdout = StreamToLogger(self._log,logging.INFO)
            sys.stderr = StreamToLogger(self._log,logging.ERROR)

    def close(self):
        for handler in self._log.handlers:
            handler.close()
            self._log.removeFilter(handler)
        while self._log.hasHandlers():
            self._log.removeHandler(self._log.handlers[0])
            
            
        sys.excepthook = self._bkexecpthook
        #for handler in self._log.handlers:
        #    handler.close()
        #    self._log.removeFilter(handler)

# logger = logging.getLogger()
# while logger.hasHandlers():
#     logger.removeHandler(logger.handlers[0])
