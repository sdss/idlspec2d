#!/usr/bin/env python3
from boss_drp.utils.chpc2html import chpc2html
from boss_drp.utils.path_to_html import path_to_html

import logging
import collections
import sys
import os.path as ptt
from os import rename, remove, getenv
from email.message import EmailMessage
import inspect
import datetime
import smtplib
import warnings
import numpy as np
import re
import linecache
splog_name = None

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

class DailyFormatter(logging.Formatter):
    def __init__(self):
        super().__init__(fmt="%(levelno)d: %(msg)s", datefmt=None, style='%')
    def format(self, record):
        format_orig = self._style._fmt
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        elif record.levelno == logging.DEBUG:
            self._style._fmt = '%(funcName)s: %(message)s'
        else:
            self._style._fmt = "%(levelname)s: %(message)s"
        result = logging.Formatter.format(self, record)
        self._style._fmt = format_orig
        return result

def build_email(subject, emails, content, from_domain, attachment, link=False):
    msg = EmailMessage()
    attachment_note = ''
    if content is None:
        content = subject
    msg['Subject'] = subject
    msg['From'] = f"BOSS Pipeline <{getenv('USER')}@{from_domain}>"
    msg['BCC'] = ', '.join(emails)
    if attachment is not None:
        attachment = np.atleast_1d(attachment)
        msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
        if link:
            attachment_note = (
                "\n\nAttachments removed:\n" + "\n".join(f"- {chpc2html(path_to_html(att))}" for att in attachment if ptt.exists(att))
            )
            msg.set_content(content + attachment_note)
        else:
            msg.set_content(content)
            for fa in attachment:
                if ptt.exists(fa):
                    with open(fa, 'rb') as fp:
                        logdata = fp.read()
                        msg.add_attachment(logdata, maintype='text', subtype='plain', filename=ptt.basename(fa))
    else:
        msg.set_content(content)
    return msg

def send_email(subject, email_file, attachment, content=None,
                from_domain="chpc.utah.edu", allemail=False):

    try:
        emails = open(email_file).read().splitlines()
    except:
        emails = []
        splog.info(email_file+' does not exist')
        
    emails = ' '.join(emails).split()
    if not allemail:
        emails = [emails[0]]
        
    msg = build_email(subject, emails, content, from_domain, attachment, link=False)
    try:
        s = smtplib.SMTP('localhost')
        s.send_message(msg)
        s.set_debuglevel(1)
        s.quit()
    except smtplib.SMTPException as e:
        msg = build_email(subject, emails, content, from_domain, attachment, link=True)
        s = smtplib.SMTP('localhost')
        s.set_debuglevel(1)
        s.send_message(msg)
        s.quit()
    return(None)

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

    def send(self, subject, email_file, log, allemail=False):
        try:
            emails = open(email_file).read().splitlines()
        except:
            log.error(email_file+' does not exist')
            emails = []
            
        try:
            send_email(subject, email_file, None, content=self.contents(), from_domain="chpc.utah.edu", allemail=allemail)

        except:
            outputs = []
            for line in self.contents():
                if 'slurm.session.Client:' in line:
                    continue
                if 'slurm.session.Client: task #' in line:
                    continue
                outputs.append(line)
            self.contents = '\n'.join(outputs)
            try:
                send_email(subject, email_file, None, content=self.contents, from_domain="chpc.utah.edu", allemail=allemail)
            except:
                self.contents = 'ERROR Building Email Log'
                send_email(subject, email_file, None, content=self.contents, from_domain="chpc.utah.edu", allemail=allemail)

        return

def backup_log(logfile):
    if not ptt.exists(logfile): return
    else:
        num = 0
        backname = logfile
        while ptt.exists(backname):
            num = num + 1
            backname = logfile + '.' + str(num)
        rename(logfile, backname)
    return

class IncrementalFormatter(logging.Formatter):

    def __init__(self, *args, pad=0, **kwargs):
        super().__init__(*args, **kwargs)
        self.pad = pad
    def format(self, record):
        # Get the current call stack depth
        stack = []
        for s in inspect.stack():
            if 'bin' in s.filename:
                stack.append(s)
            elif 'boss_drp' in s.filename:
                stack.append(s)
        stack_depth = len(stack) - 2 if len(stack) >= 3 else 0
        # Add tabs based on the stack depth
        indent = '  ' * (stack_depth+self.pad)
        formatted_message = super().format(record)
        original_message = record.getMessage()  # Extract unformatted message
        padh = ' '*(len(formatted_message)- len(original_message))
        #print(repr(original_message))
        lines = formatted_message.splitlines()
        #for l in lines: print(repr(l))
        if len(lines) > 1:
            formatted_message = lines[0] + "\n" + "\n".join(f"{indent}{padh}{line.lstrip()}" for line in lines[1:])
        
        return f"{indent}{formatted_message}"


class Splog:
    def __init__(self, name=None, sos=False, no_exception = False,  lname = None, cfg=None):
        global splog_name

        self.sos = sos
        self.name = name
        self.sec_fhandlers = None
        self.elog = None
        self.elog_handler = None
        self.set(name=name, sos=sos, no_exception=no_exception, lname = lname, cfg=cfg)

        self._pending_external_handlers = []  # Store external handlers to apply later
        self._external_loggers = None
        self._is_open = False  # Flag to track if file handlers have been set up

        self._original_warn = warnings.warn
        warnings.warn = self.Warning

    def set(self,name=None, sos=False, no_exception=False, lname = None, cfg=None):
        global splog_name
        if name is None:
            if splog_name is not None:
                name = splog_name
            if inspect.stack()[1].function == '<module>':
                name = ptt.splitext(ptt.basename(inspect.stack()[1].filename))[0]
            else:
                name = inspect.stack()[1].function
        splog_name = name
        if not sos:
            self._log = logging.getLogger(name)
            self._log.setLevel(logging.DEBUG)
            self.no_exception = no_exception
            self._formatter = IncrementalFormatter('%(funcName)s: %(message)s')#logging.Formatter('%(funcName)s: %(message)s')
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(self._formatter)
            self._log.addHandler(ch)
            self.console = ch

        else:
            self.console = None
            self._log = logging.getLogger(name)
            self.no_exception = False
            rollover = datetime.time(hour=18)
            ext = '.log' if cfg.utah else ''
        
            print("Starting to log to " + lname+ext+" ")
            h  = logging.handlers.TimedRotatingFileHandler(lname+ext,
                                                          when='midnight',
                                                          interval=1, backupCount=5,
                                                          atTime=rollover)
            hc = logging.handlers.TimedRotatingFileHandler(lname + "-error"+ext,
                                                          when='midnight',
                                                          interval=1, backupCount=5,
                                                          atTime=rollover)

            f = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s")
            h.setFormatter(f)
            hc.setFormatter(f)
            h.setLevel(cfg.logLevel)
            hc.setLevel(logging.ERROR)
            self._log.setLevel(cfg.logLevel)
            self._log.addHandler(h)
            self._log.addHandler(hc)
            
        self.debug = self._log.debug
        self.info = self._log.info
        self.log  = self._log.info
        self.warning = self._log.warning
        self.error = self._log.error
        self.critical = self._log.critical

    def set_SOS(self,name,lname, cfg):
        self.close()
        self.__init__(name=name, lname=lname, cfg=cfg, sos=True)

    def Warning(self, message, category=None, stacklevel=1):
        """
        Custom wrapper for warnings.warn that logs warnings via splog.warning
        and includes the warning class name in the log message.
        """
        for filter in warnings.filters:
            if filter[0] == 'default':
                if isinstance(filter[1], type(re.compile(''))):
                    if filter[1].search(message):
                        break
                    else:
                        continue
                # If filter category is a regex, check if it matches the warning message
                if isinstance(filter[2], type(re.compile(''))):  # cat is a regular expression
                    if filter[2].search(message):  # If regex matches the message, ignore the warning
                        break

                # If the filter category is a class, check if it's a subclass of the warning's category
                elif isinstance(filter[2], type) and isinstance(category, type):
                    if issubclass(category, filter[2]):  # If the warning's category matches, ignore it
                        break  # Ignore the warning

            elif filter[0] == 'ignore':
                if isinstance(filter[1], type(re.compile(''))):
                    if filter[1].search(message):
                        return
                    else:
                        continue
                # If filter category is a regex, check if it matches the warning message
                if isinstance(filter[2], type(re.compile(''))):  # cat is a regular expression
                    if filter[2].search(message):  # If regex matches the message, ignore the warning
                        return

                # If the filter category is a class, check if it's a subclass of the warning's category
                elif isinstance(filter[2], type) and isinstance(category, type):
                    if issubclass(category, filter[2]):  # If the warning's category matches, ignore it
                        return  # Ignore the warning
            # If the filter is set to 'error', raise an exception
            elif filter[0] == 'error':
                # If the filter category is a regex, check if it matches the warning message
                if isinstance(filter[2], type(re.compile(''))):
                    if filter[2].search(message):
                        raise category(message)  # Raise the warning as an exception
                # If the filter category is a class, check if it's a subclass of the warning's category
                elif isinstance(filter[2], type) and isinstance(category, type):
                    if issubclass(category, filter[2]):
                        raise category(message)  # Raise the warning as an exception
                elif (filter[1] is None) and (isinstance(filter[2], type)) and category is None:
                    raise filter[2](message)
        # Capture the class name of the warning (category)
        warn_class = category.__name__ if category else "Warning"
#        if (warn_class == 'DeprecationWarning'):
#            return
            
        # Get the frame corresponding to the stacklevel
        
        stack = inspect.stack()
        stack.reverse()
        call = None
        line_number = None
        line_code = None
        for s in stack:
            tc = inspect.getmodule(s[0]).__name__ if inspect.getmodule(s[0]) else ""
            if 'boss_drp' in tc:
                call = tc
                line_number = s[2]  # Line number
                frame = s[0]
                try:
                    filename = frame.f_code.co_filename
                    line_code = linecache.getline(filename, line_number).strip()  # Get the exact line without opening the file
                except: line_code =''
            elif call is not None:
                break
        if call is None:
            call = '__main__'
            frame = inspect.stack()[-1]
            line_number = frame[2]
            line_code = ''
        frame = inspect.stack()[stacklevel]
        module = inspect.getmodule(frame[0]).__name__ if inspect.getmodule(frame[0]) else "<unknown module>"
        module = module.split('.')[0]
        module = f'[{module}]' if module.lower() != __name__.split('.')[0] else ''
        # Log the warning to splog with the class name and message
        if (warn_class == 'DeprecationWarning'):
            if self._log.getEffectiveLevel() < logging.getLevelName('WARNING'):
                print(f"{module}**{warn_class}**: {message} (Line {line_number}: {line_code})", file=sys.stderr)
        else:
            self.warning(f"{module}**{warn_class}**: {message} (Line {line_number}: {line_code})",stacklevel=stacklevel + 1)


    def exception(self, exctype=None, value=None, tb=None):
        if exctype is None or value is None or tb is None:
            import sys
            exctype, value, tb = sys.exc_info()
            if tb is None:
                self._log.error("No traceback available")
                return

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
        
    def emailer(self):
        self.elog = emailLogger()
        emaillog = self.elog.log_handler
        emaillog.setLevel(logging.DEBUG)
        emaillog.setFormatter(DailyFormatter())
        self._log.addHandler(emaillog)
        self.elog_handler = emaillog

    def send_email(self, subject, email_file, allemail=False):
        if self.elog is not None:
            self.elog.send(subject, email_file, self, allemail=allemail)
            self._log.removeHandler(self.elog_handler)
            self.elog_handler.close()
            self.elog = None
            self.elog_handler = None

    def close_elogger(self):
        if self.elog is not None:
            self._log.removeHandler(self.elog_handler)
            self.elog_handler.close()
            self.elog = None
            self.elog_handler = None

    def open(self, logfile=None, logprint=False, backup=False, append=False):
        if not self.sos:
            if logfile is not None:
                if backup:
                    backup_log(logfile)
                else:
                    if append is False:
                        if ptt.exists(logfile):
                            remove(logfile)
                fh = logging.FileHandler(logfile)
                fh.setLevel(logging.DEBUG)
                fh.setFormatter(self._formatter)
                self._log.addHandler(fh)

        if not self.no_exception:
            self._bkexecpthook = sys.excepthook
            sys.excepthook = self.exception

        if logprint is True:
            sys.stdout = StreamToLogger(self._log, logging.INFO)
            sys.stderr = StreamToLogger(self._log, logging.ERROR)

        # Now that the file handlers are open, apply any deferred external handler redirections
        if self._pending_external_handlers:
            for ex_logger in self._pending_external_handlers:
                self._redirect_stream_to_file_handler(ex_logger)
            self._pending_external_handlers = []  # Clear the list after applying
        if self._external_loggers:
            for ex_logger in self._external_loggers:
                if not ex_logger.hasHandlers():
                    self._forward_to_logger(ex_logger)
        
        self._is_open = True  # Mark the logger as open

    def add_file(self, filename, mode='a'):
        f = logging.Formatter('%(funcName)s-%(asctime)s: %(message)s')
        fh = logging.FileHandler(filename, mode=mode)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(f)
        self._log.addHandler(fh)
        self.sec_fhandlers = fh

    def close_file(self):
        if self.sec_fhandlers is not None:
            self._log.removeHandler(self.sec_fhandlers)

    def pause_file(self):
        if self.sec_fhandlers is not None:
            self.sec_fhandlers.setLevel(logging.CRITICAL + 1)

    def unpause_file(self):
        if self.sec_fhandlers is not None:
            self.sec_fhandlers.setLevel(logging.DEBUG)

    def close(self):
        if not self.sos:
            for handler in self._log.handlers:
                handler.close()
                self._log.removeFilter(handler)
            while self._log.hasHandlers():
                self._log.removeHandler(self._log.handlers[0])

        if (not self.no_exception) & (hasattr(self, '_bkexecpthook')):
            sys.excepthook = self._bkexecpthook

    def add_external_handlers(self, external_logger):
        """
        Add handlers from an external logger and apply a hanging indent format.
        Attach them to the already opened file handler in Splog if available.
        """
        ex_logger = logging.getLogger(external_logger)
        ex_logger.setLevel(self._log.level)
        if self._external_loggers is None:
            self._external_loggers = {ex_logger}
        else:
            self._external_loggers.add(ex_logger)
        # Format for indented logs
        if ex_logger.hasHandlers():
            exformatter = IncrementalFormatter(f': %(funcName)s: [%(module)s] %(message)s', pad = 1)
            #exformatter = IncrementalFormatter(f': %(funcName)s: [{ex_logger.name}] %(message)s', pad = 1)
            # First, we add the external StreamHandler(s) to the external logger
            for handler in ex_logger.handlers:
                if isinstance(handler, logging.StreamHandler):
                    # Apply the formatter for indented logs
                    handler.setFormatter(exformatter)#self._formatter)

                    # If splog is open, immediately attach the stream handler to file handlers
                    if self._is_open:
                        self._redirect_stream_to_file_handler(ex_logger)
                    else:
                        # Otherwise, defer this redirection until later
                        self._pending_external_handlers.append(ex_logger)

            # Ensure that the external logger's handlers now propagate their logs into Splog's handlers
            ex_logger.propagate = False
            ex_logger.setLevel(self._log.level)
        else:
            self._forward_to_logger(ex_logger)

    def _forward_to_logger(self, ex_logger):
        ex_logger.setLevel(self._log.level)
        forward_handler = ForwardToLogger(target_logger=self._log, name=ex_logger.name)
        ex_logger.addHandler(forward_handler)
        ex_logger.setLevel(self._log.level)
        ex_logger.propagate = False

    def _redirect_stream_to_file_handler(self, ex_logger):
        """Helper method to redirect external logger's StreamHandler to Splog's file handlers"""
        for file_handler in self._log.handlers:
            if isinstance(file_handler, logging.FileHandler):
                for handler in ex_logger.handlers:
                    if isinstance(handler, logging.StreamHandler):
                        # Create a custom MultiStreamHandler that writes to both file handler stream
                        # and the external logger's own stream
                        multi_stream = MultiStreamHandler(file_handler.stream, handler.stream)
                        handler.setStream(multi_stream)  # Redirect to the custom multi-stream

# General-purpose forwarding handler
class ForwardToLogger(logging.Handler):
    def __init__(self, target_logger, name=None):
        super().__init__()
        self.target_logger = target_logger
        self.name = name

    def emit(self, record):
        # Forward the record to all handlers attached to the target logger
        if self.name:
            record.msg = f"[{self.name}] {record.msg}"

        for handler in self.target_logger.handlers:
            if handler.level <= record.levelno:  # Check the handler's level
                handler.handle(record)
                
        
class MultiStreamHandler:
    """
    A custom stream handler that writes to multiple streams at once.
    """
    def __init__(self, *streams):
        self.streams = streams

    def write(self, message):
        for stream in self.streams:
            stream.write(message)

    def flush(self):
        for stream in self.streams:
            stream.flush()

splog = Splog()
