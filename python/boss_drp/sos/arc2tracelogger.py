#!/usr/bin/env python3
import os.path as ptt
import sys
import time

class Logger(object):
    def __init__(self, filename):
        self.filename = filename

    class Transcript:
        def __init__(self, filename, cmd=None):
            self.terminal = sys.stdout
            self.log = open(filename, "a")
            if cmd is not None:
                cmd[0] = ptt.basename(cmd[0])
                self.log.write(' '.join(cmd)+'\n')
                self.log.write('\n')
                self.log.write('Opening log '+time.ctime())
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

        def flush(self):
            # this flush method is needed for python 3 compatibility.
            # this handles the flush command by doing nothing.
            # you might want to specify some extra behavior here.
            pass

    def start(self, cmd=None):
        sys.stdout = self.Transcript(self.filename, cmd=cmd)

    def stop(self):
        sys.stdout.log.close()
        sys.stdout = sys.stdout.terminal

class Logger(object):
    def __init__(self, filename):
        self.filename = filename

    class Transcript:
        def __init__(self, filename, cmd=None):
            self.terminal = sys.stdout
            try:
                self.log = open(filename, "a")
            except:
                self.log = open(filename, "w")
            if cmd is not None:
                cmd[0] = ptt.basename(cmd[0])
                self.log.write(' '.join(cmd))
                self.log.write('Opening log '+time.ctime())
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

        def flush(self):
            # this flush method is needed for python 3 compatibility.
            # this handles the flush command by doing nothing.
            # you might want to specify some extra behavior here.
            pass

    def start(self, cmd=None):
        sys.stdout = self.Transcript(self.filename, cmd=cmd)

    def stop(self):
        sys.stdout.log.close()
        sys.stdout = sys.stdout.terminal

