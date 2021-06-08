#!/usr/bin/env python

import os, sys, commands, time

""" 
killChildren:

Something of a quick hack to kill the children of a process.  One problem is that
if processes are creating new processes as we are kill them, not everything will be
cleaned up.

Written by Gary Kushner (LBL).  Sept 2010.

"""



####
def usage():
	"""Display usage and exit"""
	
	usageCMD = os.path.basename(sys.argv[0])

	print "usage:  "+ usageCMD + " pid"

	sys.exit(1)



def kill(pid):
	#	A bit of recursion will happen here
	killChildren(pid)
	
	#	Try SIGINT
	print "Sending SIGINT to process %d ..." % pid
	sys.stdout.flush()
	try:
		os.kill(int(pid), 2)
	except:
		pass
			
	#	If still running, sleep and check again
	rc = commands.getstatusoutput("ps ax -o pid | grep %d" % pid)
	if len(rc[1]) != 0:
		time.sleep(5)
		
	#	SIGKILL if still alive
	rc = commands.getstatusoutput("ps ax -o pid | grep %d" % pid)
	if len(rc[1]) != 0:
		print "Sending SIGKILL to process %d..." % pid
		sys.stdout.flush()
		try:
			os.kill(pid, 9)
		except:
			pass
		
	#	If still running, sleep and check again
	rc = commands.getstatusoutput("ps ax -o pid | grep %d" % pid)
	if len(rc[1]) != 0:
		time.sleep(5)
	rc = commands.getstatusoutput("ps ax -o pid | grep %d" % pid)
	if len(rc[1]) != 0:
		print "** Could not kill process %d **" % pid
		

def killChildren(pid):
	#	Get children
	rc = commands.getstatusoutput("ps ax -o pid,ppid | grep %d | awk \"{ if ( \$2 == %d ) { print \$1 }}\"" % (pid,pid))
	if len(rc[1]) == 0:
		print "No children of %d" % pid
		return
	
	pids = rc[1].split("\n")
	print "Killing children of %d: %s" % (pid, str(pids))
	for p in pids:
		kill(int(p))


####
def main(argv):
	if len(argv) != 1:
		usage()
	pid = int(argv[0])

	kill(pid)
	
### Start of script

if __name__=='__main__':
	main(sys.argv[1:])
	print os.path.basename(sys.argv[0]) + " finished running!"


