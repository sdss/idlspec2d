#! /bin/sh
#------------------------------------------------------------------------------
# Check to see if "aporsync.sh" is already running.  If not, then start it.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

# ps has also changed on sos.apo.nmsu.edu

if \ps -A | grep aporsync.sh  | grep -v -e grep
then
  echo "APORSYNC: Already running at "`date`
  exit
fi

aporsync.sh

