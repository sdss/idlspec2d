#! /bin/sh
#------------------------------------------------------------------------------
# Check to see if "aporsync.sh" is already running.  If not, then start it.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

if ps ax | grep aporsync.sh  | grep -v -e grep
then
  echo "APORSYNC: Already running at "`date`
  exit
fi

aporsync.sh

