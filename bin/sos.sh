#! /bin/sh
#------------------------------------------------------------------------------
# Script to launch the three cron jobs: one for guider and logs files,
# one for blue files, and one for red files.  Start each job only if it
# isn't already running.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

if \ps -A | grep aporsync_logs.sh  | grep -v -e grep
then
  echo "APORSYNC_LOGS: Already running at "`date`
else
  aporsync_logs.sh &
fi

if \ps -A | grep aporsync_blue.sh  | grep -v -e grep
then
  echo "APORSYNC_BLUE: Already running at "`date`
else
  aporsync_blue.sh &
fi

if \ps -A | grep aporsync_red.sh  | grep -v -e grep
then
  echo "APORSYNC_RED: Already running at "`date`
else
  aporsync_red.sh &
fi

