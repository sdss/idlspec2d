#! /bin/sh
#------------------------------------------------------------------------------
# Script to start the Spectro-2D processing.
#
# This script takes one argument, which contains the command-line arguments
# for the BATCH2D procedure.  For example:
#   sprobot2d.sh ",topdir='/home/data/2d_v4', nice=10"
#
# D. Schlegel, Princeton, 19 Dec 2000
#------------------------------------------------------------------------------

# Exit if this process is already running.

n=`ps -elf | grep sprobot2d.sh | grep -v grep | wc -l`
if [ X"$n" != X"" -a "$n" -gt 2 ]; then
  echo "Already running sprobot2d.sh.  Exiting."
  exit
fi

# Do not put this in the background, because we search for the "sprobot2d.sh"
# process to determine if this is already running!

echo SPROBOT2D: batch2d $1
echo batch2d $1 | idl

exit
#------------------------------------------------------------------------------
