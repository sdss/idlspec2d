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

n=`\ps -elf | grep sprobot2d.sh | grep -v grep | wc -l`
if [ X"$n" != X"" -a "$n" -gt 2 ]; then
   echo "SPROBOT2D: BATCH2D already running at "`date`
   exit
fi
echo "SPROBOT2D: Launched at "`date` UID=$UID PPID=$PPID
echo "IDLSPEC2D_DIR="$IDLSPEC2D_DIR
echo "IDLUTILS_DIR="$IDLUTILS_DIR

# Do not put this in the background, because we search for the "sprobot2d.sh"
# process to determine if this is already running!

echo SPROBOT2D: batch2d $1
echo batch2d $1 | idl

echo "SPROBOT2D: Finished at "`date` UID=$UID PPID=$PPID

exit
#------------------------------------------------------------------------------
