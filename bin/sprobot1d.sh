#! /bin/sh
#------------------------------------------------------------------------------
# Script to start the Spectro-1D processing.
#
# This script takes one argument, which contains the command-line arguments
# for the BATCH1D procedure.  For example:
#   sprobot1d.sh ",topdir='/u/dss/spectro',upsvers='v4_10_0', nice=19"
#
# D. Schlegel, Princeton, 19 Dec 2000
#------------------------------------------------------------------------------

# Exit if this process is already running.

n=`\ps -elf | grep sprobot1d.sh | grep -v grep | wc -l`
if [ X"$n" != X"" -a "$n" -gt 2 ]; then
   echo "SPROBOT1D: BATCH1D already running at "`date`
   exit
fi
echo "SPROBOT1D: Launched at "`date` UID=$UID PPID=$PPID
echo "IDLUTILS_DIR="$IDLUTILS_DIR

# Do not put this in the background, because we search for the "sprobot1d.sh"
# process to determine if this is already running!

echo SPROBOT1D: batch1d $1
echo batch1d $1 | idl

echo "SPROBOT1D: Finished at "`date` UID=$UID PPID=$PPID

exit
#------------------------------------------------------------------------------
