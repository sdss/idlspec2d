#! /bin/sh
#------------------------------------------------------------------------------
# Script to start the Spectro-1D processing.
#
# This script takes one argument, which contains the command-line arguments
# for the BATCH1D procedure.  For example:
#   sprobot1d.sh ",topdir='/home/data/2d_v4', nice=10"
#
# D. Schlegel, Princeton, 19 Dec 2000
#------------------------------------------------------------------------------

# Do not put this in the background, because we search for the "sprobot1d.sh"
# process to determine if this is already running!

echo SPROBOT1D: batch1d $1
echo batch1d $1 | idl

exit
#------------------------------------------------------------------------------
