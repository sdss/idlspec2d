#! /bin/bash
#------------------------------------------------------------------------------
# Script to create plate list and summary files once per day.
#
# D. Schlegel, Princeton, 1 Mar 2002
#------------------------------------------------------------------------------
# Generate summary files.

echo ""
echo "-------------------------------------------------------------------------------"
echo "SPROBOTLIST: Started at "`date`

cd $SPECTRO_DATA
echo "platelist, /create" | idl
echo "platemerge" | idl 2> /dev/null
echo "platemerge, /public" | idl 2> /dev/null

#------------------------------------------------------------------------------
