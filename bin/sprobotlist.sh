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
echo "IDLSPEC2D_DIR="$IDLSPEC2D_DIR
echo "IDLUTILS_DIR="$IDLUTILS_DIR

cd $SPECTRO_DATA
echo "platelist, /create" | idl
echo "platemerge" | idl 2> /dev/null
echo "platemerge, /public" | idl 2> /dev/null
echo "zplot" | idl 2> /dev/null
echo "platehist" | idl 2> /dev/null

#------------------------------------------------------------------------------
