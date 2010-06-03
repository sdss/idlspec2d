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

cd $BOSS_SPECTRO_REDUX
echo "platelist, /create" | idl
echo "platemerge" | idl 2> /dev/null
#echo "platemerge, /public" | idl 2> /dev/null
#echo "platemerge, public='EDR'" | idl 2> /dev/null
#echo "platemerge, public='DR1'" | idl 2> /dev/null
#echo "platemerge, public='DR2'" | idl 2> /dev/null
#echo "platemerge, public='DR3'" | idl 2> /dev/null
#echo "platemerge, public='DR4'" | idl 2> /dev/null
#echo "platemerge, public='DR5'" | idl 2> /dev/null
echo "zplot" | idl 2> /dev/null
echo "platehist" | idl 2> /dev/null

#------------------------------------------------------------------------------
