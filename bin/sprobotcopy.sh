#! /bin/bash
#------------------------------------------------------------------------------
# Script to copy reduced data to a machine at Fermi.
#
# The files copied are:
#   $SPECTRO_DATA/spAll*
#   $SPECTRO_DATA/*/spZ*.fits
#   $SPECTRO_DATA/*/spPlate*.fits
#   $SPECTRO_DATA/*/spSN2d*.ps
#
# D. Schlegel, Princeton, 2 Jan 2001
#------------------------------------------------------------------------------
# Set file names.

topoutdir=$SPECTRO_DATA
#destdir=fsgi03.fnal.gov:mydata/2d_v4
#destdir=sdssdp6.fnal.gov:/data/dp3.p/data/schlegel/2d_v4
destdir=sdssdp7.fnal.gov:/data/dp3.p/data/schlegel/2d_v4_7
htmldir=spectro.princeton.edu:/peyton/home/spectro/httpd/html

#------------------------------------------------------------------------------
# Generate summary files.

echo ""
echo "-------------------------------------------------------------------------------"
echo "SPROBOTCOPY: Started at "`date`

cd $topoutdir
echo "platelist, /create" | idl 2> /dev/null
# echo "platemerge" | idl 2> /dev/null

#------------------------------------------------------------------------------
# Copy plate summary list to HTML directory.

scp platelist.txt $htmldir

#------------------------------------------------------------------------------
# Copy to Fermi.

rsync -arv --rsh="ssh" \
 --include "platelist*" --exclude "*"  \
 $topoutdir/* $destdir
rsync -arv --rsh="ssh" \
 --include "spAll*" --exclude "*"  \
 $topoutdir/* $destdir
rsync -arv --rsh="ssh" --include "*/" \
 --include "*spZ*.fits" --exclude "*"  \
 $topoutdir/* $destdir
rsync -arv --rsh="ssh" --include "*/" \
 --include "*spPlate*.fits" --exclude "*"  \
 $topoutdir/* $destdir
rsync -arv --rsh="ssh" --include "*/" \
 --include "*spSN2d*.ps" --exclude "*"  \
 $topoutdir/* $destdir

#------------------------------------------------------------------------------
