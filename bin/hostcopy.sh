#! /bin/sh
#------------------------------------------------------------------------------
# Script to copy raw spectro data from "sdsshost.apo.nmsu.edu",
# where the data is assumed to be in the directories
#    /astrolog/$MJD
#    /data/spectro/$MJD
# The astrolog data is copied into the local directory
#    $ASTROLOG_DIR/$MJD
# and the image files are copied to the first disk in $localdisks
# with more than 4 Gb free, with a pointer to this from
#    $RAWDATA_DIR/$MJD -> $localdisks[i]/$MJD
#
# Once the data is copied over, build plan files under the directories
#    $SPECTRO_DATA/$PLATE
# then launch "batch2d" and "batch1d" to reduce them.
#
# D. Schlegel, Princeton, 19 Dec 2000
#------------------------------------------------------------------------------
# Set file names.

astrologdir=$ASTROLOG_DIR
toprawdir=$RAWDATA_DIR
localdisks='/peyton/scr/spectro0/data/rawdata /peyton/scr/spectro1/data/rawdata /peyton/scr/spectro3/data/rawdata'
topoutdir=$SPECTRO_DATA

#------------------------------------------------------------------------------
# Test that certain environment variables are already set.

if [ -z "$ASTROLOG_DIR" ] ; then
  echo "ASTROLOG_DIR must be set!"
  exit
fi

if [ -z "$RAWDATA_DIR" ] ; then
  echo "RAWDATA_DIR must be set!"
  exit
fi

#------------------------------------------------------------------------------
# Find raw data directories on sdsshost, and loop through them

mjdlist=''

remotedir=`ssh sdsshost.apo.nmsu.edu ls -d /data/spectro/[56789]???? | sed -n 's/\/.*\///p'`
for mjdstr in $remotedir ; do

   #----------
   # If the local directory does not exist, then create it

   localdir=`ls -d $toprawdir/$mjdstr 2> /dev/null` | head -1
   if [ -z "$localdir" ] ; then

      #----------
      # Find the first local disk with more than 4Gb free
      for fdisk in $localdisks ; do
         qgood=`df -m $fdisk | awk '{if (FNR==2 && $4>4000) {print 1}}'`
         if [ -n "$qgood" ] ; then
            if [ -z $localdir ] ; then
               localdir=$fdisk/$mjdstr
            fi
         fi
      done

      #----------
      # Create the local data directory for this night's data, and
      # create a symbolic link to that directory from the root data dir.

      if [ -n "$localdir" ] ; then
echo HOSTCOPY mkdir -p $localdir
#          mkdir -p $toprawdir/$mjdstr
         if [ $localdir != $toprawdir/$mjdstr ] ; then
echo HOSTCOPY ln -s $localdir $toprawdir/$mjdstr
#            ln -fs $localdir $toprawdir/$mjdstr
         fi
      fi
   fi

   #----------
   # Proceed only if a good local directory exists to copy the data

   if [ -n "$localdir" ] ; then
       echo HOSTCOPY /astrolog/$mjdstr $astrologdir
       echo HOSTCOPY /data/spectro/$mjdstr $localdir
#      rsync -ar --rsh="ssh -c arcfour" \
#       --rsync-path=/p/rsync/v2_4_3/rsync \
#       sdsshost.apo.nmsu.edu:/astrolog/$mjdstr $astrologdir
#      rsync -ar --rsh="ssh -c arcfour" \
#       --rsync-path=/p/rsync/v2_4_3/rsync \
#       sdsshost.apo.nmsu.edu:/data/spectro/$mjdstr $localdir

      if [ -z $mjdlist ] ; then
         mjdlist=$mjdstr
      else
         mjdlist=$mjdlist,$mjdstr
      fi
   fi

done

#------------------------------------------------------------------------------
# If no new data or if $topoutdir is not defined, then exit.

if [ -z $mjdlist ] ; then
   exit
fi

if [ -z $topoutdir ] ; then
   exit
fi

#------------------------------------------------------------------------------
# Build the plan files

echo MJDLIST $mjdlist

echo "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]"
echo "spplan1d, topoutdir='$topoutdir', mjd=["$mjdlist"]"
#echo "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]" | idl
#echo "spplan1d, topoutdir='$topoutdir', mjd=["$mjdlist"]" | idl

#------------------------------------------------------------------------------
# Start the batch processing for Spectro-2D and Spectro-1D, but only
# if not already running!!!???

echo "batch2d, topdir='$topoutdir', nice=19"
echo "batch1d, topdir='$topoutdir', nice=19"
#echo "batch2d, topdir='$topoutdir', nice=19" | idl
#echo "batch1d, topdir='$topoutdir', nice=19" | idl

exit
#------------------------------------------------------------------------------
