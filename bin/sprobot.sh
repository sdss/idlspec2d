#! /bin/sh
#------------------------------------------------------------------------------
# Script to copy raw spectro data from "sdsshost.apo.nmsu.edu",
# where the data is assumed to be in the directories
#    /astrolog/$MJD
#    /data/spectro/$MJD
# The astrolog data is copied into the local directory
#    $SPECLOG_DIR/$MJD
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

astrologdir=$SPECLOG_DIR
toprawdir=$RAWDATA_DIR
localdisks='/peyton/scr/spectro0/data/rawdata /peyton/scr/spectro1/data/rawdata /peyton/scr/spectro2/data/rawdata'
topoutdir=$SPECTRO_DATA
hostname=sdsshost.apo.nmsu.edu

#------------------------------------------------------------------------------
# Test that certain environment variables are already set.

if [ -z "$SPECLOG_DIR" ] ; then
  echo "SPECLOG_DIR must be set!"
  exit
fi

if [ -z "$RAWDATA_DIR" ] ; then
  echo "RAWDATA_DIR must be set!"
  exit
fi

echo ""
echo "-------------------------------------------------------------------------------"
echo "SPROBOT: Started at "`date`

#------------------------------------------------------------------------------
# Find raw data directories on the machine $hostname, and loop through them

mjdlist=''

remotedir=`ssh $hostname ls -d /data/spectro/[56789]???? | sed -n 's/\/.*\///p'`
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
            if [ -z "$localdir" ] ; then
               localdir=$fdisk/$mjdstr
            fi
         fi
      done

      #----------
      # Create the local data directory for this night's data, and
      # create a symbolic link to that directory from the root data dir.

      if [ -n "$localdir" ] ; then
         echo SPROBOT: mkdir -p $localdir
         mkdir -p $localdir
         if [ $localdir != $toprawdir/$mjdstr ] ; then
            echo SPROBOT: ln -s $localdir $toprawdir/$mjdstr
            ln -fs $localdir $toprawdir/$mjdstr
         fi
      fi
   fi

   #----------
   # Proceed only if a good local directory exists to copy the data

   if [ -n "$localdir" ] ; then
      echo SPROBOT: rsync "$hostname:/astrolog/$mjdstr" $astrologdir
      rsync -ar --rsh="ssh -c arcfour" \
       --rsync-path=/p/rsync/v2_4_3/rsync \
       "$hostname:/astrolog/$mjdstr" $astrologdir
      echo SPROBOT: rsync "$hostname:/data/spectro/$mjdstr/*" $localdir
      rsync -ar --rsh="ssh -c arcfour" \
       --rsync-path=/p/rsync/v2_4_3/rsync \
       "$hostname:/data/spectro/$mjdstr/*" $localdir

      if [ -z "$mjdlist" ] ; then
         mjdlist=$mjdstr
      else
         mjdlist=$mjdlist,$mjdstr
      fi
   else
      echo "SPROBOT: All disks are full!!!"
   fi

done

#------------------------------------------------------------------------------
# If no new data or if $topoutdir is not defined, then exit.

if [ -z "$mjdlist" ] ; then
   exit
fi

if [ -z "$topoutdir" ] ; then
   exit
fi

#------------------------------------------------------------------------------
# Build the plan files

echo SPROBOT: MJDLIST=$mjdlist

echo ""
echo SPROBOT: "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]"
echo "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]" | idl
echo ""
echo SPROBOT: "spplan1d, topindir='$topoutdir', mjd=["$mjdlist"]"
echo "spplan1d, topindir='$topoutdir', mjd=["$mjdlist"]" | idl

#------------------------------------------------------------------------------
# Start the batch processing for Spectro-1D if it's not already running.
# Put this parent process in the background so this script will continue.

if ps ax | grep sprobot1d.sh  | grep -v -e grep
then
   echo "SPROBOT: BATCH1D already running at "`date`
else
   echo "SPROBOT: BATCH1D started at "`date`
   sprobot1d.sh ",topdir='$topoutdir', nice=10" &
fi

#------------------------------------------------------------------------------
# Start the batch processing for Spectro-2D if it's not already running.
# Put this parent process in the background so this script will continue.

if ps ax | grep sprobot2d.sh  | grep -v -e grep
then
   echo "SPROBOT: BATCH2D already running at "`date`
else
   echo "SPROBOT: BATCH2D started at "`date`
   sprobot2d.sh ",topdir='$topoutdir', nice=10" &
fi

exit
#------------------------------------------------------------------------------
