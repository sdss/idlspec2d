#! /bin/sh
#------------------------------------------------------------------------------
# Script to copy raw spectro data from the machine $SPROBOT_HOST (currently
# sos.apo.nmsu.edu), where the data is assumed to be in the directories
#    /astrolog/$MJD
#    /data/spectro/$MJD
# The astrolog data is copied into the local directory
#    $ASTROLOG_DIR/$MJD
# and the image files are copied to the first disk in $SPROBOT_LOCALDISKS
# with more than 2.5 Gb free, with a pointer to this from
#    $RAWDATA_DIR/$MJD -> $localdisks[i]/$MJD
#
# Once the data is copied over, build plan files under the directories
#    $SPECTRO_DATA/$PLATE
# then launch "batch2d" and "batch1d" to reduce them.  Run the same version
# of the idlspec2d product as is set up for this script.  For example, even
# if v4_9 is declared current, run v4_8 if that is what's set up when this
# script is run.  So, you should run "sprobot_start" from the version of
# the code that you want to reduce the data.
#
# D. Schlegel, Princeton, 19 Dec 2000
#------------------------------------------------------------------------------
# Set file names.

astrologdir=$ASTROLOG_DIR
toprawdir=$RAWDATA_DIR
localdisks=$SPROBOT_LOCALDISKS
#localdisks='/peyton/scr/spectro3/data/rawdata /peyton/scr/spectro1/data/rawdata /peyton/scr/spectro0/data/rawdata /peyton/scr/spectro2/data/rawdata'
topoutdir=$SPECTRO_DATA
hostname=$SPROBOT_HOST
#hostname=sdsshost.apo.nmsu.edu
#hostname=sos.apo.nmsu.edu

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

if [ -z "$SPROBOT_HOST" ] ; then
  echo "SPROBOT_HOST must be set!"
  exit
fi

if [ -z "$SPROBOT_LOCALDISKS" ] ; then
  echo "SPROBOT_LOCALDISKS must be set!"
  exit
fi

# Default to using "ssh" protocol
if [ -z "$SPROBOT_RSH" ] ; then
  SPROBOT_RSH=ssh
fi

echo ""
echo "-------------------------------------------------------------------------------"
echo "SPROBOT: Launched at "`date` UID=$UID PPID=$PPID

#------------------------------------------------------------------------------
# Find raw data directories on the machine $hostname, and loop through them

mjdlist=''

remotedir=`$SPROBOT_RSH $hostname ls -d /data/spectro/[56789]???? | sed -n 's/\/.*\///p'`
for mjdstr in $remotedir ; do

   #----------
   # If the local directory does not exist, then create it

   localdir=`ls -d $toprawdir/$mjdstr 2> /dev/null | head -1`
   if [ -z "$localdir" ] ; then

      #----------
      # Find the first local disk with more than 2.5 Gb free
      for fdisk in $localdisks ; do
         qgood=`df -m $fdisk | awk '{if (FNR==2 && $4>2500) {print 1}}'`
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
         echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
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
      # Copy the astrolog files...
      echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
      echo SPROBOT: rsync "$hostname:/astrolog/$mjdstr" $astrologdir
      rsync -ar --rsh="$SPROBOT_RSH" \
       "$hostname:/astrolog/$mjdstr" $astrologdir
#       --rsync-path=/p/rsync/v2_4_3/rsync

      # Copy the raw FITS files... copy only files ending in ".fit.gz"
      echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
      echo SPROBOT: rsync "$hostname:/data/spectro/$mjdstr/*.fit.gz" $localdir
      rsync -ar --rsh="$SPROBOT_RSH" \
       "$hostname:/data/spectro/$mjdstr/*.fit.gz" $localdir
#       --rsync-path=/p/rsync/v2_4_3/rsync
#       "$hostname:/data/spectro/$mjdstr/*" $localdir

      # Copy the guider image directory...
      echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
      echo SPROBOT: rsync "$hostname:/data/spectro/$mjdstr/guider" $localdir
      rsync -ar --rsh="$SPROBOT_RSH" \
       "$hostname:/data/spectro/$mjdstr/guider/*.fits.gz" $localdir/guider

      # Compress the raw FITS files w/gzip...
#      echo SPROBOT: gzip $localdir/*.fit $localdir/*/*.fit
#      gzip $localdir/*.fit $localdir/*/*.fit

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
# If $topoutdir is not defined, then exit.

if [ -z "$topoutdir" ] ; then
   echo "SPROBOT: TOPOUTDIR must be set!"
   exit
fi

# The following would exit if no new data
if [ -z "$mjdlist" ] ; then
   echo "SPROBOT: No new data (MJDLIST)!"
   exit
fi

#------------------------------------------------------------------------------
# Build the plan files if there is new data

echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
echo SPROBOT: MJDLIST=$mjdlist

if [ -n "$mjdlist" ] ; then
   echo ""
   echo SPROBOT: "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]"
   echo "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]" | idl
   echo ""
   echo SPROBOT: "spplan1d, topindir='$topoutdir', mjd=["$mjdlist"]"
   echo "spplan1d, topindir='$topoutdir', mjd=["$mjdlist"]" | idl
fi

#------------------------------------------------------------------------------
# Decide what the current version of idlspec2d is, and run that version
# if it is a UPS-declared version.

vers=`echo "print,idlspec2d_version()" | idl 2> /dev/null`
if [ ${vers:0:5} != NOCVS ] ; then upsversion=$vers ; fi

#------------------------------------------------------------------------------
# Batch process 2D first, wait for it to complete, then batch process 1D
# in the background.  The calls to the 2d and 1d scripts will exit if
# those scripts are already running.

   echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
   sprobot2d.sh ",topdir='$topoutdir', upsversion='$upsversion', nice=19"
   cd $SPECTRO_DATA
   echo "platelist, /create" | idl
   sprobot1d.sh ",topdir='$topoutdir', upsversion='$upsversion', nice=19" &

#------------------------------------------------------------------------------
# Start the batch processing for Spectro-1D if it's not already running.
# Put this parent process in the background so this script will continue.

#if \ps -elf | grep sprobot1d.sh  | grep -v -e grep
#then
#   echo "SPROBOT: BATCH1D already running at "`date`
#else
#   echo "SPROBOT: BATCH1D started at "`date`
#   sprobot1d.sh ",topdir='$topoutdir', nice=19" &
#fi

#------------------------------------------------------------------------------
# Start the batch processing for Spectro-2D if it's not already running.
# Put this parent process in the background so this script will continue.

#if \ps -elf | grep sprobot2d.sh  | grep -v -e grep
#then
#   echo "SPROBOT: BATCH2D already running at "`date`
#else
#   echo "SPROBOT: BATCH2D started at "`date`
#   sprobot2d.sh ",topdir='$topoutdir', nice=19" &
#fi

#------------------------------------------------------------------------------

echo "SPROBOT: Finished at "`date` UID=$UID PPID=$PPID

exit
#------------------------------------------------------------------------------
