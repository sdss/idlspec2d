#! /bin/sh
#------------------------------------------------------------------------------
# Script to copy raw spectro data from the machine $SPROBOT_HOST (currently
# sos.apo.nmsu.edu), where the data is assumed to be in the directories
#    /data/spectro/astrolog/$MJD
#    /data/spectro/$MJD
# The astrolog data is copied into the local directory
#    $ASTROLOG_DIR/$MJD
# and the image files are copied to the first disk in $SPROBOT_LOCALDISKS
# with more than 2.5 Gb free, with a pointer to this from
#    $RAWDATA_DIR/$MJD -> $localdisks[i]/$MJD
#
# Once the data is copied over, list the new MJDs in a file for other programs
# to read.
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
copiedMJDs=$RAWDATA_DIR/copiedMJDs.list

#------------------------------------------------------------------------------
# Test that certain environment variables are already set.

if [ -z "$SSH_AGENT_PID" ] ; then
  echo "SSH_AGENT_PID must be set!"
  exit
fi

if [ -z "$SSH_AUTH_SOCK" ] ; then
  echo "SSH_AUTH_SOCK must be set!"
  exit
fi

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

# Sanity check
if ! `ssh-add -l | grep -q -n dssadmin@`; then
    echo "no dssadmin ssh-agent for spectro!"
    exit 1
fi

echo ""
echo "-------------------------------------------------------------------------------"
echo "SPROBOT: Launched at "`date` UID=$UID PPID=$PPID
echo "IDLSPEC2D_DIR="$IDLSPEC2D_DIR
echo "IDLUTILS_DIR="$IDLUTILS_DIR

#------------------------------------------------------------------------------
# Find raw data directories on the machine $hostname, and loop through them

mjdlist=''

# Remove the leading "/data/spectro" and the trailing "/" from the directory list.
# remotedir=`$SPROBOT_RSH $hostname ls -d /data/spectro/[5-9]???? | sed 's/\/data\/spectro\///g'  | sed 's/\///g'`

# Instead of looking for "/data/spectro/$MJD" to get an MJD list, look for
# "/data/spectro/astrolog/$MJD".  This will get the log files for MJDs
# where there are no spectroscopic data.
# The following line fails, because the directory list too long!
# So instead, "cd" into the directory first.
# remotedir=`$SPROBOT_RSH $hostname ls -d /data/spectro/astrolog/[5-9]???? | tail -7 | sed 's/\/data\/spectro\/astrolog\///g' | sed 's/\///g'`
remotedir=`$SPROBOT_RSH $hostname "(cd /data/spectro/astrolog ; ls -d [5-9]???? | tail -7 | sed 's/\///g')"`
echo REMOTEDIR=$remotedir

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
      echo SPROBOT: rsync "$hostname:/data/spectro/astrolog/$mjdstr" $astrologdir
      rsync -ar --rsh="$SPROBOT_RSH" \
       "$hostname:/data/spectro/astrolog/$mjdstr" $astrologdir

      # Normalize permissions
      find $astrologdir/$mjdstr -type f -print | xargs chmod 644
      find $astrologdir/$mjdstr -type d -print | xargs chmod 755

      # Move the important astrolog files from $ASTROLOG_DIR->$SPECLOG_DIR
      # to make a CVS-versioned copy.
      echo "SPROBOT: speclog_update $mjdstr"
      speclog_update $mjdstr

      # Copy the raw FITS files... copy only files ending in ".fit.gz"
      echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
      echo SPROBOT: rsync "$hostname:/data/spectro/$mjdstr/*.fit.gz" $localdir
      rsync -ar --rsh="$SPROBOT_RSH" \
       "$hostname:/data/spectro/$mjdstr/*.fit.gz" $localdir
#       "$hostname:/data/spectro/$mjdstr/*" $localdir

      # Copy the guider image directory...
      echo "SPROBOT: Current time "`date` UID=$UID PPID=$PPID
      echo SPROBOT: rsync "$hostname:/data/spectro/$mjdstr/guider" $localdir
      rsync -ar --rsh="$SPROBOT_RSH" \
       "$hostname:/data/spectro/$mjdstr/guider/*.fits.gz" $localdir/guider

      # Compress the raw FITS files w/gzip...
#      echo SPROBOT: gzip $localdir/*.fit $localdir/*/*.fit
#      gzip $localdir/*.fit $localdir/*/*.fit

      # Normalize permissions.
      find $localdir -type f -print | xargs chmod 644
      find $localdir -type d -print | xargs chmod 755

      if [ -z "$mjdlist" ] ; then
         mjdlist=$mjdstr
      else
         mjdlist=$mjdlist,$mjdstr
      fi
   else
      echo "SPROBOT: WARNING: All disks are full!!!"
   fi

done

#------------------------------------------------------------------------------
# If $topoutdir is not defined, then exit.

if [ -z "$topoutdir" ] ; then
   echo "SPROBOT: TOPOUTDIR must be set!"
   exit
fi

#------------------------------------------------------------------------------
# Hand off to the processing stage simply by listing the mjds which we have copied.

# The writes should be atomic, and we are the only writer.
if [ -n "$mjdlist" ] ; then
    echo $mjdlist | perl -ane 'print join("\n", split(","));' >> $copiedMJDs
fi

echo "SPROBOT: Finished at "`date` UID=$UID PPID=$PPID
echo SPROBOT: MJDLIST=$mjdlist

