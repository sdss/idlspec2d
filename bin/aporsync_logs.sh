#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "sos_start".
# It syncs files from $RAWDATA_HOST to the local machine (i.e., sos.apo.nmsu.edu).
#
# Only sync the /astrolog/$MJD directories when there is a corresponding
# /data/spectro/$MJD directory.  This is because the data directories are
# periodically purged on the host machine, whereas the astrolog directories are not.
# For reducing data with SOS, we obviously only need the astrolog directories
# that correspond to data that presently exists on disk.
#
# Note that the directory locations on the host machine are hardwired
# as /astrolog and /data/spectro.
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC_LOGS: Launched at "`date -u` UID=$UID PPID=$PPID

if [ -z "$RAWDATA_HOST" ] ; then
   echo "Abort: RAWDATA_HOST not set!"
   exit
fi

if [ -z "$RAWDATA_DIR" ] ; then
   echo "Abort: RAWDATA_DIR not set!"
   exit
fi

if [ -z "$ASTROLOG_DIR" ] ; then
   echo "Abort: ASTROLOG_DIR not set!"
   exit
fi

# This syncs /astrolog/[5-9]???? from the host machine to the local machine.
# Only consider those MJD subdirectories /astrolog/$MJD where a corresponding
# data directory exists in /data/spectro/$MJD.
# Only copy the following select set of files:
#   sdReport*.par
#   plPlugMap*.par
#   plPlugMap*.log
#   plSlitpos*.par
#   guiderMon*.par
#   exposureLog*.par
#   op*.par
#   Unplugged*.ps
#   fiberScan*.par
# Copy the following file from the local machine back to the host machine
# (since this file is created on the Son-of-Spectro machine), but only
# if that file actually exists:
#   sdHdrFix-$MJD.par

datadirs=`ssh ${RAWDATA_HOST} ls -d /data/spectro/[5-9]????`
astrologdirs=`echo $datadirs | sed -n 's/\/data\/spectro/\/astrolog/pg'`

for thisdir in $astrologdirs
do
   rsync -ar --rsh="ssh -c blowfish" \
    --include "sdReport*.par"    \
    --include "plPlugMap*.par"   \
    --include "plPlugMap*.log"   \
    --include "plSlitpos*.par"   \
    --include "guiderMon*.par"   \
    --include "exposureLog*.par" \
    --include "op*.par"          \
    --include "Unplugged*.ps"    \
    --include "fiberScan*.par"   \
    --exclude="*/*"              \
    --log-format="/astrolog/%f"  \
    ${RAWDATA_HOST}:$thisdir $ASTROLOG_DIR
   thismjd=`echo $thisdir | sed -n 's/\/.*\///p'`
   if [ -e $ASTROLOG_DIR/$thismjd/sdHdrFix-$thismjd.par ] ; then
      rsync -ar --rsh="ssh -c blowfish" \
       $ASTROLOG_DIR/$thismjd/sdHdrFix-$thismjd.par \
       ${RAWDATA_HOST}:$thisdir
   fi
done

# This syncs /astrolog/[5-9]???? from the host machine to the local machine,
# exluding the blue and red files (only include guider files).
rsync -ar --rsh="ssh -c blowfish" \
      --log-format="/data/spectro/%f" \
      --include "*/guider" \
      --include "gimg*" \
      --exclude="*/*" \
      "${RAWDATA_HOST}:/data/spectro/[5-9]????" $RAWDATA_DIR | startapo.sh

echo "APORSYNC_LOGS: Finished at "`date -u` UID=$UID PPID=$PPID

