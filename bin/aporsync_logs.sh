#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "sos_start".
# It syncs files from sdsshost to the local machine (i.e., sos.apo.nmsu.edu).
#
# Only sync the /astrolog/$MJD directories when there is a corresponding
# /data/spectro/$MJD directory.  This is because the data directories are
# periodically purged on sdsshost, whereas the astrolog directories are not.
# For reducing data with SOS, we obviously only need the astrolog directories
# that correspond to data that presently exists on disk.
#
# Note that the directory locations on sdsshost.apo are hardwired
# as /astrolog and /data/spectro.
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC_LOGS: Launched at "`date -u` UID=$UID PPID=$PPID

if [ -z "$RAWDATA_DIR" ] ; then
   echo "Abort: RAWDATA_DIR not set!"
   exit
fi

if [ -z "$ASTROLOG_DIR" ] ; then
   echo "Abort: ASTROLOG_DIR not set!"
   exit
fi

# This syncs /astrolog/[5-9]???? from sdsshost to the local machine.
# Only consider those MJD subdirectories /astrolog/$MJD where a corresponding
# data directory exists in /data/spectro/$MJD.
# Only copy the following select set of files:
#   Unplugged*.ps
#   fiberScan*.par
#   guiderMon*.par
#   op*.par
#   plPlugMap*.par
#   sdReport*.par

datadirs=`ssh sdsshost.apo.nmsu.edu ls -d /data/spectro/[5-9]????`
astrologdirs=`echo $datadirs | sed -n 's/\/data\/spectro/\/astrolog/pg'`

for dir in $astrologdirs
do
   rsync -ar --rsh="ssh -c blowfish" \
    --rsync-path=/p/rsync/v2_4_3/rsync \
    --include "Unplugged*.ps"  --include "fiberScan*.par" \
    --include "guiderMon*.par" --include "op*.par" \
    --include "plPlugMap*.par" --include "sdReport*.par" \
    --include "exposureLog*.par" \
    --exclude="*/*" \
    --include "*/" \
    --log-format="/astrolog/%f" \
    sdsshost.apo.nmsu.edu:$dir $ASTROLOG_DIR
done

# This syncs /astrolog/[5-9]???? from sdsshost to the local machine,
# exluding the blue and red files (only include guider files).
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="/data/spectro/%f" \
      --include "*/guider" \
      --include "gimg*" \
      --exclude="*/*" \
      "sdsshost.apo.nmsu.edu:/data/spectro/[5-9]????" $RAWDATA_DIR | startapo.sh

echo "APORSYNC_LOGS: Finished at "`date -u` UID=$UID PPID=$PPID

