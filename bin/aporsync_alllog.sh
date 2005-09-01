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

echo "APORSYNC_ALLLOG: Launched at "`date -u` UID=$UID PPID=$PPID

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
# Copy all log files in that directory, including files not
# copied by "aporsync_logs.sh".
# We follow symbolic links in this copy (--copy-unsafe-links),
# since it now appears that some brain surgeon decided that some
# subset of these files should live elsewhere, with only symbolic
# links in the /astrolog/$MJD directory.

# datadirs=`ssh sdsshost.apo.nmsu.edu ls -d /data/spectro/[5-9]????`
# astrologdirs=`echo $datadirs | sed -n 's/\/data\/spectro/\/astrolog/pg'`

# Instead, rsync the largest 7 MJDs in the /astrolog directory, whether
# or not there is a corresponding data directory in /data/spectro for the same MJD.
astrologdirs=`ssh sdsshost.apo.nmsu.edu 'ls -d /astrolog/[5-9]???? | sort | tail -7'`

for thisdir in $astrologdirs
do
   echo Copy all log files from $thisdir
   rsync -ar --rsh="ssh -c blowfish" \
    --copy-unsafe-links \
    sdsshost.apo.nmsu.edu:$thisdir $ASTROLOG_DIR
done

echo "APORSYNC_ALLLOG: Finished at "`date -u` UID=$UID PPID=$PPID

