#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "sos_start".
# It syncs files from sdsshost to the local machine (i.e., sos.apo.nmsu.edu).
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC_LOGS: Launched at "`date` UID=$UID PPID=$PPID

if [ -n "$RAWDATA_DIR" ]  
then 
   rawdata_dir=$RAWDATA_DIR
else
   rawdata_dir='/data/spectro'
fi

if [ -n "$ASTROLOG_DIR" ]  
then 
   speclog_dir=$ASTROLOG_DIR
else
   speclog_dir='/data/spectro/astrolog'
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

datadirs=`ssh sdsshost ls -d /data/spectro/[5-9]????`
astrologdirs=`echo $datadirs | sed -n 's/\/data\/spectro/\/astrolog/pg'`

for dir in $astrologdirs
do
   rsync -ar --rsh="ssh -c blowfish" \
    --rsync-path=/p/rsync/v2_4_3/rsync \
    --include "Unplugged*.ps" --include "fiberScan*.par" \
    --include "guiderMon*.par" --include "op*.par" \
    --include "plPlugMap*.par" --include "sdReport*.par" \
    --exclude="*" --log-format="/astrolog/%f" \
    sdsshost:$dir $speclog_dir
done

# This syncs /astrolog/[5-9]???? from sdsshost to the local machine,
# exluding the blue and red files (only include guider files).
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="/data/spectro/%f" \
      --exclude="*-b*" \
      --exclude="*-r*" \
      "sdsshost:/data/spectro/[5-9]????" $rawdata_dir

echo "APORSYNC_LOGS: Finished at "`date` UID=$UID PPID=$PPID

