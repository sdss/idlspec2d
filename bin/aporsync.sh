#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "cron.table".
# It syncs files from sdsshost to the local machine (i.e., sos.apo.nmsu.edu).
#
# We set up two jobs, one for the blue files and one for the red.
# We assume that the local machine has two processors, and we utilize both.
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC: Launched at "`date` UID=$UID PPID=$PPID

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

# This syncs /data/spectro/[5-9]???? from sdsshost to the local machine,
# exluding the red and guider files.
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --exclude="*guider*" \
      --log-format="/data/spectro/%f" --exclude="*-r*" \
      "sdsshost:/data/spectro/[5-9]????" $rawdata_dir | startapo.sh &

# Historically, we have had a sleep statement here to keep the blue side
# copying over before the red side.  It used to be that APOREDUCE created
# its summary files after the r2 CCD is reduced, and assumed that all other
# CCD's have been reduced by then.  This is no longer true -- we create the
# summary file after each file is reduced.
sleep 1

# This syncs /astrolog/[5-9]???? from sdsshost to the local machine,
# exluding the blue and guider files.
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="/data/spectro/%f" --exclude="*-b*" \
      "sdsshost:/data/spectro/[5-9]????" $rawdata_dir | startapo.sh 

echo "APORSYNC: Finished at "`date` UID=$UID PPID=$PPID

