#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "cron.table".
# It syncs files from sdsshost to plate-mapper at apo.nmsu.edu.
#
# We set up two jobs, one for the blue files and one for the red.
# We do this so that we can utilize both processors on plate-mapper.
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC: Launched at "`date`

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
   

# This syncs /astrolog/[56789]???? from sdsshost to plate-mapper, excluding
# the id* files.  Only get /data/spectro/ directories
# This is hardwired to match the sdsshost directory structure

datadirs=`ssh sdsshost ls -d /data/spectro/[5-9]????`
astrologdirs=`echo $datadirs | sed -n 's/\/data\/spectro/\/astrolog/pg'`

for dir in $astrologdirs
do

  rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --exclude="*id*" \
      --exclude="*log" --log-format="/astrolog/%f" \
      sdsshost:$dir $speclog_dir

done


#     --log-format="/astrolog/%f" "sdsshost:/astrolog/[56789]????" /astrolog/   

# This syncs /data/spectro/[56789]???? from sdsshost to plate-mapper, excluding
# the red and guider files.
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --exclude="*guider*" \
      --log-format="/data/spectro/%f" --exclude="*-r*" \
      "sdsshost:/data/spectro/[56789]????" $rawdata_dir | startapo.sh &

# This sleep is to try to keep the blue side copying over before the red
# side.  We do this because APOREDUCE creates its summary files after the r2
# CCD is reduced, and assumes that all other CCD's have been reduced by then.
sleep 1

# This syncs /astrolog/[56789]???? from sdsshost to plate-mapper, excluding
# the blue and guider files.
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="/data/spectro/%f" --exclude="*-b*" \
      "sdsshost:/data/spectro/[56789]????" $rawdata_dir | startapo.sh 

