#! /bin/sh

#------------------------------------------------------------------------------
# This routine is called by the cron daemon, as set up with "sos_start".
# It syncs files from sdsshost to the local machine (i.e., sos.apo.nmsu.edu).
#
# Note that the directory locations on sdsshost.apo is hardwired
# as /data/spectro.
#
# We need the executable code "rsync" in the default path
# (e.g., as /usr/bin/rsync).
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

echo "APORSYNC_RED: Launched at "`date` UID=$UID PPID=$PPID

if [ -z "$RAWDATA_DIR" ] ; then
   echo "Abort: RAWDATA_DIR not set!"
   exit
fi

# This syncs /astrolog/[5-9]???? from sdsshost to the local machine,
# exluding the blue and guider files.
# Many files might be passed to startapo.sh at once.
rsync -ar --rsh="ssh -c blowfish" \
      --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="/data/spectro/%f" \
      --exclude="*-b*" \
      --exclude="*guider*" \
      "sdsshost.apo.nmsu.edu:/data/spectro/[5-9]????" $RAWDATA_DIR | startapo.sh 

echo "APORSYNC_RED: Finished at "`date` UID=$UID PPID=$PPID

