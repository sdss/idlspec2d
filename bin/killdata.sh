#! /bin/sh
#------------------------------------------------------------------------------
# This is a deadly script to blow away all local copies of raw data in
# $RAWDATA_DIR/$MJD which no longer exist on sdsshost.apo.
# It is invoked from mailhtml, and called every morning.
# Also check spectrologs, and remove 1 month old files
#
# S. Burles, APO, 4 May 2000
#                25 Jan 2002, Retain last 25 MJDs of spectrologs
#------------------------------------------------------------------------------

if [ -z "$RAWDATA_DIR" ] ; then
   echo "Abort: RAWDATA_DIR not set!"
   exit
fi

if [ -z "$SPECTROLOG_DIR" ] ; then
   echo "Abort: SPECTROLOG_DIR not set!"
   exit
fi

if [ -z "$ASTROLOG_DIR" ] ; then
   echo "Abort: ASTROLOG_DIR not set!"
   exit
fi

datadirs=`ls -d $RAWDATA_DIR/[56789]????`

#------------------------------------------------------------------------------
# Find all data directories which exist locally but do not exist on
# sdsshost.apo.  Those directories will also be deleted locally.

for deaddir in $datadirs
do

    if  [ `ssh sdsshost ls -d $deaddir 2>/dev/null` ]   
    then
      echo KILLDATA: $deaddir still exists
    else
      echo KILLDATA: I should kill $deaddir
      rm -rf $deaddir
    fi

done

spectrologs=`ls -d $SPECTROLOG_DIR/[56789]???? | sort -r | tail +25`

for dir in $spectrologs
do

    introuble=`ls $dir/fflat*.fits $dir/sci*.fits $dir/tset*.fits $dir/wset*.fits $dir/*.ps $dir/logfile*.html 2>1`

    if  [ $introuble ]
    then
      echo KILLDATA: cleaning up $dir
      rm -f $introuble
    fi
done

#------------------------------------------------------------------------------
# Find all astrolog directories which exist locally but do not exist on
# sdsshost.apo.  Those directories will also be deleted locally.
#
#datadirs=`ls -d $ASTROLOG_DIR/[56789]????`
#
#for deaddir in $datadirs
#do
#    if  [ `ssh sdsshost ls -d $deaddir 2>/dev/null` ]   
#    then
#      echo KILLDATA: $deaddir still exists
#    else
#      echo KILLDATA: I should kill $deaddir
#      rm -rf $deaddir
#    fi
#done

