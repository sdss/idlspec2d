#! /bin/sh

if [ -n "$RAWDATA_DIR" ]
then
   rawdata_dir=$RAWDATA_DIR
else
   rawdata_dir='/data/spectro'
fi

if [ -n "$RAWCOPY_DIR" ]
then
   rdir=$RAWCOPY_DIR
else
   rdir='/data/spectro/rawcopy'
fi

if [ -n "$ASTROLOG_DIR" ]
then
   speclog_dir=$ASTROLOG_DIR
else
   speclog_dir='/data/spectro/astrolog'
fi


data=`ls -d $rawdata_dir/[56789]???? | tail -n 1`

offsite=schlegel@spectro.princeton.edu
target_data=`ssh ''$offsite'' echo '$RAWCOPY_DIR'`
target_log=`ssh ''$offsite'' echo '$ASTROCOPY_DIR'`

fermi=sburles@fsgi03.fnal.gov
fermi_data=`ssh -a -x -v ''$fermi'' echo '$RAWCOPY_DIR'`
fermi_log=`ssh -a -x -v ''$fermi'' echo '$ASTROCOPY_DIR'`

echo $target_data $target_log
echo Fermi: $fermi $fermi_data $fermi_log
#------------------------------------------------------------------------------
# First find all /data/spectro directories which are missing

for dir in $data
do

    echo $rawdata_dir $speclog_dir
    mjd=`echo $dir | sed -n 's/\/.*\///p'`
    astrolog=$speclog_dir/$mjd
    echo $mjd $astrolog $dir

#
#   First move to rawcopy locally...	
#   Exclude plain .fit files, assuming that .gz versions exist
#
    rsync -arv --exclude="*.fit" $dir $rdir


#######################################################
#       Now move offsite, transfer sdReport last...
#
    rsync -arv --rsh="ssh -c blowfish" $rdir/$mjd $offsite:$target_data

    rsync -arv --rsh="ssh -c blowfish" \
       --exclude="*sdReport*" $astrolog $offsite:$target_log

#
#	Second round, just to be sure everything transferred
#
    rsync -arv --exclude="*.fit" $dir $rdir

    rsync -arv --rsh="ssh -c blowfish" $rdir/$mjd $offsite:$target_data

    rsync -arv --rsh="ssh -c blowfish" $astrolog $offsite:$target_log

    #############################################################
    #########  Now ready to move to other offsites ##############

    exist=`ssh ''$fermi'' ls -d ''$fermi_data''/''$mjd'' 2>/dev/null`
    echo $exist
    if [ ! $exist ]
    then
      scp -rv $offsite:$target_data/$mjd $fermi:$fermi_data/
      scp -rv $offsite:$target_log/$mjd  $fermi:$fermi_log/
    fi
    
done

