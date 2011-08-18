#! /bin/sh

#------------------------------------------------------------------------------
# Remove until 95% full or less
#
# S. Burles, APO, 7 June 2001
#------------------------------------------------------------------------------

if [ -n "$BOSS_SPECTRO_DATA" ]
then
   rawdata_dir=$BOSS_SPECTRO_DATA
fi

if [ -n "$ASTROLOG_DIR" ]
then 
   speclog_dir=$ASTROLOG_DIR
fi

if [ $1 ]
then 
  rawdata_dir=$1
fi

if [ $2 ]
then 
  speclog_dir=$2
fi

#  Quit if rawdata is not specified

if [ ! $rawdata_dir ] 
then
  exit
fi

destdisk=`df $rawdata_dir | \
       awk '{if (FNR>1 && $5 * 1>95) {print $5 $6}}' | head -1`

echo $destdisk

if [ $destdisk ] 
then

#  first get last directory

    dir=`ls $rawdata_dir | grep -e "[5-9][0-9][0-9][0-9][0-9]" | sort | head -1`
    echo $dir
  
    echo rm -rf $rawdata_dir/$dir
    rm -rf $rawdata_dir/$dir

    if [ $speclog_dir ]
    then 
      echo rm -rf $speclog_dir/$dir
      rm -rf $speclog_dir/$dir
    fi

#   Recursively call itself to check again

    $0 $1 $2
    
fi



