#! /bin/sh
#------------------------------------------------------------------------------
# This is a deadly script to blow away all MJDs in /astrolog and /data/spectro
# which do not exist on host.  It is invoked from mailhtml, and called
# every morning
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

data=`ls -d /data/spectro/[56789]????`

#------------------------------------------------------------------------------
# First find all /data/spectro directories which are missing

for dead in $data
do

    if  [ `ssh sdsshost ls -d $dead 2>/dev/null` ]   
    then
      echo KILLDATA: $dead still exists
    else
      echo KILLDATA: I should kill $dead
      rm -rf $dead
    fi

done

#------------------------------------------------------------------------------
# Now find all /astrolog directories which are missing
#
#data=`ls -d /astrolog/[56789]????`
#
#for dead in $data
#do
#    if  [ `ssh sdsshost ls -d $dead 2>/dev/null` ]   
#    then
#      echo KILLDATA: $dead still exists
#    else
#      echo KILLDATA: I should kill $dead
#      rm -rf $dead
#    fi
#done
#
