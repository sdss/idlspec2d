#! /bin/sh
#----------------------------------------------------------------------------------
# This is a deadly script to blow away all MJDs in /astrolog and /data/spectro
#  which do not exist on host.  It is invoked from mailhtml, and called
# every morning
#
# S. Burles, APO, 4 May 2000
#----------------------------------------------------------------------------------

data=`ls -d /data/spectro/5????`

echo $data

#first find all /data/spectro directories which are missing

for dead in $data
do

    if  [ `ssh sdsshost ls -d $dead 2>/dev/null` ]   
    then
      echo $dead still exists
    else
      echo I should kill $dead
      rm -rf $dead
    fi

done

data=`ls -d /astrolog/5????`

echo $data

#now find all /astrolog/5 directories which are missing

for dead in $data
do

    if  [ `ssh sdsshost ls -d $dead 2>/dev/null` ]   
    then
      echo $dead still exists
    else
      echo I should kill $dead
      rm -rf $dead
    fi

done



