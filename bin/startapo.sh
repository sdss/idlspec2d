#! /bin/sh
#------------------------------------------------------------------------------
# The script "aporsync" calls this script.
# Try to parse 1 filename sent from stdin.
# The file should look something like 
#    /data/spectro/51666/sdR-b1-00001234.fit
# Break this up into path and simple filename.
# Guess astrolog directory, replacing "/data/spectro" with "/astrolog".
# Send all this to the IDL routine APOREDUCE.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

#
# Wait for first file to finish complete copy
#
sleep 5

# Loop through each file name that has been passed.
while read f
do 

  input=$f
  dir=`echo $input | sed -n 's/\/[^\/]*$//p'`
  filename=`echo $input | sed -n 's/\/.*\///p'`
  echo STARTAPO: Directory $dir
  echo STARTAPO: Filename $filename
  astrolog=`echo $dir | sed -n 's/\/data\/spectro/\/astrolog/p'`
  outdir=`echo $dir | sed -n 's/spectro/spectrologs/p'`
  copydir=/data/spectrologs/html/

# If we don't have the output directories, make them.
  if [ ! -d $outdir ] 
  then
    mkdir $outdir
  fi
  if [ ! -d $copydir ] 
  then
    mkdir $copydir
  fi

  good=`expr "$filename" : 'sdR'`
  if [ "$good" -gt 0 ] 
    then
     echo STARTAPO: Processing $input at `date`
     echo "aporeduce, '$filename',indir='$dir', outdir='$outdir', \
          plugdir='$astrolog', copydir='$copydir' " | nice idl >& /dev/null
  fi

done

exit
