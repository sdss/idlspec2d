#! /bin/sh
#------------------------------------------------------------------------------
# The scripts "aporsync_blue.sh" and "aporsync_red.sh" call this script.
# Try to parse 1 or many filename(s) sent from stdin.
# The file should look something like 
#    /data/spectro/51666/sdR-b1-00001234.fit
# Break this up into path and simple filename.
# Send the file name and astrolog directory to the IDL routine APOREDUCE.
# Finally, make a g-zipped copy of the input file, leaving the uncompressed
#  file there too.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

if [ -n "$ASTROLOG_DIR" ] ; then
   echo "Abort: ASTROLOG_DIR not set!"
   exit
fi

if [ -n "$SPECTROLOG_DIR" ] ; then
   echo "Abort: SPECTROLOG_DIR not set!"
   exit
fi

copydir=$SPECTROLOG_DIR/html

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

  echo STARTAPO: Directory $dir Filename $filename

  if [ ! -d $copydir ] 
  then
    mkdir $copydir
  fi

  good=`expr "$filename" : 'sdR'`
  if [ "$good" -gt 0 ] 
    then
     mjd=`echo $dir | sed -n 's/\/.*\///p'`
     astrolog=$ASTROLOG_DIR/$mjd
     outdir=$SPECTROLOG_DIR/$mjd

     if [ ! -d $outdir ] 
     then
       mkdir $outdir
     fi

     echo STARTAPO: Processing $input at `date`
#     $IDL_DIR/bin/lmutil lmstat
     echo "aporeduce, '$filename',indir='$dir', outdir='$outdir', \
          plugdir='$astrolog', copydir='$copydir' " | nice idl >& $outdir/err.$filename
#          plugdir='$astrolog', copydir='$copydir' " | nice idl >& /dev/null
  fi

  # If the file ends in exactly ".fit", then compress it with gzip,
  # and leave the uncompressed file there too.
  good=`echo $filename | sed -n 's/fit//p'`
  if [ $good ]
  then 
    gzip -c $input > $input.gz &
    chmod 664 $input.gz
  fi

done

exit
