#! /bin/sh
#------------------------------------------------------------------------------
# The scripts "aporsync_blue.sh" and "aporsync_red.sh" call this script.
# Try to parse 1 or many filename(s) sent from stdin.
# The file should look something like 
#    /data/spectro/51666/sdR-b1-00001234.fit
# Break this up into path and simple filename.
# Guess astrolog directory, replacing "/data/spectro" with "/astrolog".
# Send all this to the IDL routine APOREDUCE.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

if [ -n "$ASTROLOG_DIR" ]
then 
   astrolog_dir=$ASTROLOG_DIR
else
   astrolog_dir='/data/spectro/astrolog'
fi

if [ -n "$SPECTROLOG_DIR" ]
then 
   spectrolog_dir=$SPECTROLOG_DIR
else
   spectrolog_dir='/data/spectro/spectrologs'
fi

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
  copydir=/data/spectro/spectrologs/html/

  if [ ! -d $copydir ] 
  then
    mkdir $copydir
  fi

  good=`expr "$filename" : 'sdR'`
  if [ "$good" -gt 0 ] 
    then
     mjd=`echo $dir | sed -n 's/\/.*\///p'`
     astrolog=$astrolog_dir/$mjd
     outdir=$spectrolog_dir/$mjd

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

### GZIP in same location

  good=`echo $filename | sed -n 's/fit//p'`
  if [ $good ]
  then 
    gzip -c $input > $input.gz &
    chmod 664 $input.gz
  fi


done

exit
