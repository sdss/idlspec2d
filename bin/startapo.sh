#! /bin/sh
#
#  Try to parse 1 filename sent from stdin
#  Break up into path and simple filename
#  Guess astrolog directory and send to idl
#


input=`echo $1 | sed -n 'p'`
echo Processing ... $input 


dir=`echo $input | sed -n 's/\/[^\/]*$//p'`
filename=`echo $input | sed -n 's/\/.*\///p'`
echo Dir: $dir
echo file: $filename
astrolog=`echo $dir | sed -n 's/\/data\/spectro/\/astrolog/p'`
echo astro: $astrolog


good=`expr "$filename" : 'sdR'` 
echo $good
if [ "$good" -gt 0 ] 
then
#echo  "  echo \"aporeduce, '$filename',indir='$dir', plugmapdir='$astrolog'\" "
echo Processing ... $input 
echo  "aporeduce, '$filename',indir='$dir', plugdir='$astrolog' " 
#echo  "aporeduce, '$filename',indir='$dir', plugdir='$astrolog' " | idl
fi

exit
