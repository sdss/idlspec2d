#! /bin/sh

dest=/data/spectro/logfiles/
lock=/data/spectro/maillock

if [ ! -r $lock ] 
then
   echo "Lock file does not exists"
   exit
fi


logs=`ls /data/spectro/5*/* | grep logfile | grep html | grep -v lock`

for thislog in $logs 
do
  if newer $thislog $lock
  then 
  
    dir=`echo $thislog | sed -n 's/\/[^\/]*$//p'`
    mailfile=`echo $thislog | sed -n 's/logfile/mail/p'`
    echo $thislog $dir $mailfile

    subject=`grep LOGSHEET $thislog | grep TITLE | sed -n 's/<.[A-Z]*>//pg'`
    echo $subject
    echo "</pre>" > $mailfile
    cat $thislog | grep -v HTML | grep -v HEAD >> $mailfile
    echo "<pre>" >> $mailfile
    
    sn=`ls $dir | grep snplot | grep ps`
    for thissn in $sn
    do 
       echo $thissn
       echo "!$thissn <<EOT" >> $mailfile
       cat $dir/$thissn >> $mailfile
       echo "EOT" >> $mailfile
    done 

    mail -s "$subject" sdss-test@astro.princeton.edu < $mailfile
  fi
done


#echo "zero point for mailing"  > $lock
