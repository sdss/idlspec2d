#! /bin/sh
#------------------------------------------------------------------------------
# This is a cron job that should run once per day from
# plate-mapper.apo.nmsu.edu, currently at 7:20 am.
#
# Look for all MJD directories, "$SPECTROLOG_DIR/[56789]*".
# Loop through each such directory.  If a file "logfile*html" exists, then
# construct a message to send to the SDSS mailing list at:
#    sdss-speclog.princeton.edu
# This message is a 1-line text message that links to the HTML file, which in
# turn links to any PostScript plots that were also in that same directory.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

if [ -z "$SPECTROLOG_DIR" ] ; then
   echo "Abort: SPECTROLOG_DIR not set!"
   exit
fi
SOSLOG=$SPECTROLOG_DIR/sos.log

##
#  First, send data offsite, and then create mail messages.
##  (Don't do this now, everybody gets the data themselves)
##  (Need to at least move to rawcopy)
data_rsync.sh

logs=`find $SPECTROLOG_DIR/[56789]* -name "logfile*html" -print | grep -v current`

# The variable $thislog is the name of the HTML file with its full path.
# The variable $filename has the path stripped off.

for thislog in $logs 
do
    dir=`echo $thislog | sed -n 's/\/[^\/]*$//p'`
    filename=`echo $thislog | sed -n 's/\/.*\///p'`
    mailfile=`echo $thislog | sed -n 's/logfile/mail/p'`

    if [ ! -f $mailfile ]
    then 
      subject=`grep "APO Spectro" $thislog | grep TITLE | sed -n 's/<.[A-Z]*>//pg'`
      echo $subject
      echo "" > $mailfile
      echo '<A HREF="'$filename'">'$filename'</A>' >> $mailfile
      echo "" >> $mailfile

      echo '!'"$filename<<EOT" >> $mailfile
      cat $thislog | sed -e 's/<BODY.*>/<BODY>/' >> $mailfile
      echo "EOT" >> $mailfile
   
      sn=`find $dir -name "snplot*ps" -print` 
#     sn=`ls $dir | grep snplot | grep ps`
      for thissn in $sn
      do 
        echo MAILHTML: SNFILE $thissn
        snname=`echo $thissn | sed -n 's/\/.*\///p'`
        echo '!'"$snname <<EOT" >> $mailfile
        cat $thissn >> $mailfile
        echo "EOT" >> $mailfile
      done 

      mail -s "$subject" sdss-speclog@astro.princeton.edu < $mailfile

#     Kill almost everything in the log directory
#      rm -f $dir/fflat*.fits
#      rm -f $dir/sci*.fits
#      rm -f $dir/tset*.fits
#      rm -f $dir/wset*.fits
#      rm -f $dir/*.ps
#      rm -f $dir/logfile*.fits
#      rm -f $dir/logfile*.html
   fi
done

killdata.sh

\mv $SOSLOG.6dayold $SOSLOG.7dayold
\mv $SOSLOG.5dayold $SOSLOG.6dayold
\mv $SOSLOG.4dayold $SOSLOG.5dayold
\mv $SOSLOG.3dayold $SOSLOG.4dayold
\mv $SOSLOG.2dayold $SOSLOG.3dayold
\mv $SOSLOG.1dayold $SOSLOG.2dayold
\mv $SOSLOG         $SOSLOG.1dayold

