#! /bin/sh
#------------------------------------------------------------------------------
# This is a cron job that should run once per day from plate-mapper.apo.nmsu.edu,
# currently at 7am.
#
# Look for all MJD directories, "/data/spectro/spectrologs/5*".  Loop through
# each such directory.  If a file "logfile*html" exists, then construct a
# message to send to the SDSS mailing list sdss-test (???).
# This message is a 1-line text message that links to the HTML file, which in
# turn links to any PostScript plots that were also in that same directory.
#
# S. Burles, APO, 4 May 2000
#------------------------------------------------------------------------------

logs=`find /data/spectro/spectrologs/5* -name "logfile*html" -print | grep -v current`

#   This doesn't work below below the argument list gets too large
# logs=`ls -d /data/spectro/spectrologs/5*/* \
#          | grep logfile | grep html | grep -v lock`

# The variable $thislog is the name of the HTML file with its full path.
# The variable $filename has the path stripped off.

for thislog in $logs 
do
    dir=`echo $thislog | sed -n 's/\/[^\/]*$//p'`
    filename=`echo $thislog | sed -n 's/\/.*\///p'`
    mailfile=`echo $thislog | sed -n 's/logfile/mail/p'`

    subject=`grep LOGSHEET $thislog | grep TITLE | sed -n 's/<.[A-Z]*>//pg'`
    echo $subject
    echo "" > $mailfile
    echo '<A HREF="'$filename'">'$filename'</A>' >> $mailfile
    echo "" >> $mailfile

    echo "!$filename<<EOT" >> $mailfile
    cat $thislog | sed -e 's/<BODY.*>/<BODY>/' >> $mailfile
    echo "EOT" >> $mailfile
   
    sn=`find $dir -name "snplot*ps" -print` 
#    sn=`ls $dir | grep snplot | grep ps`
    for thissn in $sn
    do 
       echo $thissn
       snname=`echo $thissn | sed -n 's/\/.*\///p'`
       echo "!$snname <<EOT" >> $mailfile
       cat $thissn >> $mailfile
       echo "EOT" >> $mailfile
    done 

    mail -s "$subject" sdss-speclog@astro.princeton.edu < $mailfile
#    mail -s "$subject" sdss-test@astro.princeton.edu < $mailfile

#
#	Kill almost everything in the log directory
#

    rm -f $dir'/wset*.fits'
    rm -f $dir'/tset*.fits'
    rm -f $dir'/fflat*.fits'
    rm -f $dir'/sci*.fits'

done

killdata.sh

