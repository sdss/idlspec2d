#!/bin/bash

thistime="eBOSS cron autocommit on "`date`

echo $SPECLOG_DIR
if [ -z "$SPECLOG_DIR" ] ; then
  export SPECLOG_DIR='/home/eboss/software/svn.sdss.org/data/sdss/speclog/trunk'
fi
echo $SPECLOG_DIR
cd $SPECLOG_DIR
tmp_mjd=$(ls -ltr $SPECLOG_DIR | grep '^d' | tail -1)
mjd=$(echo $tmp_mjd | rev | cut  -d' ' -f1 |rev)

svn update
svn add $mjd
svn commit -m  "$thistime"
