#!/bin/bash

thistime="eBOSS cron autocommit on "`date`

if [ -z "$SPECLOG_DIR" ] ; then
  export SPECLOG_DIR='/home/sdss4/products/NULL/speclog/trunk'
fi

cd $SPECLOG_DIR
tmp_mjd=$(ls -ltr $SPECLOG_DIR | grep '^d' | tail -1)
mjd=$(echo $tmp_mjd | rev | cut  -d' ' -f1 |rev)

svn update
svn add $mjd
svn commit -m  "$thistime"
