#! /bin/sh

if ps -ax | grep aporsync.sh  | grep -v -e grep
then
  echo "aporsync already running, exiting..."
  exit
fi

aporsync.sh

