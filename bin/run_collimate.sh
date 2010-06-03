#!/bin/bash

# Test that an exposure number is given as an argument
if [ $# = 0 ] ; then
  echo "Must specify exposure number"
  exit
fi

# Launch all 4 cameras as separate processes, then wait for completion
echo "collimate,$1,docams='b1'" | idl 2> /dev/null &
echo "collimate,$1,docams='b2'" | idl 2> /dev/null &
echo "collimate,$1,docams='r1'" | idl 2> /dev/null &
echo "collimate,$1,docams='r2'" | idl 2> /dev/null &
wait

# If the log files exist, then print the last lines of each
files=`find . -name "Collimate-*$1.log" -print`

  echo ""
  echo ""
  echo ""
  echo ""

if [ -n "$files" ] ; then
  tail -n 12 Collimate-*$1.log
else
  echo "Nothing computed"
fi
exit

