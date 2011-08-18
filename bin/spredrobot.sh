#! /bin/sh
#------------------------------------------------------------------------------
# Script to reduce new spectro data, as reported by the rsync robot.
#
# Once the data is copied over, build plan files under the directories
#    $BOSS_SPECTRO_REDUX/$PLATE
# then launch "batch2d" and "batch1d" to reduce them.  Run the same version
# of the idlspec2d product as is set up for this script.  For example, even
# if v4_9 is declared current, run v4_8 if that is what's set up when this
# script is run.  So, you should run "sprobot_start" from the version of
# the code that you want to reduce the data.
#
# D. Schlegel, Princeton, 19 Dec 2000
#------------------------------------------------------------------------------
# Set file names.

astrologdir=$ASTROLOG_DIR
toprawdir=$BOSS_SPECTRO_DATA
topoutdir=$BOSS_SPECTRO_REDUX

copiedMJDs=$toprawdir/copiedMJDs.list
plannedMJDs=$topoutdir/plannedMJDs.list

#------------------------------------------------------------------------------
# Test that certain environment variables are already set.

if [ -z "$BOSS_SPECTRO_DATA" ] ; then
  echo "BOSS_SPECTRO_DATA must be set!"
  exit
fi

#------------------------------------------------------------------------------
# If $topoutdir is not defined, then exit.

if [ -z "$topoutdir" ] ; then
   echo "SPREDROBOT: TOPOUTDIR must be set!"
   exit
fi

echo ""
echo "-------------------------------------------------------------------------------"
echo "SPREDROBOT: Launched at "`date` UID=$UID PPID=$PPID
echo "IDLSPEC2D_DIR="$IDLSPEC2D_DIR
echo "IDLUTILS_DIR="$IDLUTILS_DIR

#------------------------------------------------------------------------------
# Find copied mjds which have not been reduced

mjds=`fgrep -v -f $plannedMJDs $copiedMJDs | sort -n | uniq`
mjdlist=`echo $mjds | perl -ane 'print join(",",split())'`

#------------------------------------------------------------------------------
# Build the plan files if there is new data.
# We might be called fairly often, so short circuit if we can

if [ ! "$mjds" ] ; then
    echo "SPREDROBOT: Exited at "`date` UID=$UID PPID=$PPID
    exit 0
fi

echo SPREDROBOT: MJDLIST=$mjdlist
echo ""
echo SPREDROBOT: "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]"
echo "spplan2d, topoutdir='$topoutdir', mjd=["$mjdlist"]" | idl
echo ""
echo SPREDROBOT: "spplan1d, topindir='$topoutdir', mjd=["$mjdlist"]"
echo "spplan1d, topindir='$topoutdir', mjd=["$mjdlist"]" | idl


# Register that we have consumed some mjds.
for m in $mjds; do
    echo $m >> $plannedMJDs
done

#------------------------------------------------------------------------------
# Decide what the current version of idlspec2d is, and run that version
# if it is a UPS-declared version.

vers=`echo "print,idlspec2d_version()" | idl 2> /dev/null`
if [ ${vers:0:5} != NOCVS ] ; then upsversion=$vers ; fi

#------------------------------------------------------------------------------
# Batch process 2D first, wait for it to complete, then batch process 1D
# in the background.  The calls to the 2d and 1d scripts will exit if
# those scripts are already running.

   echo "SPREDROBOT: Current time "`date` UID=$UID PPID=$PPID
   sprobot2d.sh ",topdir='$topoutdir', upsversion='$upsversion', nice=19"
   cd $BOSS_SPECTRO_REDUX
   echo "platelist, /create" | idl
   sprobot1d.sh ",topdir='$topoutdir', upsversion='$upsversion', nice=19" &

#------------------------------------------------------------------------------

echo "SPREDROBOT: Finished at "`date` UID=$UID PPID=$PPID

exit
#------------------------------------------------------------------------------
