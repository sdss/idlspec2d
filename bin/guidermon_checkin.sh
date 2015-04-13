#!/bin/bash
#
# Extract the night's guider data (seeing, etc.) from the processed guider 
# files and put the result in a .par file based on the MJD. Then check
# that file in to the speclog product.
#
# Must be run as eboss@sdss4-eboss, designed for crontab.
#
# Requires ~eboss/bin/sjd.py and guidermonfile.pro from idlspec2d.

function run_and_test {
    "$@"
    status=$?
    if [ $status -ne 0 ]; then
        # useful to have the full error on both stderr and stdout
        echo "Error with $@" | tee /dev/stderr
    fi
    return $status
}

# Need to be able to find the ssh agent in order for svn checkins to work.
export SSH_AUTH_SOCK=/home/eboss/sos/control/agent.socket

# cronjobs need the idlspec2d product
source /home/sdss4/products/eups/bin/setups.sh
setup idlspec2d trunk

# for password-less ssh
export SVN_SSH="ssh -i /home/eboss/.ssh/id_dsa-sos"

if [ $# == 0 ] ; then
    export MJD=`/home/eboss/bin/sjd.py`
else
    export MJD=$1
fi

export GUIDE_DIR=/data/gcam/
export SVN_MESSAGE="committing guiderMon for $MJD"

# guidermonfile writes the output .par to $SPECLOG_DIR/$MJD
echo "Running: guidermonfile for $MJD"
# Sadly, IDL sends its startup junk to stderr, so we have to redirect that
# and can't reliably watch IDL's stderr for actual error messages.
{
idl -e "guidermonfile, mjd=getenv('MJD')"
} 2>&1

echo "SVN commit message:" $SVN_MESSAGE
cd $SPECLOG_DIR/$MJD
run_and_test /usr/bin/svn add guiderMon-$MJD.par
run_and_test /usr/bin/svn commit -m "$SVN_MESSAGE" guiderMon-$MJD.par

exit 0
