#!/bin/bash
#
# Extract the night's guider data (seeing, etc.) from the processed guider 
# files and put the result in a .par file based on the MJD. Then check
# that file in to the speclog product.
#
# Must be run as observer@sos3, designed for crontab.
#
# Requires ~observer/bin/sjd.py and guidermonfile.pro from idlspec2d.

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
export SSH_AUTH_SOCK=/home/observer/sos/control/agent.socket

# cronjobs need the idlspec2d product
source /home/sdss3/products/eups/bin/setups.sh
setup idlspec2d

# for password-less ssh
export SVN_SSH="ssh -i /home/observer/.ssh/id_dsa-sos"

export GUIDE_DIR=/data/gcam/
export MJD=`/home/observer/bin/sjd.py`
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
run_and_test /usr/local/bin/svn add guiderMon-$MJD.par
run_and_test /usr/local/bin/svn commit -m "$SVN_MESSAGE" guiderMon-$MJD.par

exit 0
