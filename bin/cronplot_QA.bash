#!/usr/bin/env bash
#
# cronplot_QA.bash
#
# Designed to load the correct module and execute the QA plotting script
#
# usage: cronplot_QA.bash module options
#
# Created by Sean Morrison on 13 Aug 2024

function usage {
    local execName=$(basename "$0")
    cat <<EOF >&2
Usage: $execName [options]

Description:
    Load the correct module and execute the QA plotting script.

Options:
    -m module   Specify the module to load.
    -l          Use LCO observations (default is APO).
    -c          Include the --clobber_lists option.
    -n          Disable linking (default is False).
    -e          Include the --epoch option.
    -w          Generate HTML output (default is False).
    -h          Display this help message and exit.

Example:
    $execName -m myModule -l -c -n -e -w "test.html"

EOF
    exit 1
}

lco=""
tests="-t False"
obs='APO'
clobber=''
nolink='F'
epoch=''
html=''
html_name=''
while getopts m:lcnewhu: flag; do
    case "${flag}" in
        m) mod=${OPTARG} ;;
        l)
            lco="--lco"
            obs='LCO'
            ;;
        c) clobber='--clobber_lists' ;;
        n) nolink='T' ;;
        e) epoch="--epoch" ;;
        w) html="T" ;;
        u) html_name="--html_name "${OPTARG} ;;
        h) usage ;;
        *) usage ;;  # Catch invalid options
    esac
done

# Check if BOSS_QA_DIR is set
if [ -z "$BOSS_QA_DIR" ]; then
    echo "BOSS_QA_DIR is not set or is empty. Exiting."
    exit 1
fi

# Load the specified module
module purge
module load "$mod"
module list

# Check if in test mode #this is not needed since it is part of BOSS_SPECTRO_REDUX
#if [[ "$BOSS_SPECTRO_REDUX" == */test/* ]]; then
#    tests=" -t True"
#fi

# Run the plot QA script
plot_qa --run2d "$RUN2D" $tests $lco $clobber $epoch --cron

# Handle linking and HTML output
if [[ $nolink == 'F' ]]; then
    if [[ $html == 'T' ]]; then
        plot_qa --run2d "$RUN2D" $tests $lco $clobber $epoch --cron --html $html_name
    fi
    #rm -f "${BOSS_QA_DIR}/QA_$obs.png"
    #ln -s "${BOSS_SPECTRO_REDUX}/$RUN2D/spCalib_QA-$RUN2D-$obs.png" "${BOSS_QA_DIR}/QA_$obs.png"
    #rm -f "${BOSS_QA_DIR}/SN2_$obs.png"
    #ln -s "${BOSS_SPECTRO_REDUX}/$RUN2D/SN2-$RUN2D-$obs.png" "${BOSS_QA_DIR}/SN2_$obs.png"
fi