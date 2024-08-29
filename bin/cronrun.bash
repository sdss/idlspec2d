#!/usr/bin/env bash

# cronrun.bash
#
# Designed to load the correct module and execute the daily commands
#
# usage: cronrun.bash module "script"
#
# Created by Sean Morrison on 2/20/24.

function usage {
    local execName=$(basename $0)
    (
    echo "usage: $execName module 'script'"
    echo " "
    ) >&2
    exit 1
}

# Check for the -h flag
while getopts "h" flag; do
    case "$flag" in
        h) usage ;;
    esac
done

# Extract the module and script from the arguments
ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}

# Ensure both arguments are provided
if [ -z "$ARG1" ] || [ -z "$ARG2" ]; then
    usage
fi

# Load the specified module and execute the script
module purge
module load "$ARG1"
module list

export MODULE="$ARG1"
eval "$ARG2"
