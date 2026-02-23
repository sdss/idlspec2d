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


# Run the command and capture its output
JDATE=$(python -c "from boss_drp.utils import jdate; print(str(jdate.astype(str)))")
# Check if the output is empty
if [ -z "$JDATE" ]; then
  # Set to default value if no output is found
  JDATE="0"
fi

export JDATE="$JDATE"


export MODULE="$ARG1"
#eval "$ARG2"


# Read email list from file (ignore empty lines and comments)
EMAIL_FILE="$BOSS_DRP_DAILY_DIR/etc/emails"
if [[ -f "$EMAIL_FILE" ]]; then
    EMAIL_RECIPIENT=$(grep -Ev '^\s*($|#)' "$EMAIL_FILE" | head -n 1) # first address only
    #EMAIL_RECIPIENTS=$(grep -Ev '^\s*($|#)' "$EMAIL_FILE" | tr '\n' ' ') # all address
else
    echo "Warning: Email list file not found: $EMAIL_FILE" >&2
    EMAIL_RECIPIENTS=""
fi


# Set timeout duration (48 hours = 172800 seconds)
TIMEOUT_DURATION=$((48 * 60 * 60))

echo "Running command with a ${TIMEOUT_DURATION}-second timeout (48 hours)..."

if ! timeout "$TIMEOUT_DURATION" bash -c "eval \"$ARG2\""; then
    STATUS=$?

    TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
    HOSTNAME=$(hostname)
    SUBJECT="[cronrun.bash] Job Alert: ${MODULE} failed on ${HOSTNAME}"
    MESSAGE="Job started at: ${TIMESTAMP}
Module: ${MODULE}
Command: ${ARG2}
Exit code: ${STATUS}"

    if [ "$STATUS" -eq 124 ]; then
            MESSAGE="${MESSAGE}\n\nReason: Command timed out after 48 hours."
            SUBJECT="[cronrun.bash] TIMEOUT: ${MODULE} job on ${HOSTNAME}"
            echo "Error: Command timed out after 48 hours." >&2
    else
            MESSAGE="${MESSAGE}\n\nReason: Command failed with exit code ${STATUS}."
            echo "Error: Command failed with exit code ${STATUS}." >&2
    fi

    # Send email notification to all recipients
    if [[ -n "$EMAIL_RECIPIENTS" ]]; then
        echo -e "$MESSAGE" | mail -s "$SUBJECT" $EMAIL_RECIPIENTS
    else
        echo "No email recipients found, skipping notification." >&2
    fi
    exit "$STATUS"
fi
