#!/usr/bin/env bash
set -euo pipefail

usage() {
    local exec_name
    exec_name=$(basename "$0")
    cat >&2 <<EOF
usage: $exec_name module "command"
EOF
    exit 1
}

while getopts ":h" flag; do
    case "$flag" in
        h) usage ;;
        \?) usage ;;
    esac
done

shift $((OPTIND - 1))

[[ $# -ge 2 ]] || usage

module_name=$1
shift
command=$*

module purge
module load "$module_name"
module list

JDATE=$(python -c 'from boss_drp.utils import jdate; print(str(jdate.astype(str)))' || true)
JDATE=${JDATE:-0}
export JDATE
export MODULE="$module_name"

EMAIL_FILE="${BOSS_DRP_DAILY_DIR:?BOSS_DRP_DAILY_DIR not set}/etc/emails"
EMAIL_RECIPIENTS=()
if [[ -f "$EMAIL_FILE" ]]; then
    mapfile -t EMAIL_RECIPIENTS < <(grep -Ev '^\s*($|#)' "$EMAIL_FILE" || true)
else
    echo "Warning: Email list file not found: $EMAIL_FILE" >&2
fi

TIMEOUT_DURATION=$((48 * 60 * 60))
echo "Running command with a ${TIMEOUT_DURATION}-second timeout (48 hours)..."

if timeout "$TIMEOUT_DURATION" bash -lc "$command"; then
    exit 0
else
    status=$?
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    hostname=$(hostname)
    subject="[cronrun.bash] Job Alert: ${MODULE} failed on ${hostname}"
    message="Job started at: ${timestamp}
Module: ${MODULE}
Command: ${command}
Exit code: ${status}"

    if [[ $status -eq 124 ]]; then
        message+=$'\n\nReason: Command timed out after 48 hours.'
        subject="[cronrun.bash] TIMEOUT: ${MODULE} job on ${hostname}"
        echo "Error: Command timed out after 48 hours." >&2
    else
        message+=$'\n\nReason: Command failed with exit code '"$status"$'.'
        echo "Error: Command failed with exit code ${status}." >&2
    fi

    if ((${#EMAIL_RECIPIENTS[@]})); then
        printf '%s\n' "$message" | mail -s "$subject" "${EMAIL_RECIPIENTS[@]}"
    else
        echo "No email recipients found, skipping notification." >&2
    fi

    exit "$status"
fi