#!/usr/bin/env bash
set -euo pipefail

usage() {
    local exec_name
    exec_name=$(basename "$0")

    cat <<EOF >&2
Usage: $exec_name module [options]

Description:
    Load the correct module and execute the QA plotting script.

Options:
    -l          Use LCO observations (default is APO).
    -c          Include the --clobber_lists option.
    -n          Disable linking (default is linking enabled).
    -e          Include the --epoch option.
    -w          Generate HTML output.
    -u NAME     HTML output name (used with -w).
    -h          Display this help message and exit.

Example:
    $exec_name myModule -l -c -n -e -w -u test.html
EOF
    exit 1
}

[[ $# -ge 1 ]] || usage

mod=$1
shift

[[ $mod == "-h" ]] && usage

lco=""
tests="-t False"
obs="APO"
clobber=""
nolink=0
epoch=""
html=0
html_name=""

while getopts ":lcnewu:h" flag; do
    case "$flag" in
        l) lco="--lco"; obs="LCO" ;;
        c) clobber="--clobber_lists" ;;
        n) nolink=1 ;;
        e) epoch="--epoch" ;;
        w) html=1 ;;
        u) html_name=$OPTARG ;;
        h) usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
        \?) usage ;;
    esac
done

module purge
module load "$mod"
module list

: "${BOSS_QA_DIR:?BOSS_QA_DIR is not set or is empty}"
: "${BOSS_DRP_DAILY_DIR:?BOSS_DRP_DAILY_DIR is not set or is empty}"
: "${BOSS_SPECTRO_REDUX:?BOSS_SPECTRO_REDUX is not set or is empty}"
: "${RUN2D:?RUN2D is not set or is empty}"

EMAIL_FILE="$BOSS_DRP_DAILY_DIR/etc/emails"
EMAIL_RECIPIENTS=()

if [[ -f "$EMAIL_FILE" ]]; then
    mapfile -t EMAIL_RECIPIENTS < <(grep -Ev '^\s*($|#)' "$EMAIL_FILE" || true)
else
    echo "Warning: Email list file not found: $EMAIL_FILE" >&2
fi

TIMEOUT_DURATION=$((48 * 60 * 60))
echo "Running `boss_drp run Plot_QA` block with a ${TIMEOUT_DURATION}-second timeout (48 hours)..."

run_plot() {
    set -euo pipefail

    boss_drp run Plot_QA --run2d "$RUN2D" $tests $lco $clobber $epoch --cron

    if [[ "$html" -eq 1 ]]; then
        local html_out
        html_out=${html_name:-cronplot_QA.html}
        boss_drp run Plot_QA --run2d "$RUN2D" $tests $lco $clobber $epoch --cron --html "$html_out"
    fi

    if [[ "$nolink" -eq 0 ]]; then
        rm -f "${BOSS_QA_DIR}/QA_${obs}.png"
        ln -s "${BOSS_SPECTRO_REDUX}/${RUN2D}/spCalib_QA-${RUN2D}-${obs}.png" \
              "${BOSS_QA_DIR}/QA_${obs}.png"

        rm -f "${BOSS_QA_DIR}/SN2_${obs}.png"
        ln -s "${BOSS_SPECTRO_REDUX}/${RUN2D}/SN2-${RUN2D}-${obs}.png" \
              "${BOSS_QA_DIR}/SN2_${obs}.png"
    fi
}

export RUN2D tests lco clobber epoch html html_name nolink obs BOSS_QA_DIR BOSS_SPECTRO_REDUX
export -f run_plot

if timeout "$TIMEOUT_DURATION" bash -lc run_plot; then
    exit 0
else
    status=$?
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    hostname=$(hostname)
    subject="[cronplot_QA.bash] Job Alert: cronplot_QA failed on ${hostname}"
    message="Job started at: ${timestamp}
Module: ${mod}
Exit code: ${status}"

    if [[ $status -eq 124 ]]; then
        message+=$'\n\nReason: Block timed out after 48 hours.'
        subject="[cronplot_QA.bash] TIMEOUT: ${mod} job on ${hostname}"
        echo "Error: Block timed out after 48 hours." >&2
    else
        message+=$'\n\nReason: Command in block failed with exit code '"$status"$'.'
        echo "Error: Block failed with exit code ${status}." >&2
    fi

    if ((${#EMAIL_RECIPIENTS[@]})); then
        printf '%s\n' "$message" | mail -s "$subject" "${EMAIL_RECIPIENTS[@]}"
    else
        echo "No email recipients found, skipping notification." >&2
    fi

    exit "$status"
fi