#!/bin/bash

# This script runs a final cleanup (R) script and notifies the user via email
#
# Arguments:
#
# $1 script_str - R script string to run and retrieve the output for
# $2 email_subject - email subject line

echo "$HERE/r_interface.sh"

if [[ ! -f "$HERE/r_interface.sh" ]]; then
    echo "Could not find r_interface file. Exiting."
    exit 1
fi

source "$HERE/r_interface.sh"

script_str=$1
email_subject=$2

# run cleanup R script
script_output=$($script_str)
msg=$(get_r_output "$script_output")

# email user (this requires being run in subshell)
$(echo "$msg" | mail -s "$email_subject" "USERNAME")
