#!/bin/bash
# Interface with R output.

# export delimiter for use in R scripts
export r_output_delimiter=";"

# Parse R output.
#
# This function gets only desired output.
get_r_output () {
    echo "$1" | grep -o "\\[1] \"$r_output_delimiter.*" | cut -d '"' -f2 | cut -d "$r_output_delimiter" -f2
}
