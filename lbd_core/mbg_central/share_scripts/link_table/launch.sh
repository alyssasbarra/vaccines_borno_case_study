#!/bin/bash
#
# This script contains launch scripts for the various phases of link table generation:
# 1. partition
# 2. build link tables
# 3. finalize/combine link tables
# 4. check that link table was generated correctly
# 5. clean up temporary directory
# 6. validate link table
#
# Args:
#   d: shapefile directory
#   q: queue (default: geospatial.q)
#   p: project (default: proj_geo_nodes)
#   n: target number of partitions (default: 2000)
#
# Notes:
#   You MUST qlogin before running this script.


# source directory
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" || exit ; pwd -P)"

# get r interface script
source "FILEPATH"

# export link table library path for R scripts
export R_LIBS_USER="FILEPATH"

# arguments / defaults
q=geospatial.q
p=proj_geo_nodes
n=2000
while getopts "q::p::n::d:" opt; do
  case $opt in
    q) 
        q=$OPTARG      
        ;;
    p)
        p=$OPTARG
        ;;
    n)
        n=$OPTARG
        ;;
    d)
        shape_dir=$OPTARG
        ;; 
    *)
        echo "Invalid flag: $opt"
        ;; 
esac
done

if [ ! "$shape_dir" ]
then
    echo "Error: directory must be specified (-d). Exiting ..."
    exit 1
fi

if [ "$shape_dir" == "current" ]; then
   echo -e "'current' should not be modified. Please provide a different directory.\\nExiting ..."
   exit 1
fi

echo "Attempting to launch $n target jobs on $q queue"

# r shell script
r_shell=FILEPATH
r_script="FILEPATH"
sing_r=FILEPATH
sing_image=FILEPATH

echo "Getting base shapefile directory using mbg scripts ..."
# get base shapefile directory from mbg scripts
r_output=$("$sing_r" -s "$sing_image" -e s "$r_script" -c shapedir)
echo "$r_output"
base_shapefile_dir=$(get_r_output "$r_output")

if [[ -z $base_shapefile_dir ]]; then
    echo "Could not retrieve base shapefile directory. Exiting ..."
    exit 1
fi

echo "Base shapefile directory: $base_shapefile_dir"

# ensure shapefile exists in given directory
if [[ ! -f "FILEPATH" ]]; then
    echo -e "FILEPATH is not a valid shapefile directory.\\nExiting ..."
    exit 1
fi

# Partition the shapefile
temp_dir=$(mktemp --directory --tmpdir='FILEPATH')

echo "Partioning ..."
"$sing_r" -s "$sing_image" -e s "$r_script" -c partition -n "$n" -d "$temp_dir" -s "$shape_dir"

# infer the number of partitions by finding the filename with the highest number
# find: list all files created in partitioning
# sed: trim filename through underscore (polys_1.rds -> 1.rds)
# sort: sort values
# tail: get last (highest) value
# cut: remove the extension
npartitions=$(find "$temp_dir" | sed 's/^.*_//' | sort -n | tail -1 | cut -d "." -f1)

if [[ -z $npartitions ]]; then
    echo "Something went wrong in partitioning. Exiting."
    exit 1
fi

# Build link tables
job_name="build_lt"
mem="5G"
threads=1
runtime="5:00:00"                                                                                                                                                                                     
output="FILEPATH"
errors="FILEPATH"

echo "Launching array job with $npartitions partitions ..."

echo "Logging output to $output."
echo "Logging errors to $errors."

bjid=$(qsub -terse \
	-N "$job_name" \
	-t "1:$npartitions" \
	-q "$q" \
	-P "$p" \
	-l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
	-v sing_image="$sing_image",R_LIBS_USER="$R_LIBS_USER" \
	-o "$output" \
    -e "$errors" \
	"$r_shell" \
	"$r_script" --args -c build -d "$temp_dir" -p "$npartitions" \
	| cut -d "." -f 1) # cut to remove decimal from job id

echo "Build job $bjid queued."

# Combine link tables
job_name="finalize_lt"
mem="20G"
runtime="1:00:00"
threads=10

fjid=$(qsub -terse \
    -N $job_name \
    -q "$q" \
    -P "$p" \
    -l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
    -hold_jid "$bjid" \
    -v sing_image="$sing_image",R_LIBS_USER="$R_LIBS_USER" \
    -o "$output" \
    -e "$errors" \
    "$r_shell" \
    "$r_script" --args -c finalize -d "$temp_dir" -p "$npartitions")

echo "Finalize job $fjid queued."

# Check if job was successful and cleanup
job_name="check_lt"
mem="10G"
threads=1
runtime="1:00:00"
shell_script="$HERE/email.sh"

chjid=$(qsub -terse \
    -N $job_name \
    -q "$q" \
    -P "$p" \
    -l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
    -hold_jid "$fjid" \
    -v sing_image="$sing_image",HERE="$HERE",R_LIBS_USER="$R_LIBS_USER" \
    -o "$output" \
    -e "$errors" \
    "$shell_script" "$sing_r -s $sing_image -e s $HERE/run.R -c check -p $npartitions -d $temp_dir -s $shape_dir" "Link Table Status")

echo "Check job $chjid queued."

job_name="cleanup_lt"
mem="1G"
runtime="1:00:00"
shell_script="$HERE/cleanup.sh"

cljid=$(qsub -terse \
    -N $job_name \
    -q "$q" \
    -P "$p" \
    -l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
    -hold_jid "$chjid" \
    -v sing_image="$sing_image",base_shapefile_dir="$base_shapefile_dir",shape_dir="$shape_dir",temp_dir="$temp_dir",R_LIBS_USER="$R_LIBS_USER" \
    -o "$output" \
    -e "$errors" \
    "$shell_script")

echo "Cleanup job $cljid queued."

job_name="validate_lt"
mem="50G"
runtime="1:00:00:00"
shell_script="$HERE/email.sh"

vjid=$(qsub -terse \
    -N $job_name \
    -q "$q" \
    -P "$p" \
    -l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
    -hold_jid "$chjid" \
    -v sing_image="$sing_image",HERE="$HERE",R_LIBS_USER="$R_LIBS_USER" \
    -o "$output" \
    -e "$errors" \
    "$shell_script" "$sing_r -s $sing_image -e s $HERE/run.R -c validate -s $shape_dir" "Link Table Validation")

echo "Validate job $vjid queued."
