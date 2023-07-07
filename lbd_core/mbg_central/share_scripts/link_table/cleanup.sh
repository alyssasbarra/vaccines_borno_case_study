#!/bin/bash

if [[ -z $base_shapefile_dir || -z $shape_dir || -z $temp_dir ]]; then
    echo "Missing environment variables from launch script."
    exit 1
fi

LT_FILENAME=lbd_standard_link.rds

# if files are the same, remove temporary directory
if cmp -s "$temp_dir/$LT_FILENAME" "$base_shapefile_dir/$shape_dir/$LT_FILENAME";
then
   rm -r "$temp_dir"
fi
