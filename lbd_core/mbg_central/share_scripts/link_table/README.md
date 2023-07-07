# Link Table Generation in Parallel

## Overview

This script generates the link table in parallel. Most of the parameters are set in `launch.sh`.

The following pattern is followed:

1. Partition the shapefile by ADM0 geometry area.
2. Build a link table for each partition.
3. Combine the partitioned link tables together.
4. Cleanup/ensure script ran successfully.
5. Validate link table.

## Use

    ./launch.sh -d <directory> -n <number of partitions> -q <queue> -p <project>

Where `directory` is a child directory within the base admin shapefile directory.

Example:

    ./launch.sh -d DATE -n 5000

will generate a link table for FILEPATH and place the new link table and id raster there.
