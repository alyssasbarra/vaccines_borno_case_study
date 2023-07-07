#$ -S /bin/sh

export SGE_ROOT=FILEPATH
export SGE_CELL=FILEPATH

/usr/local/UGE-prod/bin/lx-amd64/qsub -e FILEPATH -o FILEPATH -cwd -l mem_free=20G -pe multi_slot 10 -P proj_geo_nodes -l geos_node=TRUE -v sing_image=default -v SET_OMP_THREADS=10 -v SET_MKL_THREADS=1 -N sync_shapefiles FILEPATH FILEPATH