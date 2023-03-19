#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

PATH="/software/singularity-v3.6.4/bin/":$PATH
MOUNT_POINT='/lustre/scratch117/cellgen/team283/tl10/NXF_WORK/'

DATE_WITH_TIME=`date "+%Y%m%d%H%M"`
TMP_NF_WORK="$MOUNT_POINT/${USER}_${DATE_WITH_TIME}_segmentation_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK LSB_DEFAULTGROUP='team283' /lustre/scratch117/cellgen/team283/tl10/nextflow/nextflow -trace nextflow.executor run /lustre/scratch117/cellgen/team283/tl10/workflow-segmentation/main.nf \
	-params-file $1 \
	-profile lsf \
	-entry run_nuc_seg \
	-resume
