#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/lustre/scratch117/cellgen/team283/tl10/NXF_WORK/'

TMP_NF_WORK="$MOUNT_POINT/${USER}_segmentation_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK NXF_VER="22.04.5" LSB_DEFAULTGROUP='team283' nextflow -trace nextflow.executor run /lustre/scratch126/cellgen/team283/tl10/workflow-segmentation/main.nf \
	-params-file $1 \
	-profile lsf \
	-entry run_nuc_seg \
	-resume
