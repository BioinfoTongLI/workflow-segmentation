#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/tmp/work/'

TMP_NF_WORK="$MOUNT_POINT/${USER}_segmentation_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run /scratch/iss_decoding/nf/workflow-segmentation/main.nf \
	-params-file $1 \
	-profile lsf \
	-entry run_nuc_seg \
	-resume
