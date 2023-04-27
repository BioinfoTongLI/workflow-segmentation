#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#


NXF_OPTS='-Dleveldb.mmap=false' \
    NXF_WORK='/tmp/work/' \
    nextflow run /lustre/scratch126/cellgen/team283/tl10/workflow-segmentation/main.nf \
	-profile local \
	-params-file $1 \
	-with-report \
	-entry sam_demo \
	-resume