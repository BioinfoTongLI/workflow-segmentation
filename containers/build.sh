#! /bin/sh
#
# build.sh
# Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the MIT license.
#

docker build -t segmenation:sam -f Dockerfile.sam .
