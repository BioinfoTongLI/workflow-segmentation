#!/usr/bin/env/ nextflow
// Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

include {run_nuc_seg; small_image_cellpose} from './workflows/segmentation'

workflow {
    run_nuc_seg()
}


workflow small {
    small_image_cellpose()
}