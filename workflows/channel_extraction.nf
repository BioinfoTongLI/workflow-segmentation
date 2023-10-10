#!/usr/bin/env/ nextflow
// Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.csv = "[path-to-template.csv]"
params.target_ch_indexes = "[1,2,3,4]" //"[4,0]"
params.out_dir = "/nfs/team283_imaging/SM_BRA/playground_Tong/Suzzana_Jimmy_hiplex/omero_roi_segmentation/"

params.docker_container = "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-segmentation:latest"
params.sif_container = "/lustre/scratch126/cellgen/team283/imaging_sifs/workflow-segmentation.sif"
params.ch_index = 0 // can only use 0 for now

include { BIOINFOTONGLI_BIOFORMATS2RAW } from '../modules/sanger/bioinfotongli/bioformats2raw/main' addParams(
    enable_conda:false,
    publish:false,
    store:true,
    out_dir:params.out_dir
)


process export_chs_from_zarr {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    storeDir params.out_dir + "/extracted_chs"

    input:
    tuple val(meta), path(zarr_in)
    val(target_ch_indexes)

    output:
    tuple val(stem), path("${stem}_target_chs.tif"), emit: tif

    script:
    stem=meta["stem"]
    """
    zarr_handler.py to_tiff --stem "$stem" --zarr_in ${zarr_in}/0 --target_ch_indexes ${target_ch_indexes}
    """
}


input_files = Channel.fromPath(params.csv)
    .splitCsv(header:true)
    .multiMap{it ->
        images: [['id':file(it.filepath).baseName], file(it.filepath)]
    }

workflow {
    extract_tif(input_files.images)
}

workflow extract_tif {
    take: img
    main:
        BIOINFOTONGLI_BIOFORMATS2RAW(img)
        export_chs_from_zarr(BIOINFOTONGLI_BIOFORMATS2RAW.out.zarr, params.target_ch_indexes)
    emit: export_chs_from_zarr.out.tif
}
