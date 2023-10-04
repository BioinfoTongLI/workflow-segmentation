#!/usr/bin/env/ nextflow
// Copyright © 2023 Tong LI <tongli.bioinfo@proton.me>

nextflow.enable.dsl=2

params.to_seg = [
    [0, '[img-path]', "[channels-to-seg, _e.g._0,0", '[serie_integer]', 20, 0.0],
]
params.out_dir = null
params.sif_container = "/lustre/scratch126/cellgen/team283/imaging_sifs/workflow-segmentation.sif"


process cellpose_3d_seg {
    debug true

    label 'gpu_normal'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "bioinfotongli/workflow-segmentation:latest":
        "bioinfotongli/workflow-segmentation:latest"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv -B /lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models:/tmp/cellpose_models -B /nfs:/nfs':'--gpus all -v /lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models:/tmp/cellpose_models'}"
    /*publishDir params.out_dir, mode: 'copy'*/
    storeDir params.out_dir + "/segmentations/"

    input:
    tuple val(id), path(img), val(chs), val(serie), val(diameter), val(cellprob_threshold)

    output:
    tuple val(id), val(meta), path("${meta['stem']}_serie_${serie}_chs_${chs}_3d_label_array.tif"), emit: segmentation

    script:
    meta = [:]
    meta['id'] = id
    meta['stem'] = img.baseName
    meta['serie'] = serie
    meta["channel"] = chs
    """
    export CELLPOSE_LOCAL_MODELS_PATH=/tmp/cellpose_models
    cellpose_3d_seg.py --stem "${meta['stem']}" --img_p ${img} --chs ${chs} --s ${serie} --diameter ${diameter} \
        cellprob_threshold ${cellprob_threshold}
    """
}


process feature_extraction {
    debug true

    label 'gpu_normal'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "bioinfotongli/workflow-segmentation:latest":
        "bioinfotongli/workflow-segmentation:latest"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode: 'copy'*/
    storeDir params.out_dir + "/features/"

    input:
    tuple val(id), path(raw), val(meta), path(img)

    output:
    tuple val(meta), path("${stem}_serie_${serie}_features.csv"), emit: centroids
    tuple val(meta), path("${stem}_serie_${serie}_raw_crops"), emit: raw_crops

    script:
    stem = meta['stem']
    serie = meta['serie']
    """
    roi_feature_extract.py \
        --stem "${stem}" \
        --label_p ${img} \
        --raw_p ${raw} \
        --s ${serie}
    """
}


process TRACKPY_TRACKING {
    debug true

    label 'cpu_only'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        "workflow-segmentation:latest"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode: 'copy'*/
    storeDir params.out_dir + "tracks/"

    input:
    tuple val(meta), path(centroids)

    output:
    path("${stem}_serie_${serie}_*.png"), emit: plots
    path("${stem}_serie_${serie}_*.csv"), emit: tracks

    script:
    stem = meta['stem']
    serie = meta['serie']
    def args = task.ext.args ?: ''
    """
    track_with_trackpy.py \
        --stem "${stem}_serie_${serie}" \
        --centroids_p ${centroids} \
        $args
    """
}


workflow Segmentation_3D {
    cellpose_3d_seg(channel.from(params.to_seg))
    Channel.from(params.to_seg)
        .map { [it[0], it[1]] }
        .join(cellpose_3d_seg.out).view()
        /*.map { [["id":file(it[0]).baseName], it[0]}*/
    feature_extraction(Channel.from(params.to_seg)
        .map { [it[0], it[1]]}
        .join(cellpose_3d_seg.out)
    )
    emit: feature_extraction.out.centroids
}

workflow Tracking {
    Segmentation_3D()
    TRACKPY_TRACKING(Segmentation_3D.out)
}
