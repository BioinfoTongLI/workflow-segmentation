#!/usr/bin/env/ nextflow
// Copyright © 2023 Tong LI <tongli.bioinfo@proton.me>

nextflow.enable.dsl=2

params.to_seg = [
    [0, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellA1_ChannelSD Red FW,SD Green FW_Seq0002.nd2", "0,0", 0],
    [1, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellA1_ChannelSD Red FW,SD Green FW_Seq0002.nd2", "0,0", 1],
    [2, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellA1_ChannelSD Red FW,SD Green FW_Seq0002.nd2", "0,0", 2],
    [3, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellB1_ChannelSD Red FW,SD Green FW_Seq0003.nd2", "0,0", 0],
    [4, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellB1_ChannelSD Red FW,SD Green FW_Seq0003.nd2", "0,0", 1],
    [5, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellB1_ChannelSD Red FW,SD Green FW_Seq0003.nd2", "0,0", 3],
    [6, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellB1_ChannelSD Red FW,SD Green FW_Seq0003.nd2", "0,0", 4],
    [7, "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/T265C NTG1/WellB1_ChannelSD Red FW,SD Green FW_Seq0003.nd2", "0,0", 5],
]
params.out_dir = "/lustre/scratch126/cellgen/team283/tl10/Freddy_3D_seg/segmentations/"
params.sif_container = "/lustre/scratch126/cellgen/team283/imaging_sifs/workflow-segmentation.sif"

process cellpose_3d_seg {
    debug true

    label 'gpu_normal'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        "workflow-segmentation:latest"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv -B /lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models:/lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models':'--gpus all -v /lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models:/lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models'}"
    /*publishDir params.out_dir, mode: 'copy'*/
    storeDir params.out_dir + "segmentations/"

    input:
    tuple val(id), path(img), val(chs), val(serie)

    output:
    tuple val(id), val(meta), path("${meta['stem']}_serie_${serie}_chs_${chs}_3d_label_array.tif"), emit: segmentation

    script:
    meta = [:]
    meta['id'] = id
    meta['stem'] = img.baseName
    meta['serie'] = serie
    meta["channel"] = chs
    """
    export CELLPOSE_LOCAL_MODELS_PATH=/lustre/scratch126/cellgen/team283/NXF_WORK/cellpose_models
    cellpose_3d_seg.py --stem "${meta['stem']}" --img_p ${img} --chs ${chs} --s ${serie}
    """
}


process feature_extraction {
    debug true

    label 'gpu_normal'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        "workflow-segmentation:latest"}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode: 'copy'*/
    storeDir params.out_dir + "features/"

    input:
    tuple val(id), path(raw), val(meta), path(img)

    output:
    tuple val(meta), path("${stem}_serie_${serie}_features.csv"), emit: centroids
    tuple val(meta), path("${stem}_serie_${serie}_raw_crops"), emit: raw_crops

    script:
    stem = meta['stem']
    serie = meta['serie']
    """
    roi_feature_extract.py --stem "${stem}" --label_p ${img} --raw_p ${raw} --s ${serie}
    """
}


process track {
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
    val(thre)

    output:
    path("${stem}_serie_${serie}_*.png"), emit: plots
    path("${stem}_serie_${serie}_*.csv"), emit: tracks

    script:
    stem = meta['stem']
    serie = meta['serie']
    """
    track_with_trackpy.py --stem "${stem}_serie_${serie}" --centroids_p ${centroids} --cell_volumn_threshold ${thre}
    """
}


workflow {
    cellpose_3d_seg(channel.from(params.to_seg))
    Channel.from(params.to_seg)
        .map { [it[0], it[1]] }
        .join(cellpose_3d_seg.out).view()
        /*.map { [["id":file(it[0]).baseName], it[0]}*/
    feature_extraction(Channel.from(params.to_seg)
        .map { [it[0], it[1]]}
        .join(cellpose_3d_seg.out)
    )
    track(feature_extraction.out.centroids, 15000)
}
