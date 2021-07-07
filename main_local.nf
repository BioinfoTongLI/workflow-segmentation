#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.img_path = "/nfs/team283_imaging/"
params.dapi_ch = 0
params.cyto_ch = 1
params.object_diameter = 70
params.flow_threshold = 0.0
params.model_type = "cyto"
params.magnification = 40
params.expand = 70
params.max_fork = 2
params.out_dir = ""


if (params.cyto_ch != ''){
    params.cellpose_ch2 = '--chan2 1'
} else{
    params.cellpose_ch2 = ''
}


process extract_channels_from_input_image {
    cache "lenient"
    container "segmentation_helper:latest"
    storeDir params.out_dir

    cpus { 4 * task.attempt }
    memory { 25.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path img

    output:
    tuple val(stem), path("${stem}_*DAPI.tif")

    script:
    stem = img.baseName
    """
    helper extract_channels ${params.dapi_ch} ${params.cyto_ch} --img_in_path=${img} --stem=${stem}
    """
}


process cellpose_cell_segmentation {
    cache "lenient"
    echo true
    /*container "gitlab-registry.internal.sanger.ac.uk/tl10/img-cellpose:latest"*/
    container "eu.gcr.io/imaging-gpu-eval/cellpose:latest"
    containerOptions "--gpus all"
    /*label "cellpose"*/
    storeDir params.out_dir

    maxForks params.max_fork

    input:
    tuple val(stem), file(ch_img)

    output:
    tuple val(stem), file("${stem}*cp_masks.tif")

    script:
    """
    python -m cellpose --dir ./ --use_gpu --diameter ${params.object_diameter} --flow_threshold ${params.flow_threshold} --chan 0 $params.cellpose_ch2 --pretrained_model ${params.model_type} --save_tif --no_npy
    """
}


process expand_labels {
    echo true
    container "eu.gcr.io/imaging-gpu-eval/expand_labels:latest"
    storeDir params.out_dir

    machineType { ['n2-highmem-8','n2-highmem-16','n2-highmem-32','n2-highmem-64'][task.attempt-1] }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 4
    maxForks 1

    input:
    tuple val(stem), file(label)

    output:
    tuple val(stem), file("${stem}*_label_expanded.tif")

    script:
    """
    expand_labels -label ${label} -prefix ${stem} -distance ${params.expand}
    """
}

process ilastik_cell_filtering {
    echo true
    container "eu.gcr.io/imaging-gpu-eval/ilastik:latest"
    storeDir params.out_dir + "/" + params.magnification + "x_cell_mask"

    machineType { ['n2-highmem-8','n2-highmem-16','n2-highmem-32','n2-highmem-64'][task.attempt-1] }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 4
    maxForks 1

    input:
    tuple val(stem), file(raw_img), file(mask)

    output:
    tuple val(stem), file("$stem*table.csv")
    tuple val(stem), file("$stem*Object Predictions.tif")

    script:
    """
    #LAZYFLOW_THREADS=10 LAZYFLOW_TOTAL_RAM_MB=60000
    bash /ilastik-1.3.3-Linux/run_ilastik.sh --headless \
        --project=/model/project.ilp \
        --readonly \
        --table_filename="./${stem}_object_features.csv" \
        --export_source="Blockwise Object Predictions" \
        --output_format="tif" \
        --raw_data=${raw_img} \
        --segmentation_image=${mask}
    """
}

process filter_labels {
    echo true
    container "eu.gcr.io/imaging-gpu-eval/filter_labels:latest"
    storeDir params.out_dir

    cpus { 4 * task.attempt }
    memory { 25.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(stem), file(label), file(mask), file(roi)

    output:
    tuple val(stem), file("$stem*filtered_labels.tif")
    tuple val(stem), file("$stem*filtered_rois.tsv")

    script:
    """
    filter_labels -label ${label} -prefix ${stem} -mask ${mask} -rois ${roi}
    """
}

workflow {
    extract_channels_from_input_image(channel.fromPath(params.img_path))
    cellpose_cell_segmentation(extract_channels_from_input_image.out)
    expand_labels(cellpose_cell_segmentation.out)
    ilastik_cell_filtering((extract_channels_from_input_image.out.join(expand_labels.out)))
    filter_labels(expand_labels.out.join(ilastik_cell_filtering.out[1]).join(ilastik_cell_filtering.out[0]))
}
