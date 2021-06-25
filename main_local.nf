#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.img_path = ""
params.dapi_ch = 0
params.cyto_ch = 1
params.object_diameter = 70
params.flow_threshold = 0.0
params.model_type = "cyto"
params.magnification = 40
params.expand = 70
params.max_fork = 4
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
    tuple val(stem), path("${stem}_*.tif")

    script:
    stem = img.baseName
    """
    helper extract_channels ${params.dapi_ch} ${params.cyto_ch} --img_in_path=${img} --stem=${stem}
    """
}


process cellpose_cell_segmentation {
    cache "lenient"
    echo true
    container "gitlab-registry.internal.sanger.ac.uk/tl10/img-cellpose:latest"
    containerOptions "--gpus all"
    /*label "cellpose"*/
    storeDir params.out_dir

    maxForks params.max_fork

    input:
    tuple val(stem), file(ch_img)

    output:
    tuple val(stem), file("${stem}*seg.npy")

    script:
    """
    ls
    python -m cellpose --dir ./ --use_gpu --diameter ${params.object_diameter} --flow_threshold ${params.flow_threshold} --chan 0 $params.cellpose_ch2 --pretrained_model ${params.model_type}
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
    tuple val(stem), file(label) from labels

    output:
    tuple val(stem), file("${stem}*_label_expanded.tif") into expanded_nuc_labels, labels_for_filtering

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
    tuple val(stem), file(raw_img), file(mask) from target_chs_for_filtering.join(expanded_nuc_labels)

    output:
    tuple val(stem), file("$stem*table.csv") into cell_roi_descriptions
    tuple val(stem), file("$stem*Object Predictions.tif") into selective_masks, selective_masks_for_intensity_quant

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
    tuple val(stem), file(label), file(mask), file(roi) from labels_for_filtering.join(selective_masks).join(cell_roi_descriptions)

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
}
