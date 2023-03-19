#!/usr/bin/env/ nextflow
// Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.csv = "[path-to-template.csv]"
params.object_diameter = [70]
params.target_ch_indexes = "[1,2,3,4]" //"[4,0]"
params.out_dir = "[path-to-output-dir]"
params.tilesize = 13000 // for tiled cell segmentation
params.target_cellpose_ch_ind = 1 //target channel to perform cell segmentation on

params.cyto_pixel_classifier = "[path-to-ilastik-cytoplasm-classifier]"

params.tissue_pixel_classifier = "[path-to-ilastik-tissue-classifier]"
params.expand_in_pixel = [10]

params.docker_container = "gitlab-registry.internal.sanger.ac.uk/tl10/workflow-segmentation:latest"
params.sif_container = "/lustre/scratch126/cellgen/team283/imaging_sifs/large_cellseg.sif"
params.max_fork = 3
params.ch_index = 0 // can only use 0 for now

include { BIOINFOTONGLI_BIOFORMATS2RAW as bf2raw } from '/lustre/scratch126/cellgen/team283/tl10/modules/modules/bioinfotongli/bioformats2raw/main.nf' addParams(
    enable_conda:false,
    publish:false,
    store:true,
    out_dir:params.out_dir
)


process slice {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/slices"

    cpus 1
    memory 130.GB
    queue "imaging"

    input:
    path(tif)
    val(tilesize)
    val(ch_index)

    output:
    tuple val(stem), path("${stem}_raw_splits"), emit: tiles
    tuple val(stem), path("${stem}_raw_splits/slicer_info.json"), emit: info

    script:
    stem = tif.baseName
    """
    slicer_runner.py -i ${tif} -o "${stem}_raw_splits" --selected_channels ${ch_index} -s ${tilesize}
    """
}


process cellpose_cell_segmentation_batch {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/sparse_segs"

    queue "gpu-normal"
    /*queue "gpu-basement"*/
    clusterOptions = {" -gpu \"mode=shared:j_exclusive=no:gmem=12000:num=1\""}
    cpus 10
    memory 64.GB

    maxForks params.max_fork

    input:
    tuple val(stem), path(tiles)
    each cell_size

    output:
    tuple val(stem), path("${stem}_label_splits_cell_size_${cell_size}"), val(cell_size), emit: labels

    script:
    """
    python -m cellpose --dir ./${tiles} --use_gpu --diameter ${cell_size} --flow_threshold 0 --chan 0 --pretrained_model cyto2 --save_tif --no_npy
    mkdir "${stem}_label_splits_cell_size_${cell_size}"
    mv ${tiles}/*cp_masks.tif "${stem}_label_splits_cell_size_${cell_size}"
    """
}


process cellpose_cell_segmentation {
    echo true

    cache "lenient"
    container params.container
    containerOptions params.contOptions
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/cellpose_seg"

    queue "gpu-normal"
    /*queue "gpu-basement"*/
    clusterOptions = "-gpu 'num=1:gmem=2000'"
    cpus 4
    memory 64.GB

    maxForks params.max_fork

    input:
    path(img)
    each cell_size
    val target_ch_ind

    output:
    tuple val(stem), path("${stem}_cp_masks.tif"), val(cell_size), emit: labels

    script:
    stem = img.baseName
    """
    python -m cellpose --dir ./ --use_gpu --diameter ${cell_size} --flow_threshold 0 --chan ${target_ch_ind} --pretrained_model cyto2 --save_tif --no_npy
    """
}


process stitch {
    debug true

    cpus 1
    /*label 'huge_mem'*/
    memory { 330.GB * task.attempt }
    /*time { 1.hour * task.attempt }*/

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    /*containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"*/
    /*publishDir params.out_dir + "/stitched_seg", mode:"copy"*/
    storeDir params.out_dir + "/stitched_seg"

    input:
    tuple val(stem), path(tiles), val(cell_size), path(slicer_json)

    output:
    tuple val(stem), path(out_tif), val(cell_size)

    script:
    out_tif = "${stem}_seg_${cell_size}.tif"
    """
    cp slicer_info.json  ${tiles}
    stitcher_runner.py -i ${tiles} -o ./ --no_cell
    mv z001_mask.ome.tiff "${out_tif}"
    """
}


process export_chs_from_zarr {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir + "/extracted_chs", mode:"copy"*/
    storeDir params.out_dir + "/extracted_chs"
    /*errorStrategy "ignore"*/

    input:
    tuple val(meta), path(zarr_in)
    val(target_ch_indexes)

    output:
    tuple val(stem), path("${stem}_target_chs.tif"), emit: tif
    /*tuple val(stem), path("${stem}_target_chs.npy"), emit: tif*/

    script:
    stem=meta["stem"]
    """
    zarr_handler.py to_tiff --stem "$stem" --zarr_in ${zarr_in}/0 --target_ch_indexes ${target_ch_indexes}
    """
}


process dapi_assisted_segmentation_improvement {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/nuclei_guided_segs"

    maxForks 2

    input:
    tuple val(stem), path(cyto_seg), path(nuc_seg)

    output:
    tuple val(stem), path("${stem}_improved_seg.tif")

    script:
    """
    nuclei_assisted_seg.py --stem $stem --cyto_seg $cyto_seg/*Simple\\ Segmentation.tif --nuc_seg ${nuc_seg}
    """
}


process Get_complementary_nuc_labels {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/nuc_seg_complementary"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    memory 150.GB

    input:
    tuple val(stem), path(cyto_seg), path(nuc_seg)

    output:
    tuple val(stem), path("${stem}_improved_seg_complementary.tif")
    tuple val(stem), path("${stem}_improved_seg_complementary_non_cell.tif")

    script:
    """
    exclude_labels.py --stem $stem --label $cyto_seg --nuc ${nuc_seg}
    """
}


process ilastik_cell_filtering {
    debug true

    /*container "eu.gcr.io/imaging-gpu-eval/ilastik:latest"*/
    container "gitlab-registry.internal.sanger.ac.uk/tl10/img-ilastik:latest"
    publishDir params.out_dir + "/classification", mode:"copy"

    /*machineType { ['n2-highmem-8','n2-highmem-16','n2-highmem-32','n2-highmem-64'][task.attempt-1] }*/
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 4
    memory 60.GB

    input:
    tuple val(stem), file(raw_img), file(mask)
    path(project_file)

    output:
    tuple val(stem), file("$stem*object_features_table.csv")
    tuple val(stem), file("$stem*Object Predictions.tif")

    script:
    """
    #LAZYFLOW_THREADS=10 LAZYFLOW_TOTAL_RAM_MB=60000
    bash /opt/ilastik/run_ilastik.sh --headless \
        --project=${project_file} \
        --readonly=yes \
        --table_filename="./${stem}_object_features.csv" \
        --export_source="Blockwise Object Predictions" \
        --output_format="tif" \
        --output_filename_format="./${stem}_{result_type}.tif" \
        --raw_data=${raw_img} \
        --segmentation_image=${mask}
    """
}


process ilastik_pixel_classification {
    debug true
    container "eu.gcr.io/imaging-gpu-eval/ilastik:latest"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir + "/pixel_prediction"

    machineType { ['n2-highmem-8','n2-highmem-16','n2-highmem-32','n2-highmem-64'][task.attempt-1] }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 4
    maxForks 2

    input:
    tuple val(stem), path(raw_img)
    path(classifier)
    val(categorie)

    output:
    /*tuple val(stem), file("$stem*table.csv")*/
    tuple val(stem), path("${stem}_${categorie}")

    script:
    """
    #LAZYFLOW_THREADS=10 LAZYFLOW_TOTAL_RAM_MB=60000
    bash /opt/ilastik/run_ilastik.sh --headless \
        --project=${classifier} \
        --readonly 1 \
        --export_source="simple segmentation" \
        --output_format="tif" \
        --raw_data=${raw_img}

    mkdir ${stem}_${categorie}
    mv *Simple\\ Segmentation.tif ${stem}_${categorie}
    """
}


process find_tissue_border {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    publishDir params.out_dir, mode:"copy"
    /*storeDir params.out_dir*/

    maxForks 2

    input:
    tuple val(stem), path(tissue_seg)

    output:
    tuple val(stem), path("${stem}_tissue_mask.tif"), path("${stem}_tissue_contour.wkt")

    script:
    """
    find_tissue_outline.py --stem $stem --tissue_seg $tissue_seg/*Simple\\ Segmentation.tif
    """
}


process expand_label_image {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.sif_container:
        params.docker_container}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir + "/expanded_label", mode:"copy"*/
    storeDir params.out_dir + "/expanded_label"

    input:
    tuple val(stem), path(nuc_label), val(cell_size), val(distance)


    output:
    tuple val(stem), path("${stem}_${cell_size}_label_expanded_by_${distance}.tif"), emit: tif

    script:
    """
    expand_labels.py --stem "$stem" --label ${nuc_label} --distance ${distance}
    mv "${stem}_label_expanded.tif" "${stem}_${cell_size}_label_expanded_by_${distance}.tif"
    """
}


input_files = Channel.fromPath(params.csv)
    .splitCsv(header:true)
    .multiMap{it ->
        images_for_bf2raw: [[stem:file(it.filepath).baseName], file(it.filepath)]
        images: file(it.filepath)
    }

workflow {
    extract_tif(input_files.images)
    /*ilastik_pixel_classification(extract_tif.out, params.cyto_pixel_classifier, "cyto")*/
    nuc_seg_only(input_files.images)
    dapi_assisted_segmentation_improvement(ilastik_pixel_classification.out.join(nuc_seg_only.out))
    /*Get_complementary_nuc_labels(dapi_assisted_segmentation_improvement.out.join(nuc_seg_only.out))*/
}

workflow extract_tif {
    take: img
    main:
        bf2raw(img)
        export_chs_from_zarr(bf2raw.out.zarr, params.target_ch_indexes)
    emit: export_chs_from_zarr.out.tif
}

workflow nuc_seg_only {
    take: img
    main:
        slice(img, params.tilesize, params.ch_index)
        cellpose_cell_segmentation_batch(slice.out.tiles, params.object_diameter)
        stitch(cellpose_cell_segmentation_batch.out.labels.join(slice.out.info))
        expand_label_image(stitch.out.combine(params.expand_in_pixel))
    emit: expand_label_image.out.tif
}

workflow run_nuc_seg {
    nuc_seg_only(input_files.images)
}

workflow run_tissue_seg {
    extract_tif(input_files.images_for_bf2raw)
    /*ilastik_pixel_classification(extract_tif.out, params.tissue_pixel_classifier, "tissue")*/
    /*find_tissue_border(ilastik_pixel_classification.out)*/
}

workflow run_cell_classification {
    nuc_seg_only(input_files.images)
    extract_tif(input_files.images_for_bf2raw)
    ilastik_cell_filtering(extract_tif.out.join(nuc_seg_only.out), params.cyto_pixel_classifier)
}
