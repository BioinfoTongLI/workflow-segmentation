#!/usr/bin/env/ nextflow
// Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.ome_tiff = [
    "/nfs/assembled_images/datasets/OB_ADR/OB_ADR/OB_ADR_01MB_DAJN76.3b_Nucleus_Adra1a_Adra1b_Adra2a_Slc1a3_Meas2_A02_F1T1_max.ome.tif",
    "/nfs/assembled_images/datasets/OB_ADR/OB_ADR/OB_ADR_01MD_MLK632.1b_Nucleus_Adra1a_Adra1b_Adra2a_Slc1a3_Meas3_A04_F1T1_max.ome.tif",
    "/nfs/assembled_images/datasets/OB_ADR/OB_ADR/OB_ADR_01MD_MLK632.1b_Nucleus_Adra1a_Adra1b_Adra2a_Slc1a3_Meas4_A04_F1T1_max.ome.tif",
    "/nfs/assembled_images/datasets/OB_ADR/OB_ADR/OB_ADR_01MC_PMGG13.1a_Nucleus_Adra1a_Adra1b_Adra2a_Slc1a3_Meas2_A03_F1T1_max.ome.tif",
    ]
/*"/nfs/team283_imaging/0HarmonyStitched/LY_BRC/LY_BRC_LY10085__2022-04-13T15_42_53-Measurement 2b_max/LY_BRC_LY10085_PD50756a_Nucleus_RPLP0_p_pan-CK_p_PTPRC_Meas2b_A03_F1T1_max.ome.tif"*/
params.object_diameter = 42
params.chs_for_cell_Seg = "[4,0]"
/*params.out_dir = "/nfs/team283_imaging/OB_ADR/playground_Tong/20220428/"*/
params.out_dir = "/nfs/team283_imaging/OB_ADR/playground_Tong/20220503/"
params.tilesize = 13000

params.container = "/lustre/scratch117/cellgen/team283/tl10/sifs/workflow_segmentation.sif"
params.cyto_pixel_classifier = "/nfs/team283_imaging/OB_ADR/playground_Tong/classifiers/cyto_detection.ilp"
params.tissue_pixel_classifier = "/nfs/team283_imaging/OB_ADR/playground_Tong/classifiers/tissue_finder.ilp"
params.max_fork = 3
/*params.expand_in_pixel = 10*/
params.expand_in_pixel = 300


process slice {

    cache "lenient"
    conda projectDir + "/slice_stitch_conda.yaml"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    path(tif)
    val(tilesize)

    output:
    tuple val(stem), path("${stem}_raw_splits"), emit: tiles
    tuple val(stem), path("${stem}_raw_splits/slicer_info.json"), emit: info

    script:
    stem = tif.baseName
    """
    slicer_runner.py -i ${tif} -o "${stem}_raw_splits" --selected_channels 0 -s ${tilesize}
    """
}


process cellpose_cell_segmentation {
    echo true

    cache "lenient"
    container params.container
    containerOptions "--nv"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    queue "gpu-normal"
    /*queue "gpu-basement"*/
    clusterOptions = "-gpu 'num=1:gmem=2000'"
    cpus 4
    memory 64.GB

    maxForks params.max_fork

    input:
    tuple val(stem), path(tiles)

    output:
    tuple val(stem), path("${stem}_label_splits"), emit: labels

    script:
    """
    python -m cellpose --dir ./${tiles} --use_gpu --diameter ${params.object_diameter} --flow_threshold 0 --chan 0 --pretrained_model cyto2 --save_tif --no_npy
    mkdir ${stem}_label_splits
    mv ${tiles}/*cp_masks.tif ${stem}_label_splits
    """
}


process stitch {

    cache "lenient"
    conda projectDir + "/slice_stitch_conda.yaml"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    tuple val(stem), path(tiles), path(slicer_json)

    output:
    tuple val(stem), path("${stem}_seg.tif")

    script:
    """
    mv slicer_info.json  ${tiles}
    stitcher_runner.py -i ./${tiles} -o ./ --no_cell
    mv z001_mask.ome.tiff ${stem}_seg.tif
    """
}


process bf2raw {

    conda "-c ome bioformats2raw"
    /*publishDir params.out_dir, mode: "copy"//, pattern: "*${params.nuc_ch}*"*/
    storeDir params.out_dir

    input:
    path(img)

    output:
    tuple val(stem), path("${stem}.zarr"), emit: img_zarr

    script:
    stem = img.baseName
    """
    bioformats2raw ${img} ${stem}.zarr --resolutions 8  --no-hcs
    """
}


process export_chs_from_zarr {
    echo true
    cache "lenient"

    container params.container
    containerOptions "--nv"
    /*conda projectDir + "/conda.yaml"*/
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    tuple val(stem), path(zarr_in)
    val(target_ch_indexes)

    output:
    tuple val(stem), path("${stem}_target_chs.tif")

    script:
    """
    zarr_handler.py to_tiff --stem $stem --zarr_in ${zarr_in}/0 --target_ch_indexes ${target_ch_indexes}
    """
}


process dapi_assisted_segmentation_improvement {
    echo true
    cache "lenient"

    container params.container
    containerOptions "--nv"
    /*conda projectDir + "/conda.yaml"*/
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

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
    echo true
    cache "lenient"

    container params.container
    containerOptions "--nv"
    /*conda projectDir + "/conda.yaml"*/
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    maxForks 1

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
    echo true
    container "eu.gcr.io/imaging-gpu-eval/ilastik:latest"
    publishDir params.out_dir, mode:"copy"

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


process ilastik_pixel_classification {
    echo true
    container "eu.gcr.io/imaging-gpu-eval/ilastik:latest"
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

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
    bash /ilastik-1.3.3-Linux/run_ilastik.sh --headless \
        --project=${classifier} \
        --readonly \
        --export_source="simple segmentation" \
        --output_format="tif" \
        --raw_data=${raw_img}

    mkdir ${stem}_${categorie}
    mv *Simple\\ Segmentation.tif ${stem}_${categorie}
    """
}


process find_tissue_border {
    echo true
    cache "lenient"

    container params.container
    containerOptions "--nv"
    /*conda projectDir + "/conda.yaml"*/
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
    echo true
    cache "lenient"

    container params.container
    containerOptions "--nv"
    /*conda projectDir + "/conda.yaml"*/
    /*publishDir params.out_dir, mode:"copy"*/
    storeDir params.out_dir

    input:
    tuple val(stem), path(nuc_label)
    val(distance)

    output:
    tuple val(stem), path("${stem}_label_expanded.tif")

    script:
    """
    expand_labels.py --stem $stem --label ${nuc_label} --distance ${distance}
    """
}


workflow {
    extract_tif(channel.fromPath(params.ome_tiff))
    ilastik_pixel_classification(extract_tif.out, params.cyto_pixel_classifier, "cyto")
    nuc_seg_only(channel.fromPath(params.ome_tiff))
    dapi_assisted_segmentation_improvement(ilastik_pixel_classification.out.join(nuc_seg_only.out))
    Get_complementary_nuc_labels(dapi_assisted_segmentation_improvement.out.join(nuc_seg_only.out))
}

workflow extract_tif {
    take: img
    main:
        bf2raw(img)
        export_chs_from_zarr(bf2raw.out, params.chs_for_cell_Seg)
    emit: export_chs_from_zarr.out
}


workflow nuc_seg_only {
    take: img
    main:
        slice(img, params.tilesize)
        cellpose_cell_segmentation(slice.out.tiles)
        stitch(cellpose_cell_segmentation.out.labels.join(slice.out.info))
        expand_label_image(stitch.out, params.expand_in_pixel)
    emit: expand_label_image.out
}

workflow run_nuc_seg {
    nuc_seg_only(channel.fromPath(params.ome_tiff))
}

workflow run_tissue_seg {
    extract_tif(channel.fromPath(params.ome_tiff))
    ilastik_pixel_classification(extract_tif.out, params.tissue_pixel_classifier, "tissue")
    find_tissue_border(ilastik_pixel_classification.out)
}
