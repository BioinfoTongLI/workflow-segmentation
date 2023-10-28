#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
import fire
import tifffile as tf
from aicsimageio import AICSImage
from cellpose import models, core
from skimage.segmentation import expand_labels
from skimage.morphology import remove_small_objects
from skimage.measure import regionprops_table
import pandas as pd
import yaml
import re


# segment cells with canonical cellpose API
def segment(
    image: str,
    mask: str,
    metadata: str,
    out_name: str,
    s=0,
    C=0,
    distance=0,
    min_object_size=600,
    **cellposeparams,
):
    with open(metadata, "r") as f:
        md_dict = yaml.safe_load(f)
    channel_label_dict = {i: item for i, item in enumerate(md_dict["channel_labels"])}

    img = AICSImage(image)
    img.set_scene(s)

    mask = tf.imread(mask)[:-1, :-1]

    dapi = img.get_image_dask_data("YX", T=0, C=C, Z=0)
    stack = img.get_image_dask_data("YXZ")
    masked_raw = dapi.compute() * mask
    segs, flows, styles, diams = model.eval(masked_raw, **cellposeparams)
    clean_lab = remove_small_objects(segs, min_object_size)
    expanded = expand_labels(clean_lab, distance=distance)
    tf.imwrite(out_name, expanded)
    props = regionprops_table(
        expanded,
        intensity_image=stack.compute(),
        properties=[
            "label",
            "centroid",
            "axis_major_length",
            "axis_minor_length",
            "area",
            "mean_intensity",
            "max_intensity",
            "min_intensity",
            "equivalent_diameter_area",
            "solidity",
            "area_convex",
            # 'image_intensity',
            # "inertia_tensor",
            # 'moments',
            # 'coords', 'image_convex'
        ],
    )
    np_dict = {}
    for k, v in props.items():
        if k == "image_intensity":
            continue
        # assuming variable `line` contains the line of interest
        if re.match(r"^(mean|max|min)_intensity", k):
            tokens = k.split("-")
            name = tokens[0] + "-" + channel_label_dict[int(tokens[-1])]
        else:
            name = k
        try:
            np_dict[name] = v.get()
        except:
            np_dict[name] = v
    pd.DataFrame(np_dict).to_csv(out_name.replace(".tif", ".csv"))


if __name__ == "__main__":
    use_GPU = core.use_gpu()
    print(">>> GPU activated? %d" % use_GPU)
    model = models.Cellpose(gpu=use_GPU, model_type="cyto2")
    fire.Fire(segment)
