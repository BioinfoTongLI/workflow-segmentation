#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
import fire
from aicsimageio import AICSImage
import pandas as pd

try:
    from cucim.skimage.measure import regionprops_table
    import cupy as xp

    print("using cucim")
except:
    from skimage.measure import regionprops_table
    import numpy as xp

    print("using skimage")

import os


def quantify(lab, raw_stack, t, s, raw_crop_dir, dims=None):
    """
    Common function to quantify objects in a 3D or 2D stack.

    :param lab: Labelled image stack.
    :type lab: dask.array
    :param raw_stack: Raw image stack.
    :type raw_stack: dask.array
    :param t: Timepoint to process.
    :type t: int
    :param s: Series to process.
    :type s: int
    :param raw_crop_dir: Directory to save raw crop images.
    :type raw_crop_dir: str
    :param dims: Z-dimension to process. If None, process entire 3D stack.
    :type dims: int, optional
    :return: Table of object properties.
    :rtype: pandas.DataFrame
    """
    if dims is None:
        # Process entire 3D stack
        lab_stack = lab.get_image_dask_data("YXZ")  # YXZ is skimage compatible order
        current_raw_stack = raw_stack.get_image_dask_data(
            "YXZC", T=t
        ).squeeze()  # YXZC is skimage compatible order
    else:
        # Process single 2D slice
        current_raw_stack = raw_stack.get_image_dask_data("YXC", T=t, Z=dims).squeeze()
        lab_stack = lab.get_image_dask_data("YX", T=t, Z=dims).squeeze()

    print(current_raw_stack.shape, dims)

    try:
        props_table = regionprops_table(
            xp.array(lab_stack),
            intensity_image=xp.array(current_raw_stack),
            properties=[
                "label",
                "centroid",
                "bbox",
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
    except:
        print(f"Regionprops error: serie_{s}_timepoint {t}")
        return None

    np_dict = {}
    for k, v in props_table.items():
        if k == "image_intensity":
            continue
        try:
            np_dict[k] = v.get()
        except:
            np_dict[k] = v
    d = pd.DataFrame(np_dict)
    d["frame"] = t
    d["Aspect_Ratio"] = d["axis_minor_length"] / d["axis_major_length"]

    # if dims is None:
    # for i, l in enumerate(result['label']):
    # raw_crop = result["image_intensity"][i]
    # tf.imwrite(f"{raw_crop_dir}/T{t}_label_{l}.tif", raw_crop.get())

    return d


def quantify_3D(lab, raw_stack, t, s, raw_crop_dir):
    """
    Quantify objects in a 3D stack.

    :param lab: Labelled image stack.
    :type lab: dask.array
    :param raw_stack: Raw image stack.
    :type raw_stack: dask.array
    :param t: Timepoint to process.
    :type t: int
    :param s: Series to process.
    :type s: int
    :param raw_crop_dir: Directory to save raw crop images.
    :type raw_crop_dir: str
    :return: Table of object properties.
    :rtype: pandas.DataFrame
    """
    return quantify(lab, raw_stack, t, s, raw_crop_dir)


def quantify_serial_2D(lab, raw_stack, t, s, raw_crop_dir):
    """
    Quantify objects in a serial 2D stack.

    :param lab: Labelled image stack.
    :type lab: dask.array
    :param raw_stack: Raw image stack.
    :type raw_stack: dask.array
    :param t: Timepoint to process.
    :type t: int
    :param s: Series to process.
    :type s: int
    :param raw_crop_dir: Directory to save raw crop images.
    :type raw_crop_dir: str
    :return: Table of object properties.
    :rtype: pandas.DataFrame
    """
    result = []
    for z in range(raw_stack.dims.Z):
        quant = quantify(lab, raw_stack, t, s, raw_crop_dir, dims=z)
        quant["z"] = z
        if quant is not None:
            result.append(quant)
    if len(result) == 0:
        return None
    return pd.concat(result)


# segment cells with canonical cellpose API
def extract(
    stem: str,
    label_p: str,
    raw_p: str,
    s: str,
    serial: bool = False,
):
    lab = AICSImage(label_p)
    print(lab.dims)
    raw_stack = AICSImage(raw_p)
    raw_stack.set_scene(s)
    print(raw_stack.dims)
    peaks_in_all_ts = []

    raw_crop_dir = f"{stem}_serie_{s}_raw_crops"
    # os.mkdir(raw_crop_dir)

    for t, scene in enumerate(lab.scenes):  # this is actually Timepoints
        lab.set_scene(scene)
        if serial:
            d = quantify_serial_2D(lab, raw_stack, t, s, raw_crop_dir)
        else:
            d = quantify_3D(lab, raw_stack, t, s, raw_crop_dir)
        peaks_in_all_ts.append(d)
    df = pd.concat(peaks_in_all_ts)
    df.to_csv(f"{stem}_serie_{s}_features.csv", index=False)


def extract_2d(
    stem: str,
    label_p: str,
    s: str,
):
    lab = AICSImage(label_p)
    print(lab.dims)

    props_of_scene = {}
    for i, scene in enumerate(lab.scenes):
        lab.set_scene(scene)
        lab_stack = lab.get_image_dask_data("YX")
        props_table = regionprops_table(
            xp.array(lab_stack),
            properties=[
                "label",
                "centroid",
                "bbox",
                "axis_major_length",
                "axis_minor_length",
                "area",
                "equivalent_diameter_area",
                "solidity",
                "area_convex",
                # "mean_intensity",
                # "max_intensity",
                # "min_intensity",
                # "inertia_tensor",
                # 'moments',
                # 'coords', 'image_convex'
            ],
        )
        np_dict = {}
        for k, v in props_table.items():
            if k == "image_intensity":
                continue
            try:
                np_dict[k] = v.get()
            except:
                np_dict[k] = v
        props_of_scene[scene] = pd.DataFrame(np_dict)
    pd.concat(props_of_scene).to_csv(f"{stem}_serie_{s}_features.csv", index=False)


if __name__ == "__main__":
    fire.Fire()
