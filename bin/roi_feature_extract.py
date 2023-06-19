#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
import fire
from aicsimageio import AICSImage
from tqdm import tqdm
import pandas as pd
from pathlib import Path
try:
    from cucim.skimage.measure import regionprops_table
    from cucim.skimage.morphology import remove_small_objects
    import cupy as xp
    print("using cucim")
except:
    from skimage.measure import regionprops_table
    import numpy as xp
    print("using skimage")
import trackpy as tp

import tifffile as tf
import numpy as np
import sys
import os


# segment cells with canonical cellpose API
def extract(stem:str, label_p:str, raw_p:str, s:str, min_object_size:int=2000):
    lab = AICSImage(label_p)
    print(lab)
    raw_stack = AICSImage(raw_p)
    raw_stack.set_scene(s)
    print(raw_stack.dims)
    peaks_in_all_ts = []

    raw_crop_dir = f"{stem}_serie_{s}_raw_crops"
    os.mkdir(raw_crop_dir)
    for t, scene in enumerate(lab.scenes): # this is actually Timepoints
        lab.set_scene(scene)
        lab_stack = lab.get_image_dask_data("YXZ", S=0) # YXZ is skimage compatible order
        lab_stack = remove_small_objects(xp.array(lab_stack),
                                         min_size=min_object_size) # remove small objects
        current_raw_stack = \
                raw_stack.get_image_dask_data("YXZC", T=t).squeeze() # YXZC is skimage compatible order
        try:
            props_table = regionprops_table(xp.array(lab_stack),
                                            intensity_image=xp.array(current_raw_stack),
                                            properties=['label', 'centroid', 'bbox',
                                                        'axis_major_length', 'axis_minor_length', 'area',
                                                        'mean_intensity', 'max_intensity', 'min_intensity',
                                                        'image_intensity',
                                                        # 'area_convex', "inertia_tensor",
                                                        # 'moments',
                                                        # 'coords', 'image_convex'
                                                        # 'solidity', 'equivalent_diameter_area'
                                                        ])
        except:
            print(f"Regionprops error: {stem}_serie_{s}_timepoint {t}")
            continue
        for i, l in enumerate(props_table['label']):
            raw_crop = props_table["image_intensity"][i]
            tf.imwrite(f"{raw_crop_dir}/T{t}_label_{l}.tif", raw_crop.get())

        np_dict = {}
        for k, v in props_table.items():
            if k ==  "image_intensity":
                continue
            try:
                np_dict[k] = v.get()
            except:
                np_dict[k] = v
        d = pd.DataFrame(np_dict)
        d["frame"] = t
        d["Aspect_Ratio"] = d["axis_minor_length"] / d["axis_major_length"]
        peaks_in_all_ts.append(d)
    df = pd.concat(peaks_in_all_ts)
    df.to_csv(f"{stem}_serie_{s}_features.csv", index=False)


if __name__ == "__main__":
    fire.Fire(extract)
