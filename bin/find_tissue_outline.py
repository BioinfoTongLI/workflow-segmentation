#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import tifffile as tf
from cucim.skimage.measure import regionprops, label
import cupy as cp
import pandas as pd
from scipy.ndimage.morphology import binary_fill_holes
from shapely.geometry import LineString
from shapely import wkt
import cv2
import numpy as np


def main(stem, tissue_seg):
    img = tf.imread(tissue_seg) == 2
    lab = label(cp.asarray(img), background=0)
    largest_area = -1
    tissue_meas = None
    meas = regionprops(lab)
    for m in meas:
        area = m.area
        if area > largest_area:
            largest_area = area
            tissue_meas = m
    tissue_mask = lab == tissue_meas.label
    tissue_mask_filled = binary_fill_holes(tissue_mask.get())
    tf.imwrite(f"{stem}_tissue_mask.tif", tissue_mask_filled)

    # outline = LineString(tissue_meas.coords.get()) //super slow since because of pixel-levle accurate
    tissue_cnt, _ = cv2.findContours(
        tissue_mask_filled.astype(np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE
    )
    assert len(tissue_cnt) == 1
    print(tissue_cnt[0].shape)
    tissue_cnt_shapely = LineString(tissue_cnt[0].squeeze())
    tissue_cnt_wkt = wkt.dumps(tissue_cnt_shapely)
    with open(f"{stem}_tissue_contour.wkt", "w") as fh:
        fh.write(tissue_cnt_wkt)


if __name__ == "__main__":
    fire.Fire(main)
