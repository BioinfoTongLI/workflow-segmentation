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
from dask_image.imread import imread
# from skimage.measure import regionprops
from cucim.skimage.measure import regionprops_table, regionprops
import tifffile as tf
import numpy as np
import dask.array as da
from scipy import ndimage as ndi
from skimage.segmentation import watershed
import cupy as cp
import pandas as pd


def main(stem, cyto_seg, nuc_seg):
    nuc = imread(nuc_seg).squeeze().compute().astype(np.uint32)
    # print(cyto, nuc)
    # meas = regionprops(nuc)
    meas = regionprops(cp.asarray(nuc))
    print(len(meas))
    # meas_df = pd.DataFrame(meas)
    # print(meas_df)
    all_centroids = []
    for m in meas:
        try:
            all_centroids.append(m.centroid)
        except:
            pass
    all_centroids = np.array(all_centroids).astype(np.uint32)
    print(all_centroids)
    nuc_centroids = cp.zeros_like(nuc, dtype=np.uint32)
    for i, c in enumerate(all_centroids):
        nuc_centroids[c[0], c[1]] = i + 1
    del meas
    del all_centroids
    nuc_centroids_np = nuc_centroids.get()
    del nuc_centroids
    cyto = imread(cyto_seg).squeeze()
    cyto = cyto == 2
    distance = ndi.distance_transform_edt(cyto.compute())
    # coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=nuc_seg)
    # mask = np.zeros(distance.shape, dtype=bool)
    # mask[tuple(coords.T)] = True
    # markers, _ = ndi.label(mask)
    # labels = watershed(-distance, markers, mask=large_masks)
    labels = watershed(-distance, nuc_centroids_np, mask=cyto)
    # large_astros = remove_small_objects(labels, 300)
    tf.imwrite(f"{stem}_improved_seg.tif", labels)



if __name__ == "__main__":
    fire.Fire(main)
