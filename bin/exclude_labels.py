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
try:
    import cupy as xp
    print("Using Cupy")
except:
    import numpy as xp
from dask_image.imread import imread
import tifffile as tf


def main(stem, label, nuc):
    lab = xp.asarray(imread(label).squeeze())
    max_seg_id = xp.max(lab)
    nuc = xp.asarray(imread(nuc).squeeze())
    print(lab, nuc)
    tf.imwrite(f"{stem}_improved_seg_complementary_non_cell.tif", (nuc * ~(lab > 0)).get())
    xp._default_memory_pool.free_all_blocks()
    nucs_to_remove = xp.unique(nuc * (lab > 0))
    del lab
    not_in_label_nuclei = nuc.copy()
    del nuc
    for id_to_remove in nucs_to_remove:
        # not_in_label_nuclei = xp.where(not_in_label_nuclei == id_to_remove, 0, not_in_label_nuclei)
        not_in_label_nuclei *= ~(not_in_label_nuclei == id_to_remove)
    # print(not_in_label_nuclei)
    not_in_label_nuclei_mask = not_in_label_nuclei > 0
    not_in_label_nuclei = (not_in_label_nuclei + max_seg_id) * not_in_label_nuclei_mask
    tf.imwrite(f"{stem}_improved_seg_complementary.tif", not_in_label_nuclei.get())


if __name__ == "__main__":
    fire.Fire(main)
