#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Expand label images
"""
import numpy as np
from skimage.segmentation import expand_labels
import tifffile as tf
import fire


def main(stem, label, distance, ilastik_mask=None):
    label = tf.imread(label).squeeze()
    if ilastik_mask:
        mask = tf.imread(ilastik_mask) == 2
        masked_label = label * mask
    else:
        masked_label = label
    expanded = expand_labels(masked_label, distance=distance)
    tf.imwrite(
        f"{stem}_label_expanded.tif",
        expanded,
        metadata={"axes": "YX"},
        # imagej=True,
        dtype=np.uint32,
    )


if __name__ == "__main__":
    fire.Fire(main)
