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
import argparse
import numpy as np
from skimage.segmentation import expand_labels
import tifffile as tf


def main(args):
    label = tf.imread(args.label).squeeze().astype(np.uint32)
    # if args.mask:
        # mask = tf.imread(args.mask)
        # masked_label = label * (mask == 2)
    # else:
        # print("Mask doesn't exist, skipped")
    masked_label = label

    expanded = expand_labels(masked_label, distance=args.distance)
    tf.imwrite("%s_label_expanded.tif" % args.prefix, expanded, dtype=np.uint32)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-label", type=str, required=True)
    # parser.add_argument("-mask", type=str, required=True)
    parser.add_argument("-prefix", type=str, required=True)
    parser.add_argument("-distance", type=int, required=True)
    # parser.add_argument("-calibration", type=int, required=True)

    args = parser.parse_args()

    main(args)
