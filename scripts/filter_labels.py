#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import argparse
import tifffile as tf
import numpy as np
import pandas as pd


def main(args):
    label_img = tf.imread(args.label)
    mask_img = tf.imread(args.mask)
    filtered_label = label_img * (mask_img == 2)
    rois = pd.read_csv(args.rois).set_index("object_id")
    filtered_cell_ids = np.unique(filtered_label)
    print(rois.loc[filtered_cell_ids].shape)
    rois.to_csv("%s_filtered_rois.tsv" %args.prefix, sep="\t")
    tf.imsave("%s_filtered_labels.tif" %args.prefix, filtered_label, dtype=np.uint32)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-label", type=str, required=True)
    parser.add_argument("-prefix", type=str, required=True)
    parser.add_argument("-mask", type=str, required=True)
    parser.add_argument("-rois", type=str, required=True)

    args = parser.parse_args()

    main(args)
