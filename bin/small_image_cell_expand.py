#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the MIT license.

"""

"""
import fire
import tifffile as tf
import numpy as np
from skimage.segmentation import expand_labels
from cellpose import models, core
use_GPU = core.use_gpu()
print('>>> GPU activated? %d'%use_GPU)


def main(image_in, mask_in, image_out, label_out, chan_ind=0, diam=70, pretrained_model='cyto2', distance=10):
    with tf.TiffFile(image_in) as tif:
        dapi = tif.pages[chan_ind].asarray()
    mask = tf.imread(mask_in)
    crop = dapi * mask
    print('>>> crop shape: ', crop.shape)

    # DEFINE CELLPOSE MODEL
    model = models.Cellpose(gpu=use_GPU, model_type='cyto2')
    channels= [[chan_ind, 0]]

    masks, flows, styles, diams = model.eval([crop], diameter=diam, channels=channels, do_3D=False)
    print('>>> masks shape: ', masks[0].shape)
    expanded = expand_labels(masks[0], distance=distance)
    tf.imwrite(label_out, expanded)
    stack = tf.imread(image_in)
    tf.imwrite(image_out, stack * mask, photometric='minisblack')


if __name__ == "__main__":
    fire.Fire(main)
