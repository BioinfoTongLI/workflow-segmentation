#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
import fire
import tifffile as tf
import numpy as np
from aicsimageio import AICSImage
from tifffile import imwrite
from cellpose import models, core
import sys


def cellpose_seg_3d(chunk, model, diam=30, chs=[2, 1]):
    masks, flows, styles, diams = model.eval(
            chunk, diameter=diam, channels=chs, do_3D=True)
    return masks


# segment cells with canonical cellpose API
def segment(stem:str, img_p:str, chs=[0, 0], s=0):
    chs_str=",".join([str(chs[0]), str(chs[1])])
    img = AICSImage(img_p)
    img.set_scene(s)
    # print(img.shape)

    with tf.TiffWriter(f"{stem}_serie_{s}_chs_{chs_str}_3d_label_array.tif",
                       append=True, bigtiff=True) as tif:
        for t in range(img.dims.T):
            # print(t)
            seg = cellpose_seg_3d(
                    img.get_image_dask_data("ZYX", T=t, C=1).compute(),
                    model, chs=chs)
            tif.write(seg)


if __name__ == "__main__":
    use_GPU = core.use_gpu()
    print('>>> GPU activated? %d'%use_GPU)
    model = models.Cellpose(gpu=True, model_type='cyto2')
    fire.Fire(segment)
