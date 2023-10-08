#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
import fire
import tifffile as tf
from aicsimageio import AICSImage
from skimage import exposure
from cellpose import models, core
import numpy as np


def cellpose_seg_3d(chunk, model, diam=30, cellprob_threshold=0.0, chs=[2, 1]):
    masks, flows, styles, diams = model.eval(
            chunk, diameter=diam, channels=chs,
            do_3D=True, cellprob_threshold=cellprob_threshold)
    return masks


def normalize_image_stack_slice_by_slice(stack):
    """
    Normalize a 3D image stack slice by slice using NumPy.

    Parameters:
    stack (numpy.ndarray): A 3D image stack.

    Returns:
    numpy.ndarray: The normalized 3D image stack.
    """
    # Get the shape of the stack
    z, y, x = stack.shape

    # Normalize the stack slice by slice
    for i in range(z):
        img = stack[i, :, :]
        mean = np.mean(img)
        std = np.std(img)
        stack[i, :, :] = (img - mean) / std
        p1, p2 = np.percentile(img, (0, 100))
        stack[i, :, :] = exposure.rescale_intensity(img, in_range=(p1, p2))
    return stack


# segment cells with canonical cellpose API
def segment(stem:str, img_p:str, chs=[0, 0], s=0,
            diameter=30, cellprob_threshold=0.0, C=0,
            Z_min=0, Z_max=1, T_min=0, T_max=-1, normalize=False):
    chs_str=",".join([str(chs[0]), str(chs[1])])
    img = AICSImage(img_p)
    img.set_scene(s)

    with tf.TiffWriter(f"{stem}_serie_{s}_chs_{chs_str}_C_{C}_3d_label_array.tif",
                       append=True, bigtiff=True) as tif:

        T_max = img.dims.T if T_max == -1 else T_max
        print(img.shape)

        for t in range(T_min, T_max):
            stack = img.get_image_dask_data(
                "ZYX", T=t, C=C, Z=np.arange(Z_min, Z_max)
            )
            if normalize:
                stack = normalize_image_stack_slice_by_slice(stack)
            seg = cellpose_seg_3d(
                    stack.compute(),
                    model, chs=chs, diam=diameter,
                    cellprob_threshold=cellprob_threshold)
            tif.write(seg)


if __name__ == "__main__":
    use_GPU = core.use_gpu()
    print('>>> GPU activated? %d'%use_GPU)
    model = models.Cellpose(gpu=use_GPU, model_type='cyto2')
    fire.Fire(segment)
