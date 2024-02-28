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
        img = (img - mean) / std
        img_norm = (img - np.min(img)) / (np.max(img) - np.min(img))

        # Rescale the dynamic range to 0 to 255
        stack[i, :, :] = (img_norm * 255).astype(np.uint8)
    return stack


# segment cells with canonical cellpose API
def segment(stem:str, img_p:str, s=0, C=0,
            Z_min=0, Z_max=-1, T_min=0, T_max=-1,
            normalize=False, **cellposeparams):
    chs_str=",".join([str(cellposeparams["channels"][0]),
                      str(cellposeparams["channels"][1])])
    img = AICSImage(img_p)
    img.set_scene(s)

    with tf.TiffWriter(f"{stem}_serie_{s}_chs_{chs_str}_C_{C}_3d_label_array.tif",
                       append=True, bigtiff=True) as tif:

        T_max = img.dims.T if T_max == -1 else T_max
        Z_max = img.dims.Z if Z_max == -1 else Z_max

        for t in range(T_min, T_max):
            stack = img.get_image_dask_data(
                "ZCYX", T=t, C=cellposeparams["channels"][0], Z=np.arange(Z_min, Z_max)
            )
            if bool(normalize):
                stack = normalize_image_stack_slice_by_slice(stack)
            print(f"Processing serie {s}, channel {chs_str}, time {t}, with shape {stack.shape}")

            segs, flows, styles, diams = model.eval(
                stack.compute(),
                **cellposeparams)
            tif.write(segs)


if __name__ == "__main__":
    use_GPU = core.use_gpu()
    print('>>> GPU activated? %d'%use_GPU)
    model = models.Cellpose(gpu=use_GPU, model_type='cyto2')
    fire.Fire(segment)
