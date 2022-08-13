#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
All py functions needed in this pipeline
"""
import fire
from aicsimageio import AICSImage
from pathlib import Path
import numpy as np
import tifffile as tf


class Helper(object):

    def __init__(self, img_in_path: str, stem: str):
        self.img_in = Path(img_in_path.replace("\\", ""))
        assert self.img_in.exists()
        self.imgs = AICSImage(self.img_in, known_dims="CYX")
        self.stem = stem

    def extract_channels(self, *ch_ind):
        ch_names = self.imgs.channel_names
        print(ch_names, self.imgs.dims, ch_ind, len(ch_ind))
        if len(ch_ind) == 1:
            final_img_dim = "YX"
            ch_ind = ch_ind[0]
            target_ch_names = ["DAPI"]
        else:
            final_img_dim = "CYX"
            target_ch_names = np.array(ch_names)[np.array(ch_ind)]
        target_chs = self.imgs.get_image_dask_data(final_img_dim,
                T=0, Z=0, S=0, C=ch_ind)
        print(target_chs)

        tf.imwrite(
            "%s_%s.tif" % (self.stem, "_".join([str(i) for i in target_ch_names])),
            np.array(target_chs).squeeze(),
            photometric="minisblack",
            imagej=True,
            metadata={"axes": final_img_dim},
            # compression=args.compression
            # tile=(2**10, 2**10),
            # bigtiff=True,
            # planarconfig='separate',
        )


    # def expand(self, distance):
        # self.imgs


if __name__ == "__main__":
    fire.Fire(Helper)
