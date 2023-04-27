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
from segment_anything import build_sam, SamAutomaticMaskGenerator
import tifffile as tf
import cv2
from typing import Any, Dict, List
import os
import numpy as np


def write_md_to_folder(masks: List[Dict[str, Any]], path: str) -> None:
    header = "id,area,bbox_x0,bbox_y0,bbox_w,bbox_h,point_input_x,point_input_y,predicted_iou,stability_score,crop_box_x0,crop_box_y0,crop_box_w,crop_box_h"  # noqa
    metadata = [header]
    for i, mask_data in enumerate(masks):
        mask = mask_data["segmentation"]
        # filename = f"{i}.png"
        # cv2.imwrite(os.path.join(path, filename), mask * 255)
        mask_metadata = [
            str(i),
            str(mask_data["area"]),
            *[str(x) for x in mask_data["bbox"]],
            *[str(x) for x in mask_data["point_coords"][0]],
            str(mask_data["predicted_iou"]),
            str(mask_data["stability_score"]),
            *[str(x) for x in mask_data["crop_box"]],
        ]
        row = ",".join(mask_metadata)
        metadata.append(row)
    metadata_path = os.path.join(path, "metadata.csv")
    with open(metadata_path, "w") as f:
        f.write("\n".join(metadata))
    return


 # save mask as label image
def save_mask(masks, path):
    mask = masks[0]["segmentation"]
    label = np.zeros((mask.shape[0], mask.shape[1]), dtype=np.uint32)
    for i, mask_data in enumerate(masks):
        if mask_data["area"] < 50000:
            mask = mask_data["segmentation"].astype(np.uint32)
            label += mask * i
    tf.imwrite(path, label)


def main(stem, image):
    mask_generator = SamAutomaticMaskGenerator(build_sam(checkpoint="/sam/model/sam_vit_h_4b8939.pth"))
    # img = cv2.imread(image)
    # print(type(img), img.shape)
    img = tf.imread(image)[:3, :, :]
    for i, ch in enumerate(img):
        min_val, max_val = np.min(ch), np.max(ch)
        img[i] = cv2.equalizeHist(cv2.convertScaleAbs((ch - min_val) / (max_val - min_val) * 255))
    # # print(img.shape)
    img = np.transpose(img, (1, 2, 0))
    masks = mask_generator.generate(cv2.convertScaleAbs(img, alpha=(255.0/65535.0)))
    save_mask(masks, os.sep.join([stem, "label.tif"]))
    write_md_to_folder(masks, stem)

if __name__ == "__main__":
    fire.Fire(main)
