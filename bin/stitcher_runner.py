#! /usr/bin/env python
import argparse
import json
from pathlib import Path
from pprint import pprint
from typing import Any, Dict

from stitcher import stitch_masks

Report = Dict[str, Dict[str, Any]]


def make_dir_if_not_exists(dir_path: Path):
    if not dir_path.exists():
        dir_path.mkdir(parents=True)


def load_slicer_info_from_file(slicer_info_path: Path):
    with open(slicer_info_path, "r") as s:
        slicer_info = json.load(s)
    padding = slicer_info["padding"]
    overlap = slicer_info["overlap"]
    num_channels = slicer_info["num_channels"]
    return padding, overlap, num_channels


def read_slicer_info(
    img_dir: Path, slicer_info_path: Path, padding_str: str, overlap: int, no_cell
):
    img_dir_slicer_path = img_dir.joinpath("slicer_info.json")
    if slicer_info_path is not None:
        padding, overlap, num_channels = load_slicer_info_from_file(slicer_info_path)
    elif img_dir_slicer_path.exists():
        padding, overlap, num_channels = load_slicer_info_from_file(img_dir_slicer_path)
    else:
        overlap = overlap
        padding_list = [int(i) for i in padding_str.split(",")]
        padding = {
            "left": padding_list[0],
            "right": padding_list[1],
            "top": padding_list[2],
            "bottom": padding_list[3],
        }
        if no_cell:
            num_channels = 1
        else:
            num_channels = 2
    return padding, overlap, num_channels


def main(
    img_dir: Path,
    out_dir: Path,
    overlap: int,
    padding_str: str,
    no_cell: bool,
    slicer_info_path: Path,
):
    padding, overlap, num_channels = read_slicer_info(
        img_dir, slicer_info_path, padding_str, overlap, no_cell
    )
    make_dir_if_not_exists(out_dir)

    mask_out_name_template = "z{z:03d}_mask.ome.tiff"
    if num_channels == 1:
        no_cell = True

    mask_report = stitch_masks(
        img_dir, out_dir, mask_out_name_template, overlap, padding, no_cell
    )
    pprint(mask_report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", type=Path, required=True, help="path to directory with images"
    )
    parser.add_argument("-o", type=Path, required=True, help="path to output directory")
    parser.add_argument(
        "-v", type=int, default=0, help="overlap size in pixels, default 0"
    )
    parser.add_argument(
        "-p",
        type=str,
        default="0,0,0,0",
        help="image padding that should be removed, 4 comma separated numbers: left, right, top, bottom."
        + "Default: 0,0,0,0",
    )
    parser.add_argument(
        "--no_cell",
        action="store_true",
        default=False,
        help="images do not contain cell marker",
    )
    parser.add_argument(
        "--slicer_info",
        default=None,
        help="path to information from slicer."
        + "By default will check for it in input image directory",
    )
    args = parser.parse_args()

    main(args.i, args.o, args.v, args.p, args.no_cell, args.slicer_info)
