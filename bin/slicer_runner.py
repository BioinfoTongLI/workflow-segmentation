#! /usr/bin/env python
from pathlib import Path
import argparse
import json

import numpy as np
import tifffile as tif
import dask

from slicer import split_by_size


def path_to_str(path: Path):
    return str(path.absolute().as_posix())


def split_plane(in_path, pages, zplane, tile_size, overlap):
    channels = []
    for p in pages:
        channels.append(tif.imread(path_to_str(in_path), key=p))
    channel_stack = np.stack(channels, axis=0)
    plane_split, plane_tile_names, slicer_info = split_by_size(
        channel_stack, zplane, tile_size, tile_size, overlap
    )
    return plane_split, plane_tile_names, slicer_info


def save_slicer_info(out_dir: Path, slicer_info: dict):
    with open(out_dir.joinpath("slicer_info.json"), "w") as s:
        json.dump(slicer_info, s, sort_keys=False, indent=4)


def split_tiff(
    in_path: Path,
    out_dir: Path,
    tile_size: int,
    overlap: int,
    nzplanes: int,
    selected_channels: list,
):
    with tif.TiffFile(in_path) as TF:
        npages = len(TF.pages)

    for z in range(0, nzplanes):
        pages = []
        for c in selected_channels:
            page = c * nzplanes + z
            pages.append(page)
        print("pages", [p + 1 for p in pages], "/", npages)
        this_plane_tiles, this_plane_tile_names, slicer_info = split_plane(
            in_path, pages, z, tile_size, overlap
        )
        task = []
        for i, tile in enumerate(this_plane_tiles):
            out_path = path_to_str(out_dir.joinpath(this_plane_tile_names[i]))
            task.append(
                dask.delayed(tif.imwrite)(out_path, tile, photometric="minisblack")
            )

        dask.compute(*task, scheduler="threads")
    return slicer_info


def main(
    in_path: Path,
    out_dir: Path,
    tile_size: int,
    overlap: int,
    nzplanes: int,
    nchannels: int,
    selected_channels: list,
):

    if in_path.suffix not in (".tif", ".tiff"):
        raise ValueError("Only tif, tiff input files are accepted")

    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    if selected_channels is None:
        selected_channels = list(range(0, nchannels))
    else:
        selected_channels = [ch_id for ch_id in selected_channels if ch_id < nchannels]
    if len(selected_channels) > 2:
        msg = (
            f"Too many channels selected: {len(selected_channels)}."
            + " Select at most 2 channels."
        )
        raise ValueError(msg)

    slicer_info = split_tiff(
        in_path, out_dir, tile_size, overlap, nzplanes, selected_channels
    )
    slicer_info["num_channels"] = len(selected_channels)
    slicer_info["num_zplanes"] = nzplanes
    save_slicer_info(out_dir, slicer_info)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Split image into number of tiles")
    parser.add_argument("-i", type=Path, help="path to image file")
    parser.add_argument("-o", type=Path, help="path to output dir")
    parser.add_argument(
        "-s",
        type=int,
        default=1000,
        help="size of tile, default 1000x1000, if set to 0, then -n parameter used instead",
    )
    parser.add_argument(
        "-v", type=int, default=50, help="size of overlap, default 50, 0 = no overlap"
    )
    parser.add_argument(
        "--nzplanes", type=int, default=1, help="number of z-planes, default 1"
    )
    parser.add_argument(
        "--nchannels", type=int, default=1, help="number of channels, default 1, max 2"
    )
    parser.add_argument(
        "--selected_channels",
        type=int,
        nargs="+",
        default=None,
        help="space separated ids of channels you want to slice, e.g. 0 1 3, default all",
    )

    args = parser.parse_args()
    main(
        args.i,
        args.o,
        args.s,
        args.v,
        args.nzplanes,
        args.nchannels,
        args.selected_channels,
    )
