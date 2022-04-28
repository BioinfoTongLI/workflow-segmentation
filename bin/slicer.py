#! /usr/bin/env python
import numpy as np

Image = np.ndarray


def get_tile(
    big_image: Image, hor_f: int, hor_t: int, ver_f: int, ver_t: int, overlap=0
):
    hor_f -= overlap
    hor_t += overlap
    ver_f -= overlap
    ver_t += overlap

    # check if tile is over image boundary
    left_check = hor_f
    top_check = ver_f
    right_check = hor_t - big_image.shape[-1]
    bot_check = ver_t - big_image.shape[-2]

    left_pad_size = 0
    top_pad_size = 0
    right_pad_size = 0
    bot_pad_size = 0

    if left_check < 0:
        left_pad_size = abs(left_check)
        hor_f = 0
    if top_check < 0:
        top_pad_size = abs(top_check)
        ver_f = 0
    if right_check > 0:
        right_pad_size = right_check
        hor_t = big_image.shape[-1]
    if bot_check > 0:
        bot_pad_size = bot_check
        ver_t = big_image.shape[-2]

    if len(big_image.shape) > 2:
        extra_dim_size = len(big_image.shape) - 2
        extra_dim_slice = tuple([slice(None)] * extra_dim_size)
        tile_slice = extra_dim_slice + (slice(ver_f, ver_t), slice(hor_f, hor_t))
        extra_dim_pad = tuple([(0, 0)] * extra_dim_size)
        padding = extra_dim_pad + (
            (top_pad_size, bot_pad_size),
            (left_pad_size, right_pad_size),
        )
    else:
        tile_slice = (slice(ver_f, ver_t), slice(hor_f, hor_t))
        padding = ((top_pad_size, bot_pad_size), (left_pad_size, right_pad_size))
    tile = big_image[tile_slice]
    if max(padding) > (0, 0):
        tile = np.pad(tile, padding, mode="constant")
    return tile


def split_by_size(arr: np.ndarray, zplane: int, tile_w: int, tile_h: int, overlap: int):
    """Splits image into tiles by size of tile.
    tile_w - tile width
    tile_h - tile height
    """
    x_axis = -1
    y_axis = -2
    arr_width, arr_height = arr.shape[x_axis], arr.shape[y_axis]

    x_ntiles = (
        arr_width // tile_w if arr_width % tile_w == 0 else (arr_width // tile_w) + 1
    )
    y_ntiles = (
        arr_height // tile_h if arr_height % tile_h == 0 else (arr_height // tile_h) + 1
    )

    tiles = []
    img_names = []

    # row
    for i in range(0, y_ntiles):
        # height of this tile
        ver_f = tile_h * i
        ver_t = ver_f + tile_h

        # col
        for j in range(0, x_ntiles):
            # width of this tile
            hor_f = tile_w * j
            hor_t = hor_f + tile_w

            tile = get_tile(arr, hor_f, hor_t, ver_f, ver_t, overlap)

            tiles.append(tile)
            name = "X{x:03d}_Y{y:03d}_Z{z:03d}.tif".format(
                x=j + 1, y=i + 1, z=zplane + 1
            )
            img_names.append(name)

    tile_shape = [tile_h, tile_w]
    ntiles = dict(x=x_ntiles, y=y_ntiles)
    padding = dict(left=0, right=0, top=0, bottom=0)
    if arr_width % tile_w == 0:
        padding["right"] = 0
    else:
        padding["right"] = tile_w - (arr_width % tile_w)
    if arr_height % tile_h == 0:
        padding["bottom"] = 0
    else:
        padding["bottom"] = tile_h - (arr_height % tile_h)
    info = dict(tile_shape=tile_shape, ntiles=ntiles, overlap=overlap, padding=padding)
    return tiles, img_names, info
