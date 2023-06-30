#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2023 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
import fire
import pandas as pd
import trackpy as tp
import matplotlib.pyplot as plt


def track(stem:str, centroids_p:str, cell_volumn_threshold:float=15000,
          search_range=25, memory=4,
          stub_threshold=15,
          pos_columns=['centroid-0', 'centroid-1', 'centroid-2']):
    df = pd.read_csv(centroids_p)
    cells_df = df[df["area"] > cell_volumn_threshold] # this is actually 3D volumn
    tp.linking.Linker.MAX_SUB_NET_SIZEj = 1000
    tracked_df = tp.link(cells_df, search_range=search_range,
                         memory=memory, pos_columns=pos_columns)
    plt.figure()
    tp.plot_traj(tracked_df, pos_columns=['centroid-1', 'centroid-0'],
                #  superimpose=raw_stack.get_image_data("ZYX", S=s, T=0).squeeze().transpose(2, 3, 0, 1).max(axis=0),
                #  label=True
                 )
    plt.savefig(f"{stem}_before_filtering.png")
    tracked_df.to_csv(f"{stem}_before_filtering.csv")
    filtered_tracked_df = tp.filtering.filter_stubs(tracked_df, threshold=stub_threshold)
    plt.figure()
    tp.plot_traj(filtered_tracked_df, pos_columns=['centroid-1', 'centroid-0'],
                #  superimpose=raw_stack.get_image_data("ZYX", S=s, T=0).squeeze().transpose(2, 3, 0, 1).max(axis=0),
                #  label=True
                 )
    plt.savefig(f"{stem}_after_filtering.png")
    filtered_tracked_df.to_csv(f"{stem}_after_filtering.csv")


if __name__ == "__main__":
    fire.Fire(track)
