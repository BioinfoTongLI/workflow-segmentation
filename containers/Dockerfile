FROM nvidia/cuda:11.2.0-cudnn8-runtime-ubuntu20.04

RUN apt-get -y update && apt-get -y install --no-install-recommends procps python3 python3-dev python3-pip python-is-python3 libgomp1 ffmpeg libsm6 libxext6

RUN pip install numpy>=1.21 scipy tifffile fire \
     dask-image "dask[distributed]" imagecodecs \
     cellpose==1.0.2 cucim scipy scikit-image cupy-cuda112 ome-zarr==0.4.1 trackpy shapely aicsimageio dask-image \
     opencv-contrib-python-headless --no-cache-dir

RUN mkdir /opt/.cellpose && \
    chmod -R a+rw /opt/.cellpose && \
    ln -s /opt/.cellpose $HOME/.cellpose &&\
    python -c "from cellpose import models;models.Cellpose(gpu=True, model_type='cyto2')"
