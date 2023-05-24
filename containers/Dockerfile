#FROM nvidia/cuda:11.2.0-cudnn8-runtime-ubuntu20.04
FROM nvidia/cuda:11.7.1-cudnn8-runtime-ubuntu22.04

RUN apt-get -y update && apt-get -y install --no-install-recommends procps python3 python3-dev python3-pip python-is-python3 libgomp1 ffmpeg libsm6 libxext6

RUN pip install numpy==1.24.3 scipy==1.10.1 tifffile==2023.4.12 fire==0.5.0 \
     dask-image==2023.3.0 imagecodecs==2023.3.16 \
     cellpose==2.2.1 cucim==23.4.1 scikit-image==0.20.0 cupy-cuda11x==12.0.0 ome-zarr==0.7.1 trackpy==0.6.1 shapely==2.0.1 aicsimageio[nd2]==4.11.0 \
     opencv-contrib-python-headless==4.7.0.72 --no-cache-dir

RUN mkdir /opt/.cellpose && \
    chmod -R a+rw /opt/.cellpose && \
    ln -s /opt/.cellpose $HOME/.cellpose &&\
    python -c "from cellpose import models;models.Cellpose(gpu=True, model_type='cyto2')"
