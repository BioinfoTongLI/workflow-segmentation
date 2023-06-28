FROM nvidia/cuda:11.7.1-cudnn8-runtime-ubuntu22.04

RUN apt-get -y update && apt-get -y install --no-install-recommends procps python3 python3-dev python3-pip python-is-python3 libgomp1 ffmpeg libsm6 libxext6

COPY requirements.txt /tmp/requirements.txt

RUN pip install -r /tmp/requirements.txt --no-cache-dir

RUN mkdir /opt/.cellpose && \
    chmod -R a+rw /opt/.cellpose && \
    ln -s /opt/.cellpose $HOME/.cellpose &&\
    python -c "from cellpose import models;models.Cellpose(gpu=True, model_type='cyto2')"
