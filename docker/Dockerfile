ARG DUSTMAPS_CONFIG_DIR
FROM python:bullseye
LABEL description="Base Image for Teglon"
LABEL maintainer="Dave Coulter (dcoulter@ucsc.edu)"

SHELL ["/bin/bash", "-c"]

WORKDIR /app

# set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

COPY ./requirements.txt /app/requirements.txt
COPY ./dust_maps_config.py /app/dust_maps_config.py

RUN apt-get -y update
RUN apt-get -y install default-libmysqlclient-dev
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN python ./dust_maps_config.py --data_dir {DUSTMAPS_CONFIG_DIR}
RUN apt-get clean && rm -rf /opt/conda/pkgs/*
RUN rm -rf /var/lib/apt/lists/*