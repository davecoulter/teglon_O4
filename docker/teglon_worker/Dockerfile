FROM python:bullseye
LABEL description="Base Image for Teglon"
LABEL maintainer="Dave Coulter (dcoulter@ucsc.edu)"

SHELL ["/bin/bash", "-c"]

WORKDIR /app

# set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

COPY ./requirements.txt /app/requirements.txt

RUN apt-get -y update
RUN apt-get -y install default-libmysqlclient-dev
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN apt-get clean && rm -rf /opt/conda/pkgs/*
RUN rm -rf /var/lib/apt/lists/*