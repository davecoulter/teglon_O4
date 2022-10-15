FROM davecoulter/astromini3:1.5
LABEL description="Image for teglon debug"
LABEL maintainer="Dave Coulter (dcoulter@ucsc.edu)"

SHELL ["/bin/bash", "-c"]

WORKDIR /app

# set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

COPY ./requirements.api.txt /app/requirements.txt

# If there are problems, remove astLib and htop... they aren't necessary right now
RUN apt-get update

RUN pip install --upgrade pip

RUN pip install -r requirements.txt

RUN apt-get clean && rm -rf /opt/conda/pkgs/*

RUN rm -rf /var/lib/apt/lists/*
