FROM ubuntu:16.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    python-pip \
    python-dev \
    libxml2-dev \
    libxslt-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install \
    bump2version==0.5.10 \
    Cython \
    wheel \
    setuptools \
    numpy
