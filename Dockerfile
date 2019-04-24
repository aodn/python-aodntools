FROM ubuntu:16.04

RUN apt-get update
ENV HDF5_DIR=/usr/include/hdf5/
RUN apt-get install -y git python-pip libxml2-dev libxslt-dev python-dev libhdf5-dev libnetcdf-dev libudunits2-dev

RUN pip install Cython wheel setuptools numpy

