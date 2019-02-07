#!/usr/bin/env bash

env PYTHONPATH=/sw/chef/src/aodn-netcdf-tools/ python rottnest.py
ncdump rottnest.nc >rottnest.cdl
