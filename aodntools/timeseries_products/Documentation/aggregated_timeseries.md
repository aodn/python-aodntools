# Non-Velocity Aggregated Product

- [Objective](#objective)
- [Scope](#scope)
- [Input](#input)
- [Output](#output)
- [Method](#method)
- [Distribution](#distribution)



## Objective

This product project will provide:

A comprehensive set of aggregated time-series files for each mooring site and parameter (excluding current velocity), 
without any interpolation or filtering, except for the exclusion of the out-of-water data. The Quality Control (QC) 
flags will be preserved. All the (python) code used for the generation of the products openly available on GitHub.
Clear documentation on the code and the output file format.

## Scope

Input Files and variables to be aggregated by this project will meet the following requirements:

- File contains data from only one deployment of one instrument;
- File is a delayed-mode, quality-controlled product (file version label “FV01”);
- File is compliant with CF-1.6 and IMOS-1.4 conventions;
- File contains, at the minimum, the variable of interest (VoI), and variables TIME, LATITUDE, LONGITUDE, and either variable NOMINAL_DEPTH or global attribute instrument_nominal_depth;
- All files to be aggregated are from the same site, and have the same site_code attribute;
- Variables to be aggregated have TIME as their only dimension (or if LATITUDE and LONGITUDE are included as dimensions, they have size 1);
- Variables to be aggregated do not represent current velocity (these will be included in another product); and
- The in-water data are bounded by the global attributes time_deployment_start and time_deployment_end.

## Input

The aggregation function will accept a list of input files, and the name of the variable of interest (VoI).

The code will be able to aggregate variables and files that meet the following requirements:

- File contains data from only one deployment of one instrument;
- File is a delayed-mode product, with file version label FV01;
- File is compliant with CF-1.6 and IMOS-1.4 conventions;
- File contains, at the minimum, the VoI, and variables TIME, LATITUDE, LONGITUDE, and   NOMINAL_DEPTH;
- All files to be aggregated are from the same site, and have the same site_code attribute;
- Variables to be aggregated have TIME as their only dimension (or if LATITUDE and LONGITUDE are included as dimensions, they have size 1);
- Variables to be aggregated do not represent current velocity;
- The in-water data are bounded by the global attributes time_deployment_start and time_deployment_end.

This project will not make any corrections to or implement work-arounds for input files that do not meet these requirements. A report on files that cannot be aggregated will be produced. If any files are re-processed by the moorings facility to meet the requirements, these will be incorporated into future versions of the product.

While the code will be more generally applicable, the set of product files created as part of this project will be limited to the following variables (standard/long name in brackets):

- TEMP (sea_water_temperature);
- PSAL (sea_water_salinity);
- CPHL, CHLF, CHLU (mass_concentration_of_inferred_chlorophyll_from_relative_fluorescence_units_in_sea_water) - considered as separate variables for now;
- TURB (sea_water_turbidity);
- DOX1, DOX1_2, DOX1_3 (mole_concentration_of_dissolved_molecular_oxygen_in_sea_water);
- DOX2 (moles_of_oxygen_per_unit_mass_in_sea_water);
- DOXS (fractional_saturation_of_oxygen_in_sea_water).
- PAR (downwelling_photosynthetic_photon_flux)

## Method

Generating function: 

```
usage: aggregated_timeseries.py [-h] -var VARNAME -site SITE_CODE -files
                                FILENAMES [-path OUTPUT_PATH]

Concatenate ONE variable from ALL instruments from ALL deployments from ONE
site

optional arguments:
  -h, --help         show this help message and exit
  -var VARNAME       name of the variable to concatenate. Like TEMP, PSAL
  -site SITE_CODE    site code, like NRMMAI
  -files FILENAMES   name of the file that contains the source URLs
  -path OUTPUT_PATH  path where the result file will be written. Default ./

```




### Input file selection

The input files are accessed using the OPeNDAP protocol via the AODN THREDDS server.  A list of OPeNDAP URLs can be produced by the utility function get_moorings_url(), which queries an index of all published moorings files to obtain the files relevant for a given site and variable of interest (VoI). The list can also be filtered by feature type (e.g time-series or profile) temporal range, file version (FV01), processing mode (real-time/delayed) or a regular expression matched against the file name.

### Input file validation

Before proceeding to the aggregation, each input file will be checked to ensure it meets the requirements (as specified above under Inputs). This will prevent errors during the aggregation process. If any files fail to meet the requirements, a list of all non-compliant files is output, and the code stops. An option can be provided to continue processing only the compliant files.

### Dimensions

The dimensions of the resulting file  are determined as follows:
`OBSERVATION`:  the total number of observation records, excluding out-of-the-water data, in all input files;
`INSTRUMENT`: the number of instruments (i.e. number of files);
`string256`: a fixed dimension of length 256 for character array variables.

### Variables

The VoI is produced by sequentially concatenating the individual values in each of the input files. The resulting variable has dimension `OBSERVATION`. The VoI’s ancillary_variables, in particular the corresponding quality-control flags, are also included, with dimension OBSERVATION. If the quality control variable is absent in any input file, the corresponding flags in the output file will  be set to 0 (“no QC performed”) and a warning will be raised.

The variable `TIME` from input files is concatenated into a variable TIME(OBS). This could result in a non-uniform time interval.

The `DEPTH` variable from input files is concatenated into a variable `DEPTH(OBS)`. If not present, fill values are stored. `DEPTH_quality_control`, if present, is also included. Where the `DEPTH` variable is absent, the corresponding `DEPTH_quality_control` values are set to 9 (“missing value”).
All output variables with the OBS dimension are sorted in chronological order.

In order to keep track of the provenance of VoI in the aggregated file, accessory variables are created:

`instrument_index(OBSERVATION)`: index [0:number of files] of the instrument used, referencing the `INSTRUMENT` dimension.

`source_file(INSTRUMENT, string256)`: URLs of the files used.

`instrument_id(INSTRUMENT, string256)`: concatenated deployment_code, instrument and instrument_serial_number from the global attributes of each file.

`LATITUDE(INSTRUMENT)`: LATITUDE per instrument.

`LONGITUDE(INSTRUMENT)`: LONGITUDE per instrument.

`NOMINAL_DEPTH(INSTRUMENT)`: nominal depth per instrument.  


### Attributes

The variable attributes will comply with the IMOS metadata standards.

The global metadata will be a set of IMOS standard attributes. Fixed attributes are read from a JSON file that contains the {key:value} pairs for each of them. See the contents of this file at the end of this document.

Attributes specific to each aggregated product, are added as follows:

`site_code`, `local_time_zone`: obtained from the input files (should be the same in all of them);

`time_coverage_start`, `time_coverage_end`: set to the full range of `TIME` values in the aggregated file;

`geospatial_vertical_min/max`: set to the full range of `DEPTH` values in the aggregated file;

`geospatial_lat_min/max`: set to the full range of `LATITUDE` values in the aggregated file;

`geospatial_lon_min/max`: set to the full range of `LONGITUDE` values in the aggregated file;

`date_created`: set to the date/time the product file is created;

`history`: set to “<date_created>: Aggregated file created.”;

`keywords`: set to a comma-separated list of the main variable names (“<VoI>, TIME, DEPTH, LATITUDE, LONGITUDE”);

`lineage`: a statement about how the file was created, including a link to the code used, and any input parameters (except the input files, which are listed in the source_file variable);

`title`: "Aggregated time-series product: <VoI> at <sit````e_code> between <time_coverage_start> and <time_coverage_end>"

## Output


The output from a single run of the code will be an aggregated file of all available measurements of a single non-velocity variable at one mooring site.

The product will be delivered, in netCDF4 format, compliant with the CF-1.6 and IMOS-1.4 conventions, and
structured according to the [indexed ragged array representation](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_indexed_ragged_array_representation).

## Sample Files

Sample files could be downloaded from AODN THREDDS(add link) server.


## Distribution

### Discovery and access

All product files are made available via the AODN THREDDS server. This enables preview of file metadata and data subsetting via the OPeNDAP protocol, and file download via HTTP. The files will be discoverable under specific site/product directories.
For ANMN files, the path within the THREDDS catalogue follows the template "IMOS/ANMN/<sub-facility>/<site_code>/aggregated_timeseries/" for each site, for example http://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/NRS/NRSMAI/aggregated_timeseries/catalog.html
For ABOS, the products are not separated by site, so the path is "IMOS/ABOS/<sub-facility>/aggregated_timeseries/", e.g.  http://thredds.aodn.org.au/thredds/catalog/IMOS/ABOS/DA/aggregated_timeseries/catalog.html

### Maintenance

Product files will be re-processed (where necessary) to include any new input data.
