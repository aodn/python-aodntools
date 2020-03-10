# Velocity Aggregated Time Series Product

- [Objective](#objective)
- [Input](#input)
- [Method](#method)
- [Output](#output)




## Objective

This product provides aggregated U, V, and W velocity time-series files for each mooring site, without any interpolation or filtering, except for the exclusion of the out-of-water data. For the profiling (ADCP) instruments, the absolute depth of the measuring cell is calculated using the `DEPTH` measured at the instrument and the `HEIGHT_ABOVE_SENSOR`, The Quality Control (QC) flags are preserved.


## Input

The aggregation function will accept a list of input files, and the code of the mooring site (`site_code`), in addition to arguments that identify the path of input and output files.

The code aggregates variables and files that meet the following requirements:

- File contains data from only one deployment of one instrument;
- File is a delayed-mode, quality-controlled product (file version label “FV01”);
- File is compliant with CF-1.6 and IMOS-1.4 conventions;
- File contains, at the minimum, the components of current velocity (`UCUR`, `VCUR`), and variables `TIME`, `DEPTH`, `LATITUDE`, `LONGITUDE`, and `HEIGHT_ABOVE_SENSOR` in the case of ADCPs;
- All files to be aggregated are from the same site, and have the same `site_code` attribute; 
- Variables to be aggregated have `TIME` and (optionally) `HEIGHT_ABOVE_SENSOR` as their only dimensions (or if `LATITUDE` and `LONGITUDE` are included as dimensions, they have size 1);
- The in-water data are bounded by the global attributes `time_deployment_start` and `time_deployment_end`;
- The `TIME` variable has an attribute `seconds_to_middle_of_measurement` to indicate the offset from each recorded timestamp to the centre of the averaging period.


The code is able to access the input files either locally, or remotely via the OPeNDAP protocol. 

## Method

Generating function: 

```
usage: velocity_aggregated_timeseries.py [-h] -site SITE_CODE -files FILENAMES
                                         [-indir INPUT_DIR]
                                         [-outdir OUTPUT_DIR]
                                         [-download_url DOWNLOAD_URL]
                                         [-opendap_url OPENDAP_URL]

Concatenate X,Y,Z velocity variables from ALL instruments from ALL deployments
from ONE site

optional arguments:
  -h, --help                    show this help message and exit
  -site SITE_CODE               site code, like NRMMAI
  -files FILENAMES              name of the file that contains the source URLs
  -indir INPUT_DIR              base path of input files
  -outdir OUTPUT_DIR            path where the result file will be written. Default ./
  -download_url DOWNLOAD_URL    path to the download_url_prefix
  -opendap_url OPENDAP_URL      path to the opendap_url_prefix


```



### Input file validation

Before proceeding to the aggregation, each input file will be checked to ensure it meets the requirements (as specified above under Inputs). Any input files that fail to meet the requirements will be excluded from the aggregation, and their URL listed in a global attribute `rejected_files`.

### Dimensions

The dimensions of the resulting file  are determined as follows:

- `OBSERVATION`:    the total number of observation records, excluding out-of-the-water data, in all input files;
- `INSTRUMENT`:     the number of instruments (i.e. number of files);
- `strlen`:         a fixed dimension of length 256 for character array variables.

### Variables

The velocity variables are produced by flattening, then concatenating the arrays in each of the input files. The resulting variable has dimension `OBSERVATION`. Each variable’s ancillary_variables, in particular the corresponding quality-control flags, are also included, with dimension `OBSERVATION`. If the quality control variable is absent in any input file, the corresponding flags in the output file will  be set to 0 (“no QC performed”).

The variable `TIME` from input files is re-shaped to match the flattened velocity variables, then concatenated into a variable `TIME(OBSERVATION)`. This will result in a non-uniform time interval and repeated timestamps.

The `DEPTH` variable from input files is concatenated into a variable `DEPTH(OBSERVATION)`. In the case of ADCP instruments, the `HEIGH_ABOVE_SENSOR`  is converted to absolute depth by subtracting each of the height values from the depth measurements at the instrument. `DEPTH_quality_control`, if present, is also included. 

All output variables with the `INSTRUMENT` dimension are sorted in chronological order, and the input files aggregated chronologically, according to the global attribute time_deployment_start.

In order to keep track of the provenance of the aggregated file, accessory variables are created:

- `instrument_index(OBSERVATION)`: index [0:number of files] of the instrument used, referencing the `INSTRUMENT` dimension.
- `source_file(INSTRUMENT, strlen)`: URLs of the files used
- `instrument_id(INSTRUMENT, strlen)`: concatenated deployment_code, instrument and instrument_serial_number from the global attributes of each file
- `LATITUDE(INSTRUMENT)`: LATITUDE per instrument.
- `LONGITUDE(INSTRUMENT)`: LONGITUDE per instrument.
- `NOMINAL_DEPTH(INSTRUMENT)`: nominal depth per instrument, from the input file’s variable `NOMINAL_DEPTH` or global attribute instrument_nominal_depth.
- `SECONDS_TO_MIDDLE(INSTRUMENT)`:  offset from the timestamp to the middle of the measurement window for each deployment
- CELL_INDEX(OBSERVATION): index of the corresponding measuring cell



### Attributes

The variable attributes will comply with the IMOS metadata standards.

The global metadata will be a set of IMOS standard attributes. Fixed attributes are read from a [JSON file](../velocity_aggregated_timeseries_template.json) that contains the {key:value} pairs for each of them.

Attributes specific to each aggregated product, are added as follows:

- `site_code`: obtained from the input files (should be the same in all of them);
- `time_coverage_start`, `time_coverage_end`: set to the full range of TIME values in the aggregated file;
- `geospatial_vertical_min`, `geospatial_vertical_max`: set to the full range of DEPTH values in the aggregated file;
- `geospatial_lat_min`, `geospatial_lat_max` : set to the full range of LATITUDE values in the aggregated file;
- `geospatial_lon_min`, `geospatial_lon_max`: set to the full range of LONGITUDE values in the aggregated file;
- `date_created`: set to the date/time the product file is created;
- `history`: set to “<date_created>: Aggregated file created.”;
- `keywords`: set to a comma-separated list of the main variable names (“UCUR, VCUR, WCUR, DEPTH, AGGREGATED”);
- `lineage`: a statement about how the file was created, including a link to the code used; 
- `title`: "Long Timeseries Velocity Aggregated product: UCUR, VCUR, WCUR, DEPTH at <site_code>  between <time_coverage_start> and <time_coverage_end>"; 
- `rejected_files`: a list of URLs for files that were in the input files list, but did not meet the input requirements. 


## Output

The output from a single run of the code will be an aggregated file of all available current velocity measurements at one mooring site.

The product will be delivered, in netCDF4 format, compliant with the CF-1.6 and IMOS-1.4 conventions, and structured according to the [indexed ragged array representation](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_indexed_ragged_array_representation).


