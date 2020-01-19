# Hourly Time Series Product (non-velocity)

- [Scope](#scope)
- [Input](#input)
- [Method](#method)
- [Distribution](#distribution)
- [Output](#output)



## Scope

The code will provide aggregated files for TEMPERATURE at each IMOS mooring site binned into one hour intervals and interpolated at defined target depths. The aggregation process excludes out-of-water data and records not flagged as “good” or “probably good” in the input file. QC flags will not be included. 


## Input

The input file is the hourly timeseries aggregated file for the site. 

Files and the variable to be aggregated by this code will meet the following requirements:

- File contains data from only one deployment of one instrument;
- File is a delayed-mode product, with file version label FV01;
- File contains, at the minimum, one of the variables to be aggregated (listed below), and variables TIME, LATITUDE, LONGITUDE, and either NOMINAL_DEPTH, or the global attribute instrument_nominal_depth;
- File is in netCDF format, compliant with CF-1.6 and IMOS-1.4 conventions;
- All files to be aggregated are from the same site, and have the same site_code attribute;
- Sea Water Temperature values to be aggregated have TIME as their only dimension in the netCDF file  (or if LATITUDE and LONGITUDE are included as dimensions, they have size 1);
- The in-water data are bounded by the global attributes time_deployment_start and time_deployment_end.
- Temperature values have been quality-controlled and flagged using the IMOS standard flags (see [IMOS NetCDF Conventions v1.4](https://s3-ap-southeast-2.amazonaws.com/content.aodn.org.au/Documents/IMOS/Conventions/IMOS_NetCDF_Conventions.pdf), p43, Reference Table B)


## Method

### Dimensions

The dimensions of the resulting file  are determined as follows:

- `TIME`:  the 1-hour timestamps from the hourly-timeseries product;
- `DEPTH`: the target depth where the temperature values have been interpolated. 


### Variables

The product will only aggregate Sea Water Temperature `TEMP` when two or more lectures exist at different depth at any timestamp. 

In order to keep track of the provenance of VoIs in the aggregated file, accessory variables are created:

- `instrument_index(OBSERVATION)`: index [0:number of files] of the instrument used, referencing the INSTRUMENT dimension.
- `source_file(INSTRUMENT, string256)`: URLs of the files used
- `instrument_id(INSTRUMENT, string256)`: concatenated deployment_code, instrument and instrument_serial_number from the global attributes of each file
- `LATITUDE(INSTRUMENT)`: LATITUDE per instrument.
- `LONGITUDE(INSTRUMENT)`: LONGITUDE per instrument.
- `NOMINAL_DEPTH(INSTRUMENT)`: nominal depth per instrument, from the input file’s variable NOMINAL_DEPTH or global attribute instrument_nominal_depth.


### Attributes

The variable attributes will comply with the IMOS metadata standards.

The global metadata will be a set of IMOS standard attributes. Fixed attributes are read from a [JSON file](https://github.com/aodn/python-aodntools/blob/master/aodntools/timeseries_products/hourly_timeseries_template.json) that contains the {key:value} pairs for each of them.

Attributes specific to each aggregated product, are added as follows:

- `site_code`: obtained from the input files (should be the same in all of them);
- `time_coverage_start`, `time_coverage_end`: set to the full range of TIME values in the aggregated file;
- `geospatial_vertical_min/max`: set to the full range of DEPTH values in the aggregated file;
- `geospatial_lat_min/max`: set to the full range of LATITUDE values in the aggregated file;
- `geospatial_lon_min/max`: set to the full range of LONGITUDE values in the aggregated file;
- `date_created`: set to the date/time the product file is created;
- `history`: set to “<date_created>: Aggregated file created.”;
- `keywords`: set to a comma-separated list of the main variable names (“<VoI>, TIME, DEPTH, LATITUDE, LONGITUDE”);
- `lineage`: a statement about how the file was created, including a link to the code used, and any input parameters (except the input files, which are listed in the source_file variable)
- `title`: "Long time series Hourly Aggregated product: all available non-velocity variables at <site_code> between <time_coverage_start> and <time_coverage_end>"
- `included_values_flagged_as`: a list of the quality-control flags accepted for inclusion into the binning

## Output

The output from a single run of the code will be an aggregated file of all available temperature measurements at one mooring site interpolated at target depths specific for each site.

e.g. *IMOS_ANMN-NSW_TZ_20091029_PH100_FV02_TEMP-gridded-timeseries_END-20190828_C-20200108.nc* 

The product will be delivered, in netCDF4 format, compliant with the CF-1.6 and IMOS-1.4 conventions, and structured according to the [indexed ragged array representation](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_indexed_ragged_array_representation).
