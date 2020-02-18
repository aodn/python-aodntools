# Gridded Time Series Product (Temperature)

- [Scope](#scope)
- [Input](#input)
- [Method](#method)
- [Output](#output)



## Scope

The code will provide aggregated files for TEMPERATURE at each IMOS mooring site binned into one hour intervals and linearly interpolated at defined target depths. The aggregation process excludes out-of-water data and records not flagged as “good” or “probably good”. QC flags will not be included. 


## Input

The input file is the hourly time series aggregated file for the site.

Only data in the hourly product will be included in the gridded product.   (see [here](https://github.com/aodn/python-aodntools/blob/master/aodntools/timeseries_products/Documentation/Hourly_timeseries.md) for detailed requirements)

## Method

### Dimensions

The dimensions of the resulting file  are determined as follows:

- `TIME`:  the 1-hour timestamps from the hourly time series product;
- `DEPTH`: the target depths where the temperature values have been interpolated. 


### Variables

The code will only aggregate the variable `TEMP` (sea water temperature) when two or more measurements exist at different depths at the same timestamp, and the separation of the sensors is not greater than a maximum distance specified for the site.  

The product will also contain: 

- `LATITUDE`: LATITUDE.
- `LONGITUDE`: LONGITUDE.
- `TEMP_count`: number of observations in the water column used for the interpolation at every timestamp

### Attributes

The variable attributes will comply with the IMOS metadata standards.

The global metadata will be a set of IMOS standard attributes. Fixed attributes are read from a [JSON file](https://github.com/aodn/python-aodntools/blob/master/aodntools/timeseries_products/gridded_timeseries_template.json) that contains the {key:value} pairs for each of them.

Attributes specific to each aggregated product, are added as follows:

- `site_code`: obtained from the input file;
- `time_coverage_start`, `time_coverage_end`: set to the full range of TIME values in the aggregated file;
- `geospatial_vertical_min/max`: set to the full range of DEPTH values in the aggregated file;
- `geospatial_lat_min/max`: set to the full range of LATITUDE values in the aggregated file;
- `geospatial_lon_min/max`: set to the full range of LONGITUDE values in the aggregated file;
- `date_created`: set to the date/time the product file is created;
- `history`: set to “<date> Hourly aggregated file created. <date>: Gridded file created.”;
- `included_values_flagged_as`: "Good_data, Probably_good_data"  
- `keywords`: set to a comma-separated list of the main variable names (“TEMP, DEPTH, HOURLY, GRIDDED”);
- `lineage`: a statement about how the file was created, including a link to the code used, and any input parameters (except the input files, which are listed in the source_file global_attribute variable)
- `title`: "Gridded Time Series Product: TEMP interpolated at <site_code> to fixed target depths at 1-hour time intervals, between <time_coverage_start> and <time_coverage_end> and <mininum target depth> and <maximum target depth> meters.


## Output

The output from a single run of the code will be an aggregated file of all available temperature measurements at one mooring site interpolated at target depths specific for each site.

The file is named according to [IMOS NETCDF file naming convention](https://s3-ap-southeast-2.amazonaws.com/content.aodn.org.au/Documents/IMOS/Conventions/IMOS_NetCDF_Conventions.pdf). 

e.g. *IMOS_ANMN-NSW_TZ_20091029_PH100_FV02_TEMP-gridded-timeseries_END-20190828_C-20200108.nc* 

The product will be delivered, in netCDF4 format, compliant with the CF-1.6 and IMOS-1.4 conventions, and structured according to the [Orthogonal multidimensional array representation](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_orthogonal_multidimensional_array_representation).
