# Hourly Time Series Product (non-velocity)

- [Objective](#objective)
- [Input](#input)
- [Method](#method)
- [Output](#output)



## Objective

The code will provide aggregated files for each IMOS mooring site parameters (excluding current velocity data) binned into one hour intervals, excluding out-of-water data and records not flagged as “good” or “probably good” in the input files. QC flags will not be included. Statistics related to the averaging process will be stored as variables (standard deviation, minimum and maximum values, number of records binned). A version of the product can be generated with the inclusion of the non quality controlled variables (qc flag 0 -- no_QC_performed).


## Input


Files and variables to be aggregated by this code will meet the following requirements:

- File contains data from only one deployment of one instrument;
- File is a delayed-mode product, with file version label FV01;
- File contains, at the minimum, one of the variables to be aggregated (listed below), and variables TIME, LATITUDE, LONGITUDE, and either NOMINAL_DEPTH, or the global attribute instrument_nominal_depth;
- File is in netCDF format, compliant with CF-1.6 and IMOS-1.4 conventions;
- All files to be aggregated are from the same site, and have the same site_code attribute;
- Variables to be aggregated have TIME as their only dimension in the netCDF file  (or if LATITUDE and LONGITUDE are included as dimensions, they have size 1);
- Variables to be aggregated do not represent current velocity.
- The in-water data are bounded by the global attributes time_deployment_start and time_deployment_end.
- Variables have been quality-controlled and flagged using the IMOS standard flags (see [IMOS NetCDF Conventions v1.4](https://s3-ap-southeast-2.amazonaws.com/content.aodn.org.au/Documents/IMOS/Conventions/IMOS_NetCDF_Conventions.pdf), p43, Reference Table B)


## Method

### Input file validation

Before proceeding to the aggregation, each input file will be checked to ensure it meets the requirements (as specified above under Inputs). This will prevent errors during the aggregation process. Any input files that fail to meet the requirements will be excluded from the aggregation, and their URL listed in a global attribute rejected_files. 

### Dimensions

The dimensions of the resulting file  are determined as follows:

- `OBSERVATION`:  the total number of observation records, excluding out-of-the-water data, in all input files;
- `INSTRUMENT`: the number of instruments (i.e. number of files); 
- `string256`: a fixed dimension of length 256 for string variables.


### Variables

The product will aggregate the following variables (standard/long name in brackets):

`CPHL`, `CHLF`, `CHLU` (mass_concentration_of_inferred_chlorophyll_from_relative_fluorescence_units_in_sea_water) - considered as separate variables for now;  
`DEPTH` (depth);  
`DOX`   (volume_concentration_of_dissolved_molecular_oxygen_in_sea_water);   
`DOX1`  (mole_concentration_of_dissolved_molecular_oxygen_in_sea_water);   
`DOX2`  (moles_of_oxygen_per_unit_mass_in_sea_water);  
`DOXS`  (fractional_saturation_of_oxygen_in_sea_water);  
`DOXY`  (mass_concentration_of_oxygen_in_sea_water);  
`PAR`   (downwelling_photosynthetic_photon_flux_in_sea_water);  
`PRES`  (sea_water_pressure);  
`PRES_REL`  (sea_water_pressure_due_to_sea_water);  
`PSAL`  (sea_water_salinity);  
`TEMP`  (sea_water_temperature);  
`TURB`  (sea_water_turbidity);   
`TURBF`  (sea_water_turbidity_in_FTU).

Each Variable of Interest (VoI) is produced by selecting only the “good” and “probably good” data (QC flags 1 and 2), binning the values to a one-hour interval, and concatenating the resulting values into a single file. The resulting variables have dimension OBSERVATION. The VoI's ancillary_variables, in particular the corresponding quality-control flags, are not included. 

The code can generate a version of the product that also includes values flagged "no QC performed". In this case and for each variable, the percentage of QCed values will be indicated as a variable attribute.

The binning intervals will be one hour long, centred on the hour (i.e. HH:00:00). The binning method will depend on the type of variable being:

| Variables | Method | 
|-----------|--------|
| TEMP, CNDC, PRES, PRES_REL,  DEPTH, PSAL, DOXY, DOX, DOX1, DOX1_1, DOX1_2, DOX1_3, DOX2, DOX2_1, DOXS | Mean | 
| FLU2, CPHL, CHLU, CHLF, TURB, TURBF, PAR | Median | 

The binned `TIME` values are concatenated into a variable `TIME(OBSERVATION)`.  As input files can have overlapping or non-contiguous time coverage, the TIME variable may have duplicate values or gaps.

The binned `DEPTH` values from input files are concatenated into a variable `DEPTH(OBSERVATION)`. If not present, fill values are stored. `PRES`, `PRES_REL` and any other variables of interest are treated in a similar way. A variable will be excluded from the aggregated file if no valid values were found in any of the input files.

Additional metadata variables (ancillary variables) will be included, derived from the binning process and the source files: 
`VoI_std`: standard deviation per averaging bin;
`VoI_min`: minimum value per bin;
`VoI_max`: maximum value per bin;
`VoI_count`: number of observations in bin;

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

The output from a single run of the code will be an aggregated file of all available measurements of all non-velocity variable at one mooring site.

Two different files could be produced according to the values retained for the aggregation. One with values flagged as "Good_data" or "Probably_good_data" will have *"hourly-timeseries"* as product type. 
e.g. *IMOS_DWM-DA_STZ_20120421_EAC1_FV02_hourly-timeseries_END-20130823_C-20191007.nc* 

The product that also includes "No_QC_performed" values will have *"hourly-timeseries-including-non-QC"* as product type. 
e.g. *IMOS_DWM-DA_STZ_20120421_EAC1_FV02_hourly-timeseries-including-non-QC_END-20130823_C-20191007.nc*


The product will be delivered, in netCDF4 format, compliant with the CF-1.6 and IMOS-1.4 conventions, and
structured according to the [indexed ragged array representation](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_indexed_ragged_array_representation).
