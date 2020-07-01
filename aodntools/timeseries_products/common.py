"""Code shared by all timeseries product generating code"""


class NoInputFilesError(Exception):
    """Exception raised if there are no valid input files to aggregate"""
    pass


def get_qc_variable_names(nc):
    """
    Return a list of the variables in a file with names ending in '_quality_control'

    :param nc: open xarray Dataset
    :return: list of variable names
    """
    varlist = list(nc.variables)
    return [v for v in varlist if v.endswith('_quality_control')]

def check_imos_flag_conventions(nc, varnames=None):
    """
    Check that the given QC variables are using the IMOS flag conventions (according to the
    quality_control_conventions variable attribute).
    If no variable names given, check all variables with names ending in '_quality_control'.
    Returns a list of error messages (empty if all ok).

    :param nc: open xarray Dataset
    :param varnames: list of variables to check, or None to check all QC variables
    :return: list of error messages
    """
    if varnames is None:
        varnames = get_qc_variable_names(nc)
    if isinstance(varnames, str):
        varnames = [varnames]

    # accept two variants on the convention name, used in versions 1.3 and 1.4 of the
    # IMOS NetCDF Conventions document
    accepted_conventions = {"IMOS standard flags", "IMOS standard set using the IODE flags"}
    errors = set()
    for var in varnames:
        if var not in nc.variables:
            errors.add('variable {var} not in file'.format(var=var))
            continue
        conventions = getattr(nc[var], 'quality_control_conventions', None)
        if conventions is None:
            errors.add('variable {var} missing quality_control_conventions'.format(var=var))
            continue
        if conventions not in accepted_conventions:
            errors.add('unexpected quality_control_conventions: "{conventions}"'.format(conventions=conventions))

    return sorted(errors)


def check_file(nc, site_code, variables_of_interest):
    """
    Check that a file meets the requirements for inclusion in a product.
    Return a list of errors
    
    Checks applied:
    * Correct site_code
    * file_version is FV01
    * Coordinate variables TIME, LATITUDE, LONGITUDE (& DEPTH for velocity) are present
    * NOMINAL_DEPTH is present as variable or attribute (instrument_nominal_depth)
    * At least one variable of interest is present
    * All variables of interest have only the allowed dimensions 
    * If LATITUDE or LONIGUTDE are dimension, they have length 1
    * Global attributes time_deployment_start and time_deployment_end exist
    * QC flag variables use the IMOS flag conventions

    :param nc: open xarray dataset
    :param site_code: code of the mooring site
    :param variables_of_interest: variable name or list of names to be included in the product
                                  (the string 'velocity' implies checks specific to velocity files,
                                  with variabls 'UCUR' AND 'VCUR' required)
    :return: list of failed tests
    """

    if variables_of_interest == 'velocity':
        variables_of_interest = {'UCUR', 'VCUR', 'WCUR'}
        required_variables = {'TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE', 'UCUR', 'VCUR'}
        allowed_dimensions = {'TIME', 'LATITUDE', 'LONGITUDE', 'HEIGHT_ABOVE_SENSOR'}
    else:
        required_variables = allowed_dimensions = {'TIME', 'LATITUDE', 'LONGITUDE'}

    if isinstance(variables_of_interest, str):
        variables_of_interest = {variables_of_interest}

    attributes = list(nc.attrs)
    variables = list(nc.variables)
    required_attributes = {'time_deployment_start', 'time_deployment_end'}
    error_list = []

    if site_code != nc.attrs.get('site_code', '[missing]'):
        error_list.append('Wrong site_code: ' + nc.site_code)

    nc_file_version = nc.attrs.get('file_version', '[missing]')
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    for var in required_variables:
        if var not in variables:
            error_list.append('{var} variable missing'.format(var=var))

    variables_to_aggregate = set(variables_of_interest) & set(variables)
    if not variables_to_aggregate:
        error_list.append('no variables to aggregate')

    for var in variables_to_aggregate:
        dims = set(nc[var].dims)
        if 'TIME' not in dims:
            error_list.append('no TIME dimension for {}'.format(var))
        if 'LATITUDE' in dims and len(nc.LATITUDE) > 1:
            error_list.append('more than one LATITUDE')
        if 'LONGITUDE' in dims and len(nc.LONGITUDE) > 1:
            error_list.append('more than one LONGITUDE')
        other_dims = dims - set(allowed_dimensions)
        if other_dims:
            error_list.append(
                'dimension(s) {other_dims} not allowed for {var}'.format(other_dims=other_dims, var=var)
            )

    if 'NOMINAL_DEPTH' not in variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    for attr in required_attributes:
        if attr not in attributes:
            error_list.append('no {} attribute'.format(attr))

    # check qc flag conventions for VoI and depth/pressure
    error_list.extend(check_imos_flag_conventions(nc))

    return error_list