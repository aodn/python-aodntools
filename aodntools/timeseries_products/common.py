"""Code shared by all timeseries product generating code"""


class NoInputFilesError(Exception):
    """Exception raised if there are no valid input files to aggregate"""
    pass


def check_file(nc, site_code, variables_of_interest):
    """
    Check that a file meets the requirements for inclusion in a product.
    Return a list of errors
    
    Checks applied:
    * Correct site_code
    * file_version is FV01
    * Coordinate variables TIME, LATITUDE and LONGITUDE are present
    * NOMINAL_DEPTH is present as variable or attribute (instrument_nominal_depth)
    * At least one variable of interest is present
    * All variables of interest have only the allowed dimensions 
    * If LATITUDE or LONIGUTDE are dimension, they have length 1
    * Global attributes time_deployment_start and time_deployment_end exist

    :param nc: open xarray dataset
    :param site_code: code of the mooring site
    :param variables_of_interest: variable name or list of names to be included in the product
    :return: list of failed tests
    """

    if isinstance(variables_of_interest, str):
        variables_of_interest = [variables_of_interest]
    attributes = list(nc.attrs)
    variables = list(nc.variables)
    allowed_dimensions = {'TIME', 'LATITUDE', 'LONGITUDE'}
    required_coordinates = {'TIME', 'LATITUDE', 'LONGITUDE'}
    required_attributes = {'time_deployment_start', 'time_deployment_end'}
    error_list = []

    if site_code != nc.attrs.get('site_code', '[missing]'):
        error_list.append('Wrong site_code: ' + nc.site_code)

    nc_file_version = nc.attrs.get('file_version', '[missing]')
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    for coord in required_coordinates:
        if coord not in variables:
            error_list.append('{coord} variable missing'.format(coord=coord))

    variables_to_aggregate = set(variables_of_interest) & set(variables)
    if not variables_to_aggregate:
        error_list.append('no variable to aggregate')

    for var in variables_to_aggregate:
        dims = set(nc[var].dims)
        if 'TIME' not in dims:
            error_list.append('no TIME dimension for {}'.format(var))
        if 'LATITUDE' in dims and len(nc.LATITUDE) > 1:
            error_list.append('more than one LATITUDE')
        if 'LONGITUDE' in dims and len(nc.LONGITUDE) > 1:
            error_list.append('more than one LONGITUDE')
        other_dims = dims - allowed_dimensions
        if other_dims:
            error_list.append(
                'dimension(s) {other_dims} not allowed for {var}'.format(other_dims=other_dims, var=var)
            )

    if 'NOMINAL_DEPTH' not in variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    for attr in required_attributes:
        if attr not in attributes:
            error_list.append('no {} attribute'.format(attr))

    return error_list