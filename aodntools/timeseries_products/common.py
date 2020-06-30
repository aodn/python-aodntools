"""Code shared by all timeseries product generating code"""


class NoInputFilesError(Exception):
    """Exception raised if there are no valid input files to aggregate"""
    pass


def check_file(nc, site_code, VoI):
    """
    Return list of errors found in the file.
    Checks applied:
    * Correct site_code
    * Variables of interest is present
    * Variables TIME, LATITUDE and LONGITUDE are present
    * NOMINAL_DEPTH is present as variable or attribute
    * file_version is FV01
    * if LATITUDE or LONIGUTDE are dimension, they have length 1
    * Global attributes time_deployment_start and time_deployment_end exist

    :param nc: xarray dataset
    :param site_code: code of the mooring site
    :param VoI: string. Variable of Interest
    :return: list of failed tests
    """

    attributes = list(nc.attrs)
    variables = list(nc.variables)
    allowed_dimensions = ['TIME', 'LATITUDE', 'LONGITUDE']
    error_list = []

    if site_code != nc.site_code:
        error_list.append('Wrong site_code: ' + nc.site_code)

    nc_file_version = nc.file_version
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    if 'TIME' not in variables:
        error_list.append('TIME variable missing')

    if 'LATITUDE' not in variables:
        error_list.append('LATITUDE variable missing')

    if 'LONGITUDE' not in variables:
        error_list.append('LONGITUDE variable missing')

    if VoI not in variables:
        error_list.append(VoI + ' not in file')
    else:
        VoIdimensions = list(nc[VoI].dims)
        if 'TIME' not in VoIdimensions:
            error_list.append('TIME is not a dimension')
        if 'LATITUDE' in VoIdimensions and len(nc.LATITUDE) > 1:
            error_list.append('more than one LATITUDE')
        if 'LONGITUDE' in VoIdimensions and len(nc.LONGITUDE) > 1:
            error_list.append('more than one LONGITUDE')
        for d in range(len(VoIdimensions)):
            if VoIdimensions[d] not in allowed_dimensions:
                error_list.append('not allowed dimensions: ' + VoIdimensions[d])
                break

    if 'NOMINAL_DEPTH' not in variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    if 'time_deployment_start' not in attributes:
        error_list.append('no time_deployment_start attribute')
    if 'time_deployment_end' not in attributes:
        error_list.append('no time_deployment_end attribute')

    return error_list