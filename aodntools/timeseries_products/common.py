"""Code shared by all timeseries product generating code"""


class NoInputFilesError(Exception):
    """Exception raised if there are no valid input files to aggregate"""
    pass
