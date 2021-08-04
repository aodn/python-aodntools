""" A wrapper class to write netcdf files directly from a dictionary.

written by: Hugo Oliveira ocehugo@gmail.com
"""

# TODO how we handle groups!?
# TODO cleanup the user precedence rules of fill_values
# TODO check_var too complex
# TODO createVariables too complex
# TODO Allow for attribute types to be specified in JSON

import json
from collections import OrderedDict
from copy import deepcopy
from warnings import warn

import netCDF4
import numpy as np

from .schema import validate_dimensions, validate_variables, validate_global_attributes, ValidationError


def metadata_attributes(attr):
    """
    Helper function to extract the metadata attributes (to be written to netCDF file) from a dictionary.
    The metadata attributes are those with keys beginning with a letter (upper or lower case).

    :param attr: Attribute dictionary (any dictionary-like class)
    :return: Dictionary of metadata attributes
    :rtype: same type as `attr`
    """
    meta = attr.__class__()
    for k, v in attr.items():
        if k[0].isalpha():
            meta[k] = v

    return meta


def special_attributes(attr):
    """
    Helper function to extract the special attributes (defining netCDF file structure) from a dictionary.
    Special attributes are those with keys beginning with an underscore character. The underscore is removed
    in the keys of the output dictionary.

    :param attr: Attribute dictionary (any dictionary-like class)
    :return: Dictionary of special attributes
    :rtype: same type as `attr`
    """
    meta = attr.__class__()
    for k, v in attr.items():
        if k[0].startswith('_'):
            meta[k[1:]] = v

    return meta


class NetCDFGroupDict(object):
    def __init__(self,
                 dimensions=None,
                 variables=None,
                 global_attributes=None,
                 **kwargs):
        """ A dictionary to hold netCDF groups
            It consist of a generic class holding 3 different dictionaries:
            dimensions is a <key:int>  dict
            variables is <key:[str,class,list,dict,int]> dict
            global_attributes is a <key:int> dict

            This class has __add__ to combine variables/dimensions/global attributes
            from :NetCDFGroupDict: instances.

            Example:
                dmn = {'lon':360,'lat':210}
                var = {}
                var['water'] = {'_datatype':'double','_dimensions':['lat','lon']}
                w1 = NetCDFGroupDict(dimensions=dmn,variables=var)

                dmn2 = {'time':300,'lon':720,'lat':330}
                var2 = {}
                var2['temp'] = {'_datatype':'double','_dimensions':['time','lat','lon']}
                w2 = NetCDFGroupDict(dimensions=dmn2,variables=var2)

                w3 = w1+w2
                #w3.variables.keys() = ['water','temp']
                #w3.dimensions = {'time':300,'lon':360,'lat':210}
        """
        self._dimensions = None
        self._variables = None
        self._global_attributes = None

        self.dimensions = dimensions or OrderedDict()
        self.variables = variables or OrderedDict()
        self.global_attributes = global_attributes or OrderedDict()

    def __add__(self, other):
        self_copy = deepcopy(self)
        self_copy.dimensions.update(other.dimensions)
        self_copy.global_attributes.update(other.global_attributes)
        for k, v in other.variables.items():
            if k in self_copy.variables:
                self_copy.variables[k].update(v)
            else:
                self_copy.variables[k] = v
        return self_copy

    @property
    def dimensions(self):
        """Property to store the dictionary mapping dimension names to their sizes."""
        return self._dimensions

    @dimensions.setter
    def dimensions(self, value):
        validate_dimensions(value)
        self._dimensions = value

    @property
    def variables(self):
        """Property to store dictionary of variables. Keys are variable names, values are dictionaries of variable
        properties (dimensions, type, attributes, etc...)
        """
        return self._variables

    @variables.setter
    def variables(self, value):
        validate_variables(value)
        self._variables = value

    @property
    def global_attributes(self):
        """Property to store dictionary of global attributes"""
        return self._global_attributes

    @global_attributes.setter
    def global_attributes(self, value):
        validate_global_attributes(value)
        self._global_attributes = value

    def validate_template_schema(self):
        """Validate the template against the full schema, raising ValidationError if not valid.
        This checks all the three dictionaries: dimensions, variables and global attributes.
        """
        validate_dimensions(self.dimensions)
        validate_variables(self.variables)
        validate_global_attributes(self.global_attributes)


class DatasetTemplate(NetCDFGroupDict):
    """Template object used for creating netCDF files"""

    STRUCTURAL_ATTRIBUTES = {'datatype', 'dimensions', 'zlib', 'complevel', 'shuffle', 'fletcher32', 'contiguous',
                             'chunksizes', 'endian', 'least_significant_digit'}
    FILL_VALUE_ALIASES = {'fill_value', 'FillValue'}

    def __init__(self, *args, **kwargs):
        super(DatasetTemplate, self).__init__(*args, **kwargs)
        self.outfile = None
        self.ncobj = None

    @classmethod
    def from_json(cls, path):
        """Load template from a JSON file"""

        with open(path) as f:
            try:
                template = json.load(f, object_pairs_hook=OrderedDict)
            except ValueError as e:
                raise ValueError("invalid JSON file '{path}' ({e})".format(path=path, e=e))

        return cls(dimensions=template.get('_dimensions'),
                   variables=template.get('_variables'),
                   global_attributes=metadata_attributes(template)
                   )

    def ensure_completeness(self):
        """Ensure that all variables have all the necessary information to create a netCDF file.
        If _data is not None and not a numpy array, attempt to convert it into a numpy array.
        If _datatype is missing, try to guess based on value assigned to _data.

        Raises ValidationError if
        * the "_data" value is missing for a variable; or
        * the "_datatype" value is missing and it can't be guessed from the "_data" value.
        """
        for name, var in self.variables.items():
            if "_dimensions" not in var:
                var["_dimensions"] = []

            if "_data" not in var:
                raise ValidationError("No data specified for variable '{name}'".format(name=name))
            if var["_data"] is not None and not isinstance(var["_data"], np.ndarray):
                var["_data"] = np.array(var["_data"])

            if "_datatype" not in var:
                datatype = getattr(var["_data"], 'dtype', None)
                if datatype is not None:
                    warn("Guessed data type '{datatype}' for variable '{name}'".format(datatype=datatype, name=name))
                    var["_datatype"] = datatype
                else:
                    raise ValidationError("No data type information for variable '{name}'".format(name=name))

    def ensure_consistency(self):
        """For each variable, ensure that the specified dimensions and data arrays are consistent with each other, and
        with the global dimensions defined for the template.

        For dimensions set to None in the template, update their size to match the size of corresponding data arrays,
        if possible to do so in a self-consistent way.

        Raises ValidationError if any variable has dimensions not defined in the template.

        Raises ValueError if unable to make dimensions consistent because
        * A data array's shape doesn't match the number of dimensions defined for the variable; or
        * A dimension that already has a defined size is not consistent with variable array sizes.

        Prerequisite: run ensure_completeness() first!
        """
        # TODO: check _datatype against type of array in _data?
        # TODO: check for "unused" dimensions?
        for name, var in self.variables.items():
            # check dimensions exist
            var_dims = var['_dimensions'] or []
            inconsistent_dims = set(var_dims).difference(self.dimensions)
            if inconsistent_dims:
                raise ValidationError("Variable '{name}' has undefined dimensions "
                                      "{inconsistent_dims}".format(name=name, inconsistent_dims=inconsistent_dims)
                                      )

            # if we have no data array, can't do any more
            values = var.get('_data')
            if values is None:
                continue

            # check number of dimensions
            var_shape = values.shape
            if len(var_shape) != len(var_dims):
                raise ValueError(
                    "Variable '{name}' has {ndim} dimensions, but value array has {nshape} dimensions.".format(
                        name=name, ndim=len(var_dims), nshape=len(var_shape))
                )

            # adjust dimension size if not already set
            for dim, size in zip(var_dims, var_shape):
                if self.dimensions[dim] is None:
                    self.dimensions[dim] = size

            # check that shape is now consistent
            template_shape = tuple(self.dimensions[d] for d in var_dims)
            if var_shape != template_shape:
                raise ValueError(
                    "Variable '{name}' has dimensions {var_dims} and shape {var_shape}, inconsistent with dimension "
                    "sizes defined in template {template_shape}".format(
                        name=name, var_dims=var_dims, var_shape=var_shape, template_shape=template_shape
                    )
                )

    def _create_var_opts(self, vname, vdict):
        """Return a dictionary of attributes required for the creation of variable
        defined by :vdict:. This includes creation/special options like `zlib`,
        `least_significant_digit`, `fill_value`, etc...

        Also check that only one fill value attribute is provided, and map it to
        the `fill_value` option for createVariable. If both `fill_value` and `FillValue`
        are specified:
        * If they have exactly the same value, print a warning and proceed;
        * If they have different values, raise ValueError.

        :vname: is only used for reporting the variable name in error/warning message
        """
        special_dict = special_attributes(vdict)
        struct_keys = self.STRUCTURAL_ATTRIBUTES.intersection(special_dict.keys())
        fill_aliases = self.FILL_VALUE_ALIASES.intersection(special_dict.keys())

        if len(fill_aliases) > 1:
            fill_dict = {f: special_dict[f] for f in fill_aliases}
            unique_values = set(fill_dict.values())
            msg = "Multiple fill value options specified for variable '{vname}': {fill_dict}.".format(
                vname=vname, fill_dict=fill_dict
            )
            if len(unique_values) == 1:
                warn(msg)
            else:
                raise ValueError(msg)

        var_opts = {k: special_dict[k] for k in struct_keys}
        if fill_aliases:
            var_opts['fill_value'] = special_dict[fill_aliases.pop()]
        return var_opts

    def create_dimensions(self):
        """Create the dimensions on the netcdf file"""
        for dname, dval in self.dimensions.items():
            self.ncobj.createDimension(dname, dval)

    def create_variables(self, **kwargs):
        """Create all variables for the current class
        **kwargs are included here to overload all options for all variables
        like `zlib` and friends.
        """

        # variable attributes to convert to the same type as the variable
        # datatype
        varattrs_to_convert_to_datatype = ['valid_min', 'valid_max', 'valid_range']

        for varname, varattr in self.variables.items():
            if not varattr['_dimensions']:  # no kwargs in createVariable
                ncvar = self.ncobj.createVariable(varname, varattr['_datatype'])
            else:
                var_c_opts = self._create_var_opts(varname, varattr)
                var_c_opts.update(kwargs)

                ncvar = self.ncobj.createVariable(varname, **var_c_opts)

            # add variable values
            if varattr['_data'] is not None:
                ncvar[:] = varattr['_data']

            # convert some variables attribute to variable datatype
            for varattr_to_convert in varattrs_to_convert_to_datatype:
                if varattr_to_convert in varattr.keys():
                    varattr[varattr_to_convert] = np.array(varattr[varattr_to_convert], dtype=varattr['_datatype'])

            # add variable attributes
            ncvar.setncatts(metadata_attributes(varattr))

    def create_global_attributes(self):
        """Add the global attributes for the current class"""
        for att in self.global_attributes.keys():
            self.ncobj.setncattr(att, self.global_attributes[att])

    def to_netcdf(self, outfile, var_args=None, **kwargs):
        """
        Create a netCDF file according to all the information in the template.
        See netCDF4 package documentation for additional arguments.

        :param outfile: Path for the output file (clobbered by default if it already exists!)
        :param var_args: Additional arguments passed on to  netCDF4.Dataset.createVariables()
        :param kwargs: Additional arguments for netCDF4.Dataset()
        :return: None
        """
        self.outfile = outfile
        _var_args = var_args or {}

        self.validate_template_schema()
        self.ensure_completeness()
        self.ensure_consistency()

        try:
            self.ncobj = netCDF4.Dataset(self.outfile, mode='w', **kwargs)
            self.create_dimensions()
            self.create_variables(**_var_args)
            self.create_global_attributes()
            self.ncobj.sync()
        except Exception:
            raise
        finally:
            self.ncobj.close()

        self.ncobj = netCDF4.Dataset(self.outfile, 'a')

    def get_data_range(self, varname):
        """
        Return the minimum and maximum values assigned to the given variable, ignoring fill values and NaNs.
        Raises ValueError if the variable doesn't exist, or has no valid data attached to it.

        :param varname: Name of the variable
        :return: (min, max) value
        :rtype: tuple
        """
        var = self.variables.get(varname)
        if var is None:
            raise ValueError("Variable '{varname}' does not exist".format(varname=varname))
        data = var.get('_data', [])

        # mask out the fillvalues
        fill_value = var.get('_FillValue') or var.get('_fill_value')
        if fill_value is None:
            mask = np.isnan(data)
        else:
            mask = np.logical_or(data == fill_value, np.isnan(data))
        data_masked = np.ma.array(data, mask=mask)
        if data_masked.mask.all():
            raise ValueError("No valid data for variable '{varname}'".format(varname=varname))

        return data_masked.min(), data_masked.max()
