""" A wrapper class to write netcdf files directly from a dictionary.

written by: Hugo Oliveira ocehugo@gmail.com
"""

# TODO how we handle groups!?
# TODO cleanup the user precedence rules of fill_values
# TODO is_dim_consistent too complex
# TODO check_var too complex
# TODO createVariables too complex
# TODO Allow for attribute types to be specified in JSON

import json
from collections import OrderedDict
from copy import deepcopy

import netCDF4

from .schema import validate_dimensions, validate_variables, validate_attributes


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
                var['water'] = {'type':'double','dimensions':['lat','lon']}
                w1 = NetCDFGroupDict(dimensions=dmn,variables=var)

                dmn2 = {'time':300,'lon':720,'lat':330}
                var2 = {}
                var2['temp'] = {'type':'double','dimensions':['time','lat','lon']}
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

        self.check_var(self.variables)
        self.check_consistency(self.dimensions, self.variables)

    def __add__(self, other):
        self_copy = deepcopy(self)
        self_copy.dimensions.update(other.dimensions)
        self_copy.variables.update(other.variables)
        self_copy.global_attributes.update(other.global_attributes)
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
        validate_attributes(value)
        self._global_attributes = value

    def is_dim_consistent(self):
        """Check if the variable dictionary is consistent with current dimensions"""
        checkdims = set()
        for k in self.variables.keys():
            try:
                for d in self.variables[k]['dimensions']:
                    checkdims.add(d)
            except KeyError:
                print("Variable %s missing dimension information `dims`" % k)

            except TypeError:
                if self.variables[k]['dimensions'] is None:
                    continue

                missing = ['dimensions']

                try:
                    self.variables['k']['vtype']
                except KeyError:
                    missing += ['type']

                try:
                    self.variables['k']['attributes']
                except KeyError:
                    missing += ['attributes']

                errstr = "Variable %s is missing information for: "
                for _ in missing:
                    errstr += '%s, '
                errtuple = tuple([k] + missing)
                print(errstr % errtuple)

        if checkdims != set(self.dimensions.keys()):
            print("Consistent dimensions are: %s" % checkdims)
            return False
        else:
            return True

    @classmethod
    def check_var(cls, vardict, name=None):
        """Check if the dictionary have all the required fields to be defined as variable """
        if name is None:
            name = 'input'

        vkeys = vardict.keys()
        have_dims = 'dimensions' in vkeys
        have_type = 'type' in vkeys
        have_att = 'attributes' in vkeys
        have_one = have_dims | have_type | have_att
        have_none = not have_one

        if have_none:
            for k in vkeys:
                cls.check_var(vardict[k], name=k)

        if have_dims:
            notnone = vardict['dimensions'] is not None
            notlist = vardict['dimensions'] is not list
            if notnone and notlist:
                ValueError(
                    "Dim for %s should be a None or a list object" % name)

        if have_att:
            notdict = vardict['attributes'] is not dict
            if notdict:
                ValueError("Attr for %s should be a dictionary object" % name)
        if have_type:
            notstr = vardict['type'].__class__ is not str
            nottype = vardict['type'].__class__ is not type
            notcompound = vardict['type'].__class__ is not netCDF4.CompoundType
            notvl = vardict['type'].__class__ is not netCDF4.VLType
            if notstr and nottype and notcompound and notvl:
                ValueError(
                    "Type for %s should be a string or type object" % name)

    @staticmethod
    def check_consistency(dimdict, vdict):
        """Check that all dimensions referenced by variables in :vdict: are defined in the :dimdict:"""
        alldims = dimdict.keys()
        allvars = vdict.keys()
        for k in allvars:
            vardims = vdict[k].get('dimensions')
            if vardims is None:
                continue
            else:
                missing = [x for x in vardims if x not in alldims]
                if missing:
                    raise ValueError("Variable %s has undefined dimensions: %s"
                                     % (k, missing))


class DatasetTemplate(NetCDFGroupDict):
    """Template object used for creating netCDF files"""

    def __init__(self, *args, **kwargs):
        super(DatasetTemplate, self).__init__(*args, **kwargs)
        self.cattrs = {'zlib', 'complevel', 'shuffle', 'fletcher32', 'contiguous', 'chunksizes', 'endian',
                       'least_significant_digit'}
        self.fill_aliases = {'fill_value', 'missing_value', 'FillValue', '_FillValue'}
        self.outfile = None
        self.ncobj = None

    @classmethod
    def from_json(cls, path):
        """Load template from a JSON file"""

        with open(path) as f:
            template = json.load(f, object_pairs_hook=OrderedDict)

        # TODO: validate_template(template)
        #       - just check that the only top-level properties are "dimensions", "variables" and "global_attribute"

        return cls(dimensions=template.get('dimensions'),
                   variables=template.get('variables'),
                   global_attributes=template.get('global_attributes')
                   )

    def _create_var_opts(self, vdict):
        """Return a list with attribute names required for the creation of variable
        defined by :vdict: This include creation/special options like:
            `zlib`
            `least_significant_digit`
            `dimensions`
            etc
        """
        vset = set(list(vdict.keys()))
        inside = vset.intersection(self.cattrs)
        aliases = vset.intersection(self.fill_aliases)

        if len(aliases) > 1:
            raise ValueError('You can only provide one missing value alias!')
        else:
            inside = inside.union(aliases)
        return list(inside)

    def update_dimensions(self):
        """Update the sizes of dimensions to be consistent with the arrays set as variable values, if possible.
        Otherwise raise ValueError. Also raise ValueError if a dimension that already has a non-zero size is not
        consistent with variable array sizes.
        """
        for name, var in self.variables.items():
            values = var.get('data')
            if values is None:
                continue

            var_shape = values.shape
            var_dims = var.get('dimensions', [])
            if len(var_shape) != len(var_dims):
                raise ValueError(
                    "Variable '{name}' has {ndim} dimensions, but value array has {nshape} dimensions.".format(
                        name=name, ndim=len(var_dims), nshape=len(var_shape)
                    )
                )

            for dim, size in zip(var_dims, var_shape):
                template_dim = self.dimensions[dim]
                if template_dim is None or template_dim == 0:
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

    def create_dimensions(self):
        """Create the dimensions on the netcdf file"""
        for dname, dval in zip(self.dimensions.keys(),
                               self.dimensions.values()):
            self.ncobj.createDimension(dname, dval)

    def create_variables(self, **kwargs):
        """Create all variables for the current class
        **kwargs are included here to overload all options for all variables
        like `zlib` and friends.
        """
        for varname, var in self.variables.items():
            datatype = var['type']
            dimensions = var['dimensions']
            cwargs = kwargs.copy()
            if dimensions is None:  # no kwargs in createVariable
                ncvar = self.ncobj.createVariable(varname, datatype)
            else:
                var_attr = var.get('attributes', {})
                var_c_keys = list(self._create_var_opts(var_attr))

                var_c_opts = {x: var_attr[x] for x in var_c_keys}

                ureq_fillvalue = [
                    x for x in cwargs.keys() if x in self.fill_aliases
                ]

                vreq_fillvalue = [
                    x for x in var_c_opts.keys() if x in self.fill_aliases
                ]

                var_c_opts.update(cwargs)

                # user precedence
                if ureq_fillvalue and vreq_fillvalue:
                    [var_c_opts.pop(x) for x in vreq_fillvalue]
                    fv_val = [var_c_opts.pop(x) for x in ureq_fillvalue]
                    var_c_opts['fill_value'] = fv_val[-1]
                elif ureq_fillvalue and not vreq_fillvalue:
                    fv_val = [var_c_opts.pop(x) for x in ureq_fillvalue]
                    var_c_opts['fill_value'] = fv_val[-1]
                else:
                    fv_val = [var_c_opts.pop(x) for x in vreq_fillvalue]
                    if fv_val:
                        var_c_opts['fill_value'] = fv_val[-1]

                ncvar = self.ncobj.createVariable(
                    varname, datatype, dimensions=dimensions, **var_c_opts)

            # add variable values
            if 'data' not in var:
                raise ValueError('No data specified for variable {varname}'.format(varname=varname))
            if var['data'] is not None:
                ncvar[:] = var['data']

            # add variable attributes
            if var.get('attributes'):
                attrs = var['attributes'].copy()
                for not_attr in self._create_var_opts(attrs):
                    attrs.pop(not_attr)
                ncvar.setncatts(attrs)

    def create_global_attributes(self):
        """Add the global attributes for the current class"""
        for att in self.global_attributes.keys():
            self.ncobj.setncattr(att, self.global_attributes[att])

    def to_netcdf(self, outfile, var_args={}, **kwargs):
        """
        Create a netCDF file according to all the information in the template.
        See netCDF4 package documentation for additional arguments.

        :param outfile: Path for the output file (clobbered by default if it already exists!)
        :param var_args: Additional arguments passed on to  netCDF4.Dataset.createVariables()
        :param kwargs: Additional arguments for netCDF4.Dataset()
        :return: None
        """
        self.outfile = outfile

        self.update_dimensions()
        if not self.is_dim_consistent():
            raise ValueError("Dimensions.")

        self.ncobj = netCDF4.Dataset(self.outfile, mode='w', **kwargs)
        self.create_dimensions()
        self.create_variables(**var_args)
        self.create_global_attributes()
        self.ncobj.sync()
        self.ncobj.close()
        self.ncobj = netCDF4.Dataset(self.outfile, 'a')
