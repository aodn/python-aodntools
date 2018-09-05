""" A wrapper class to write netcdf files directly from a dictionary.

written by: Hugo Oliveira ocehugo@gmail.com
"""

# TODO write tests
# TODO how we handle groups!?
# TODO cleanup the user precedence rules of fill_values
# TODO is_dim_consistent too complex
# TODO change_time too complex
# TODO check_var too complex
# TODO createVariables too complex
# TODO implement from_file/from_cdl/from_json kwarg!?

from copy import deepcopy

import netCDF4
import numpy as np


class NetCDFGroupDict(object):
    def __init__(self,
                 dimensions={},
                 variables={},
                 global_attributes={},
                 title='NetCDFGroupDict',
                 **kwargs):
        """ A dictionary to hold netCDF groups
            It consist of a generic class holding 3 different dictionaries:
            dimensions is a <key:int>  dict
            variables is <key:[str,class,list,dict,int]> dict
            global_attributes is a <key:int> dict

            This class has __add__ to growth variables/dims/global attrs
            and __sub__ to remove unwanted variables from
            other :NetCDFGroupDict: instances.

            Example:
                dmn = {'lon':360,'lat':210}
                var = {}
                var['water'] = {'type':'double','dims':['lat','lon']}
                w1 = NetCDFGroupDict(dimensions=dmn,variables=var)
                dmn2 = {'time':300,'lon':720,'lat':330}
                var2 = {}
                var2['temp'] = {'type':'double','dims':['time','lat','lon']}
                w2 = NetCDFGroupDict(dimensions=dmn2,variables=var2)
                w3 = w1+w2
                #w3.variables.keys() = ['water','temp']
                #w3.dimensions = {'time':300,'lon':360,'lat':210}
                w4 = w2-w1
                #w4.variables.keys() = ['temp']
                #w4.dimensions = {'lon':720,'lat':330,'time':300}
        """
        self.title = title
        self.dimensions = dimensions
        self.global_attributes = global_attributes
        self.variables = variables

        if self.is_dim_consistent:
            self.rdimensions = dict((x, True) if y is -1 else (x, False)
                                    for x, y in zip(self.dimensions.keys(),
                                                    self.dimensions.values()))
        else:
            raise TypeError("Correct the dimensions.")

        notstr = self.title.__class__ is not str
        if notstr:
            raise TypeError("Title is not a str object")

        self.check_dims(self.dimensions)
        self.check_var(self.variables)
        self.check_global_attributes(self.global_attributes)
        self.check_consistency(self.dimensions, self.variables)

    def __add__(self, other):
        self_copy = deepcopy(self)
        self_copy.dimensions.update(other.dimensions)
        self_copy.variables.update(other.variables)
        self_copy.global_attributes.update(other.global_attributes)
        self_copy.title = "{t1} + {t2}".format(t1=self.title, t2=other.title)
        return self_copy

    def is_dim_consistent(self):
        """Check if the variable dictionary
        is consistent with current dimensions"""
        checkdims = set()
        for k in self.variables.keys():
            try:
                for d in self.variables[k]['dims']:
                    checkdims.add(d)
            except KeyError:
                print("Variable %s missing dimension information `dims`" % k)

            except TypeError:
                if self.variables[k]['dims'] is None:
                    continue

                missing = ['dims']

                try:
                    self.variables['k']['vtype']
                except KeyError:
                    missing += ['type']

                try:
                    self.variables['k']['attr']
                except KeyError:
                    missing += ['attr']

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

    def search_time_in_vars(self):
        """Check all vars for specific time variables associated with them"""
        tvars = set()
        for v in self.variables:
            try:
                tvars.add(self.variables[v]['attr']['time']['value'])
            except KeyError:
                None

        isnone = tvars == set()
        if isnone:
            return None
        else:
            return tvars

    def change_time(self, var, timevar):
        """Change the time dimension associated with variable :var:
            :var: a list or str
                Ex: 'zeta'
                   ['zeta','u']
                   ['u','v']
                   ['Ptracer1','Ptracer2']
            :timevar: a list or str
                Ex: 'bry_time'
                   ['zeta_time','uv_time']
                   ['uv_time']
                   ['ptime1','ptime2']
        """

        if var.__class__ is str:
            var = [var]
        if timevar.__class__ is str:
            timevar = [timevar]

        if len(var) == 1 and len(timevar) > 1:
            raise ValueError('Invalid input')
        elif len(var) > 1 and len(timevar) == 1:
            timevar = [timevar for x in range(len(var))]

        for v, t in zip(var, timevar):
            vargroup = set(self.variables.keys())
            dimgroup = set(self.dimensions.keys())
            v_included = v in vargroup
            t_included = t in vargroup and t in dimgroup

            # varname should match dimname for time info
            if not t_included:
                raise ValueError('Time variable:', t, 'not present!')
            if not v_included:
                for k in self.variables.keys():
                    if v in k:
                        self.variables[k]['dims'][0] = t
                        self.variables[k]['attr']['time']['value'] = t
            else:
                self.variables[v]['dims'][0] = t
                self.variables[v]['attr']['time']['value'] = t

    @classmethod
    def check_dims(self, dimdict):
        """ Check the dictionary """
        for d in dimdict:
            notint = dimdict[d].__class__ is not int
            if notint:
                ValueError("Dimension %s is not an integer object" % d)

    @classmethod
    def check_var(self, vardict, name=None):
        """ Check if the dictionary have all the reuqired fields
        to be defined as variable"""
        if name is None:
            name = 'input'

        vkeys = vardict.keys()
        have_dims = 'dims' in vkeys
        have_type = 'type' in vkeys
        have_att = 'attr' in vkeys
        have_one = have_dims | have_type | have_att
        have_none = not have_one

        if have_none:
            for k in vkeys:
                self.check_var(vardict[k], name=k)

        if have_dims:
            notnone = vardict['dims'] is not None
            notlist = vardict['dims'] is not list
            if notnone and notlist:
                ValueError(
                    "Dim for %s should be a None or a list object" % name)

        if have_att:
            notdict = vardict['attr'] is not dict
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

    @classmethod
    def check_global_attributes(self, gadict):
        """ Check the dictionary """
        for g in gadict:
            notstr = gadict[g].__class__ is not str
            if notstr:
                ValueError("Global Attr %s is not an integer object" % g)

    @classmethod
    def check_consistency(self, dimdict, vdict):
        """ Check the dictionary """
        alldims = dimdict.keys()
        allvars = vdict.keys()
        for k in allvars:
            vardims = vdict[k]['dims']
            if vardims is None:
                continue
            else:
                missing = [x for x in vardims if x not in alldims]
                if missing:
                    raise ValueError("Variable %s has undefined dimensions: %s"
                                     % (k, missing))


class DatasetTemplate(NetCDFGroupDict):
    def __init__(self, *args, **kwargs):
        super(DatasetTemplate, self).__init__(*args, **kwargs)
        self.cattrs = set([
            'zlib', 'complevel', 'shuffle', 'fletcher32', 'contiguous',
            'chunksizes', 'endian', 'least_significant_digit'
        ])
        self.fill_aliases = set(
            ['fill_value', 'missing_value', 'FillValue', '_FillValue'])

    def set_output(self, outfile, mode='w', **kwargs):
        """Create the dataset """
        self.outfile = outfile
        self.ncobj = netCDF4.Dataset(self.outfile, mode=mode, **kwargs)

    def _create_var_opts(self, vdict):
        """Return a list with attribute names required for the creation of variable
        defined by :vdict: This include creation/special options like:
            `zlib`
            `least_significant_digit`
            `dimensions`
            etc"""
        vset = set(list(vdict.keys()))
        inside = vset.intersection(self.cattrs)
        aliases = vset.intersection(self.fill_aliases)

        if len(aliases) > 1:
            raise ValueError('You can only provide one missing value alias!')
        else:
            inside = inside.union(aliases)
        return list(inside)

    def createDimensions(self):
        """Create the dimensions on the netcdf file"""
        for dname, dval in zip(self.dimensions.keys(),
                               self.dimensions.values()):
            self.ncobj.createDimension(dname, dval)

    def createVariables(self, **kwargs):
        """Create all variables for the current class
        **kwargs are included here to overload all options for all variables
        like `zlib` and friends.
        """
        for v in self.variables.keys():
            varname = v #self.variables[v]['name']
            datatype = self.variables[v]['type']
            dimensions = self.variables[v]['dims']

            var_c_opts = {}
            cwargs = kwargs.copy()
            if dimensions is None:  # no kwargs in createVariable
                self.ncobj.createVariable(varname, datatype)
            else:
                var_c_keys = list(self._create_var_opts(self.variables[v]))

                var_c_opts = dict(
                    (x, self.variables[v][x]) for x in var_c_keys)

                ureq_fillvalue = [
                    x for x in cwargs.keys() if x in self.fill_aliases
                ]

                vreq_fillvalue = [
                    x for x in var_c_opts.keys() if x in self.fill_aliases
                ]

                var_c_opts.update(cwargs)

                # user precendence
                if (ureq_fillvalue and vreq_fillvalue):
                    [var_c_opts.pop(x) for x in vreq_fillvalue]
                    fv_val = [var_c_opts.pop(x) for x in ureq_fillvalue]
                    var_c_opts['fill_value'] = fv_val[-1]
                elif (ureq_fillvalue and not vreq_fillvalue):
                    fv_val = [var_c_opts.pop(x) for x in ureq_fillvalue]
                    var_c_opts['fill_value'] = fv_val[-1]
                else:
                    fv_val = [var_c_opts.pop(x) for x in vreq_fillvalue]
                    if fv_val:
                        var_c_opts['fill_value'] = fv_val[-1]

                self.ncobj.createVariable(
                    varname, datatype, dimensions=dimensions, **var_c_opts)

            if 'attr' in self.variables[v].keys():
                attrs = self.variables[v]['attr'].copy()
                for not_attr in self._create_var_opts(attrs):
                    attrs.pop(not_attr)

                for attname in attrs.keys():
                    var = self.ncobj.variables[varname]
                    value = np.array(attrs[attname]['value']).astype(
                        attrs[attname]['type'])
                    var.setncattr(attname, value)

    def createGlobalAttrs(self):
        """Add the global attributes for the current class"""
        for att in self.global_attributes.keys():
            self.ncobj.setncattr(att, self.global_attributes[att])

    def create(self, **kwargs):
        """Create in the dimensions/variable and attributes and fill with
        basic information"""
        self.createDimensions()
        self.createVariables(**kwargs)
        self.createGlobalAttrs()
        self.ncobj.sync()
        self.ncobj.close()
        self.ncobj = netCDF4.Dataset(self.outfile, 'a')
        pass
