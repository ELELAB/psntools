#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    _util.py
#
#    Generic utility functions.
#
#    *** PRIVATE *** not part of the public API
#
#    Copyright (C) 2020 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Matteo Tiberti 
#                       <matteo.tiberti@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.


# standard library
import os
# third-party packages
import matplotlib.font_manager as fm
import yaml
# psntools
from .defaults import CONFIG_PLOT_DIR



############################## PRIVATE ################################



def _recursive_traverse(data, \
                        actions, \
                        keys = None, \
                        func = None):
    """Recursively traverse a dictionary performing actions on
    its items. It is used to traverse and modify the dictionary
    of options retrieved when parsing a YAML configuration file.

    Actions than can be performed are:
    
    - pop_empty : removal of keys associated to `None` 
                  values.
    - substitute : substitution of values of specific
                   keys with a function of those values.
    - substitute_dict : substitution of an entire dictionary with
                        a the return value of a function taking 
                        as arguments the items in the dictionary.
    """

    # if data is a dictionary
    if isinstance(data, dict):
        
        # keys of items on which the actions will be
        # performed. If no keys are passed, all keys
        # in the dictionary will be considered.
        selkeys = keys if keys else data.keys()
        
        # for each key, value pair
        for k, v in list(data.items()):

            # if value is None
            if v is None:
                # if removal of None values has been requested
                if "pop_empty" in actions:
                    # remove the key from the dictionary
                    data.pop(k)
                    continue

            # if value is a dictionary
            elif isinstance(v, dict):
                # if the susbtitution concerns the entire
                # dictionary
                if "substitute_dict" in actions:
                    # if the key is in the selected keys
                    if k in selkeys:
                        # substitute the value with a function
                        # of the current value, with the function
                        # taking as inputs the key, values pairs
                        # in the dictionary
                        data[k] = func(**v)
                
                # recursively check the sub-dictionaries
                # of the current dictionary
                _recursive_traverse(data = v, \
                                    actions = actions, \
                                    keys = keys, \
                                    func = func)
        
            # if value is something else
            else:
                # if substitution of the current value
                # has been requested
                if "substitute" in actions:
                    # if the key is in the list of selected keys
                    if k in selkeys:
                        # substitute the value with a function
                        # of it, assuming the function takes
                        # the key and the value as arguments
                        data[k] = func(k,v)

        # return the dictionary
        return data


def _get_config_plot_version_1(config):
    """Parse the configuration file for plotting,
    version 1.
    """
    
    # create a copy of the configuration
    newconfig = dict(config)
    # substitute the font properties definitions
    # with the corresponding FontProperties instances
    _recursive_traverse(data = newconfig, \
                       actions = ["substitute_dict"], \
                       func = fm.FontProperties, \
                       keys = {"fontproperties"})
    # return the configuration
    return newconfig



############################### PUBLIC ################################



def get_config_plot(configfile):
    """Get the plotting configuration."""

    # get the name of the configuration file for plotting
    # (without extension)
    configfilename = os.path.splitext(os.path.basename(configfile))[0]

    # if the configuration file is a name without extension
    if configfile == configfilename:
        # assume it is a file in the directory where
        # configuration files for plotting are stored
        configfile = os.path.join(CONFIG_PLOT_DIR, \
                                  configfilename + ".yaml")
    # otherwise assume it is a file name/file path
    else:
        configfile = get_abspath(configfile)
    
    # load the configuration from the file
    config = yaml.safe_load(open(configfile, "r"))

    # check the version of the configuration file
    if config["version"] == 1:
        # return the configuration written in version 1 format
        return _get_config_plot_version_1(config = config)
    else:
        raise ValueError("Only version 1 configuration files " \
                         "are supported for now.")


def get_abspath(path):
    """Given a path, return its absolute path. Return
    `None` if the path given is `None`.
    """

    return os.path.abspath(path) if path is not None else path


def get_items(d, keys, default = None):
    """Similar to operator.itemgetter but defaults to
    a specific value if the key is not found (instead of
    throwing an exception).
    """

    return [d.get(k, default) for k in keys]