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
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
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


def generate_ticks_positions(values, config):
    """Generate the positions that the ticks
    will have on a plot axis/colorbar/etc.
    """
    
    # get the configurations
    inttype = config.get("type")
    rtn = config.get("round_to_nearest")
    top = config.get("top")
    bottom = config.get("bottom")
    steps = config.get("steps")
    spacing = config.get("spacing")
    ciz = config.get("center_in_zero")
    
    # if no rounding has been specified and the
    # interval is continuous
    if not rtn and inttype == "continuous":
        # default to rounding to the nearest 0.5
        rtn = 0.5

    # if the maximum of the ticks interval has not
    # been specified
    if not top:
        if inttype == "discrete":
            # default top value is the maximum
            # of the values provided
            top = int(max(values))
        elif inttype == "continuous":
            # default top value is the rounded up 
            # maximum of the values
            top = np.ceil(max(values)*(1/rtn))/(1/rtn)

    # if the minimum of the ticks interval has not
    # been specified
    if not bottom:
        if inttype == "discrete":
            # default bottom value is the minimum
            # of the values provided
            bottom = int(min(values))
        elif inttype == "continuous":
            # default bottom value is the rounded down 
            # minimum of the values
            bottom = np.floor(min(values)*(1/rtn))/(1/rtn)

    # if the number of steps the interval should have
    # has not been specified
    if not steps:
        if inttype == "discrete":
            # default number of steps is lenght of
            # the integer range between the bottom
            # value and the top value
            steps = len(list(range(bottom, top)))
        elif inttype == "continuous":
            # default is 5 steps
            steps = 5

    # if the interval spacing has not been specified
    if not spacing:
        if inttype == "discrete":
            # default spacing is the one between two steps,
            # rounded up
            spacing = int(np.ceil(np.linspace(bottom, \
                                              top, \
                                              steps, \
                                              retstep = True)[1]))
        elif inttype == "continuous":
            # default spacing is the one between two steps
            spacing = np.linspace(bottom, \
                                  top, \
                                  steps, \
                                  retstep = True)[1]

    # if the two extremes of the interval coincide
    if top == bottom:
        # return only one value
        return np.array([bottom])

    # if the interval needs to be centered in zero
    if ciz:
        # get the highest absolute value
        absval = np.ceil(top) if top > bottom \
                 else np.floor(bottom)
        # top and bottom will be opposite numbers with
        # absolute value equal to absval
        top, bottom = absval, -absval
        # return an evenly spaced interval
        # between top and bottom values
        return np.linspace(bottom, top, steps)

    # return the ticks interval
    return np.arange(bottom, top + spacing, spacing)


def generate_mask_nancells(ax, cells, config):
    """Generate a mask to mark differently cells
    in a heatmap containing NaN values.
    """

    # for each NaN cell
    for y,x in cells:
        # add a rectangulat patch over the cell
        ax.add_patch(mpatches.Rectangle(xy = (x,y), **config))

    # return the ax the patches have been plotted on
    return ax


def generate_heatmap_annotations(df, config):
    """Generate the annotations to be plotted on 
    a heatmap (each cell is annotated with the
    corresponding value).
    """

    # if the configuration is empty
    if config == dict():
        # return a tuple filled with None values
        return (None, None)

    # get the configuration for the style of the annotations
    # and for the number of decimals to be kept in the
    # annotations
    annot = config.get("annot")
    ndecimals = config.get("ndecimals", 2)
    # if no annotation is requested, leave the dictionary
    # for the annotation properties empty
    annotkws = {}

    # if annotations are requested
    if config.get("annot"):
        # create a function to set the annotations 
        # to the desired precision
        annotfunc = lambda x : np.around(x, ndecimals) 
        # vectorize the function
        annottransform = np.vectorize(annotfunc)
        # create annotations for all cells of the heatmap
        annot = annottransform(df.values)
        # get the style of the annotations
        annotkws = config["style"]

    # return annotations and style
    return (annot, annotkws)


def generate_colorbar(mappable, \
                      ticks, \
                      config):
    """Generate a colorbar associated to a mappable
    (for example, a heatmap).
    """
 
    # plot the colorbar
    cbar = plt.colorbar(mappable, **config["colorbar"])

    # if there is an axis label (horizontal orientation)
    if config["label"].get("xlabel"):
        # set the colorbar label  
        cbar.ax.set_xlabel(**config.get("label"))
    # if there is an axis label (vertical orientation)
    elif config["label"].get("ylabel"):
        # set the colorbar label     
        cbar.ax.set_ylabel(**config.get("label"))

    # set the colorbar ticks and ticks labels
    # setting ticks on cbar.ax raises a UserWarning, but setting
    # tick labels does not
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(ticks, **config.get("ticklabels"))

    # return the colorbar
    return cbar


def set_axis(ax, \
             axis, \
             config, \
             ticks = None, \
             ticklabels = None):
    """Set up the x- or y-axis."""
    
    if ticks is None:
        if axis == "x":
            # default to the tick locations already present
            ticks = plt.xticks()[0]
        elif axis == "y":
            # default to the tick locations already present
            ticks = plt.yticks()[0]
    
    if ticklabels is None:
        # default to the string representations
        # of the ticks' locations
        ticklabels = [str(t) for t in ticks]

    # set the tick labels
    ticklabelsconfig = {}
    if config.get("ticklabels"):
        ticklabelsconfig = config["ticklabels"]
    
    # if it is the x-axis
    if axis == "x":    
        # if there is an axis label
        if config.get("label"):
            # set the axis label
            ax.set_xlabel(**config["label"])        
        # set the ticks
        ax.set_xticks(ticks = ticks)        
        # set the tick labels
        ax.set_xticklabels(labels = ticklabels, \
                           **ticklabelsconfig)
        if ticks != []:      
            # set the axis boundaries
            ax.spines["bottom"].set_bounds(ticks[0], \
                                           ticks[-1])
    
    # if it is the y-axis
    elif axis == "y":        
        # if there is an axis label
        if config.get("label"):
            # set the axis label
            ax.set_ylabel(**config["label"])        
        # set the ticks
        ax.set_yticks(ticks = ticks)        
        # set the tick labels
        ax.set_yticklabels(labels = ticklabels, \
                           **ticklabelsconfig)
        if ticks != []:       
            # set the axis boundaries
            ax.spines["left"].set_bounds(ticks[0], \
                                         ticks[-1])

    # if a configuration for the tick parameters was provided
    if config.get("tick_params"):
        # apply the configuration to the ticks
        ax.tick_params(axis = axis, \
                       **config["tick_params"])

    # return the axis
    return ax


def get_chunk_indexes(df, n):
    """Yield the indexes of successive n-sized 
    chunks from a dataframe.
    """
    for i in range(0, len(df), n):
        yield (i,i + n)


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