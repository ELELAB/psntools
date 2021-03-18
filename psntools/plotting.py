#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    plotting.py
#
#    Utilities to plot results from PSN analyses.
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



# third-party packages
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# psntools
from .defaults import CONFIG_PLOT_HEATMAP_NODES
from ._util import (
    get_config_plot,
    get_items
    )


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
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(df), n):
        yield (i,i + n)


#--------------------------------- Plot ------------------------------#



def plot_heatmap_nodes(df, \
                       outfile, \
                       nodes_per_page = 20, \
                       psn_labels = None, \
                       configfile = CONFIG_PLOT_HEATMAP_NODES):
    """Plot a heatmap with selected nodes on the x-axis and node
    metrics of interest for such nodes on the y-axis.
    """

    df = df.T

    #------------------------- Configuration -------------------------#


    config = get_config_plot(configfile = configfile)
    
    config_out = config["output"]
    
    config_heat, config_cbar, config_nan, \
    config_xaxis, config_yaxis = \
        get_items(config["plot"]["options"], \
                       ("heatmap", "colorbar", "nancells", \
                        "xaxis", "yaxis"), \
                       {})
    
    config_heatmap, config_annot = \
        get_items(config_heat, ("heatmap", "annot"), {})
    
    config_int = \
        get_items(config_cbar, ("interval",), {})[0]



    # get the colorbar ticks positions
    c_ticks = generate_ticks_positions(values = df.values.flatten(), \
                                       config = config_int)
            
    # get maximum and minimum value from the interval
    vmin, vmax = c_ticks[0], c_ticks[-1]


    with PdfPages(outfile) as pdf:

        chunk_ixs = get_chunk_indexes(df.T, nodes_per_page)

        for start_ix, end_ix in chunk_ixs:


            #-------------------- Data processing --------------------#


            # get a sub-dataframe of the original dataframe
            sub_df = df.iloc[:,start_ix:end_ix]

            # x-axis tick labels will be the column names
            xticklabels = sub_df.columns.values.tolist()
            
            # y-axis tick labels will be the row names
            yticklabels = psn_labels if psn_labels is not None \
                          else sub_df.index.values.tolist()
            
            # flatten the array so that we are dealing only with
            # a list of values
            values = sub_df.values.flatten()
            
            # drop NaN values
            yvalues = values[~np.isnan(values)]
            
            # get the cells where the value is NaN (nodes for which
            # the value is not available, i.e. if you PSNs with a
            # different number of nodes)
            nan_cells = np.argwhere(np.isnan(sub_df.values))


            #------------------------- Plot --------------------------#


            # create a new figure
            plt.figure()

            # generate the heatmap annotations
            annots = generate_heatmap_annotations(df = sub_df, \
                                                  config = config_annot)

            # generate the heatmap
            ax = sns.heatmap(sub_df, \
                             cbar = False, \
                             annot = annots[0], \
                             annot_kws = annots[1], \
                             vmin = vmin, \
                             vmax = vmax, \
                             center = (vmax+vmin)/2, \
                             **config_heatmap)

            # add a mask to the NaN cells
            generate_mask_nancells(ax = ax, \
                                   cells = nan_cells, \
                                   config = config_nan)

            # add the colorbar to the heatmap
            generate_colorbar(mappable = ax.get_children()[0], \
                              ticks = c_ticks, \
                              config = config_cbar)

            # set the x-axis
            set_axis(ax = ax, \
                     axis = "x", \
                     ticklabels = xticklabels, \
                     config = config_xaxis)

            # set the y-axis
            set_axis(ax = ax, \
                     axis = "y", \
                     ticks = None, \
                     ticklabels = yticklabels, \
                     config = config_yaxis)
        
            # save the figure to the PDF page
            pdf.savefig(**config_out)
            # clear the figure
            plt.clf()
            # close the current figure window
            plt.close()


def plot_2Dmap_nodes(dfs, \
                     x_metric, \
                     y_metric, \
                     hue_metric = None, 
                     psn_labels = None, \
                     selected_nodes = None):
    """Generate a 2D plot where two different node metric
    are plotted against each other to identify patterns or
    correlation between them.

    Parameters
    ----------
    dfs : `pandas.DataFrame`
    """

    #df = pd.DataFrame({"df1" : df1.values.flatten().tolist(), \
    #                   "df2" : df2.values.flatten().tolist()})
    
    df = pd.DataFrame()

    metrics = [x_metric, y_metric]
    
    if len(dfs) == 3:
        hue = hue_metric
        metrics += [hue_metric]

    # clear the figure
    plt.clf()
    # close the current figure window
    plt.close()

    for i, psn_col in enumerate(dfs[0].columns):

        dfs_to_concat = [d[psn_col] for d in dfs]

        psn_df = pd.concat(dfs_to_concat, \
                           axis = 1, \
                           keys = metrics)

        if selected_nodes is not None:
            psn_df = psn_df.loc[selected_nodes,:]

        ax = sns.scatterplot(data = psn_df, \
                             x = x_metric, \
                             y = y_metric, \
                             hue = hue_metric)

        psn_df["psn"] = psn_col

        df = pd.concat([df, psn_df], axis = 0)

        plt.savefig(f"{x_metric}_vs_{y_metric}_{psn_col}.pdf")

        # clear the figure
        plt.clf()
        # close the current figure window
        plt.close()


def plot_node_metric_distribution(self, \
                                  metric):

    df = self.get_nodes_df(metrics = [metric])
    sns.displot(df, x = metric)
    plt.savefig("prova.pdf")

