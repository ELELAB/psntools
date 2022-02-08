#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-


# Standard library
import os
# Third-party packages
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



def plot_ccs_imins(csv_file,
                   outfile,
                   left,
                   right,
                   xspace,
                   bottom,
                   top,
                   yspace,
                   axvlines,
                   colors,
                   config):
    
    """Plot the size of the largest connected component as a function
    of the interaction strength cut-off.
    """

    #------------------------- Configuration -------------------------#


    # Get the configurations for all plot elements and
    # for the output
    config_xaxis = config["xaxis"]
    config_yaxis = config["yaxis"]
    config_legend = config["legend"]
    config_spines = config["spines"]
    config_axvlines = config["axvlines"]
    config_tickparams = config["tickparams"]
    config_output = config["output"]


    #----------------------------- Plot ------------------------------#


    # Load the data frame
    df = pd.read_csv(csv_file,
                     sep = ",",
                     index_col = 0)

    # Generate the plot
    ax = df.plot(kind = "line",
                 legend = False,
                 color = colors)


    #------------------------- Axes' limits --------------------------#


    # If the user did not provide axes' limits, set the defaults
    bottom = bottom if bottom is not None else 0
    top = top if top is not None else df.values.max()
    left = left if left is not None else 0.0
    right = right if right is not None else df.index.max()


    #---------------------------- X-axis -----------------------------#


    # Set and plot the x-axis label
    ax.set_xlabel("Interaction strength cut-off",
                  **config_xaxis["label"])

    # If no tick spacing on the x-axis has been specified
    if xspace is None:

        # Generate a linear interval of 10 steps between the
        # lowest and highest value and the spacing between
        # the steps
        x_ticks_interval, x_ticks_distance = \
            np.linspace(left, right, num = 10, retstep = True)

        # Round up the spacing to the nearest greater integer 
        x_ticks_distance = int(x_ticks_distance)

        # Generate the ticks for the x-axis
        x_ticks = np.arange(left,
                            x_ticks_distance*10,
                            x_ticks_distance)
    
    # If the tick spacing on the x-axis was specified
    else:

        # Generate a range between the lowest and highest
        # value, spaced according to what the user provided
        x_ticks = np.arange(left, right+xspace, xspace)

    # Set the ticks for the x-axis
    ax.set_xticks(x_ticks)

    # Plot x-ticks labels
    ax.set_xticklabels(x_ticks,
                       **config_xaxis["ticklabels"])


    #---------------------------- Y-axis -----------------------------#


    # Set and plot the y-axis label
    y_label = \
        "Size of the most populated\nconnected component\n" \
        "(number of nodes)"
    
    ax.set_ylabel(y_label, **config_yaxis["label"])

    # If no tick spacing on the y-axis has been specified
    if yspace is None:

        # Generate a linear interval of 10 steps between the
        # lowest and highest value and the spacing between
        # the steps
        y_ticks_interval, y_ticks_distance = \
            np.linspace(bottom, top, num = 10, retstep = True)

        # Round up the spacing to the nearest greater integer 
        y_ticks_distance = int(y_ticks_distance)

        # Generate the ticks for the y-axis
        y_ticks = np.arange(left,
                            y_ticks_distance*10,
                            y_ticks_distance)
    
    # If the tick spacing on the y-axis was specified
    else:

        # Generate a range between the lowest and highest
        # value, spaced according to what the user provided
        y_ticks = np.arange(bottom, top+yspace, yspace)

    # Set the ticks on the y-axis
    ax.set_yticks(y_ticks)

    # Plot y-ticks labels
    ax.set_yticklabels(y_ticks, **config_yaxis["ticklabels"])


    #---------------------------- Spines -----------------------------#

    
    # Set spine properties
    for spine, spine_settings in config_spines.items():
        ax.spines[spine].set(**spine_settings)

    # Set bounds for the x-axis
    ax.spines["bottom"].set_bounds(x_ticks[0], x_ticks[-1])

    # Set bounds for the y-axis
    ax.spines["left"].set_bounds(y_ticks[0], y_ticks[-1])


    #------------------------ Tick parameters ------------------------#


    # Set tick parameters for both axes
    ax.tick_params(**config_tickparams)


    #------------------------ Vertical lines -------------------------#


    if axvlines is not None:
        for axvline, sys in zip(axvlines, df.columns):
            plt.axvline(x = axvline,
                        ymin = y_ticks[0],
                        ymax = df.loc[axvline,sys]/y_ticks[-1],
                        **config_axvlines)


    #---------------------------- Legend -----------------------------#


    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles = handles,
              labels = labels,
              **config_legend)


    #---------------------------- Output -----------------------------#


    plt.savefig(outfile, **config_output)



if __name__ == "__main__":

    import argparse

    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add arguments
    i_help = "Input CSV file."
    parser.add_argument("-i", "--input-csv",
                        type = str,
                        required = True,
                        help = i_help)

    o_help = "Output file containing the plot."
    parser.add_argument("-o", "--output-plot",
                        type = str,
                        required = True,
                        help = o_help)

    colors_help = "Colors to be used for each system."
    parser.add_argument("-c", "--colors",
                        type = str,
                        default = None,
                        nargs = "*",
                        help = colors_help)

    imins_vlines_help = \
        "Interaction strength cut-offs where to " \
        "draw a vertical line (to highlight the " \
        "cut-off)."
    parser.add_argument("--imins-vlines",
                        type = str,
                        default = None,
                        nargs = "*",
                        help = imins_vlines_help)

    imin_min_help = \
        "Minimum interaction cut-off to be shown " \
        "on the x-axis. If no minimum is provided, it " \
        "will be taken from the data."
    parser.add_argument("--imin-min",
                        type = float,
                        default = None,
                        help = imin_min_help)

    imin_max_help = \
        "Maximum interaction cut-off to be shown " \
        "on the x-axis. If no maximum is provided, it " \
        "will be taken from the data."
    parser.add_argument("--imin-max",
                        type = float,
                        default = None,
                        help = imin_max_help)

    imin_spacing_help = \
        "Spacing for ticks on the x-axis showing the " \
        "interaction cut-offs. If no spacing is " \
        "provided, it will be taken directly from the data."
    parser.add_argument("--imin-spacing",
                        type = float,
                        default = None,
                        help = imin_spacing_help)

    cc_size_min_help = \
        "Minimum connected component size to be shown " \
        "on the y-axis. If no minimum is provided, it " \
        "will be 0."
    parser.add_argument("--cc-size-min",
                        type = int,
                        default = 0,
                        help = cc_size_min_help)

    cc_size_max_help = \
        "Maximum connected component size to be shown " \
        "on the y-axis. If no maximum is provided, it " \
        "will be taken from the data."
    parser.add_argument("--cc-size-max",
                        type = int,
                        default = None,
                        help = cc_size_max_help)

    cc_size_spacing_help = \
        "Spacing for ticks on the y-axis showing the " \
        "connected components' sizes. If no spacing is " \
        "provided, it will be taken directly from the data."
    parser.add_argument("--cc-size-spacing",
                        type = int,
                        default = None,
                        help = cc_size_spacing_help)

    # Parse the arguments
    args = parser.parse_args()

    # Convert the color definitions into matplotlib-compatible
    # colors
    colors = []
    for color in args.colors:
        # RGB triplet
        if color[0].isdigit():
            colors.append(tuple([float(c) for c in color.split(",")]))
        # Named color
        else:
            colors.append(color)

# Set the font properties
fname = "fonts/Nexa-Font/NexaLight.otf"
fp_ticklabels = fm.FontProperties(fname = fname, size = 24)
fp_axlabels = fm.FontProperties(fname = fname, size = 28)
fp_legend = fm.FontProperties(fname = fname, size = 24)

# Configuration of other plot aesthetics
config = \
    {"spines" : \
        {"bottom" : {"linewidth" : 1},
         "left" : {"linewidth" : 1},
         "top" : {"visible" : False},
         "right" : {"visible" : False}},
     "axvlines" : \
        {"linestyle" : "--",
         "linewidth" : 2,
         "color" : "#333333",
         "dashes" : (2,4),
         "dash_capstyle" : "round",
         "dash_joinstyle" : "round"},
     "legend" : \
        {"frameon" : False,
         "bbox_to_anchor" : (1, 0.5),
         "loc" : "center left",
         "prop" : fp_legend},
     "tickparams" : \
        {"axis" : "both", 
         "which" : "major",
         "pad" : 10,
         "grid_linewidth" : 1,
         "length" : 10},
     "xaxis" : \
        {"label" : {"fontproperties" : fp_axlabels,
                    "labelpad" : 30},
         "ticklabels" : {"fontproperties" : fp_ticklabels}},
     "yaxis" : \
        {"label" : {"fontproperties" : fp_axlabels,
                    "labelpad" : 30},
         "ticklabels" : {"fontproperties" : fp_ticklabels}},
     "output" : \
        {"bbox_inches" : "tight",
         "dpi" : 900,
         "transparent" : True}}

# Generate the plot
plot_ccs_imins(\
        csv_file = args.input_csv,
        outfile = args.output_plot,
        axvlines = args.imins_vlines,
        colors = args.colors,
        left = args.imin_min,
        right = args.imin_max,
        xspace = args.imin_spacing,
        bottom = args.cc_size_min,
        top = args.cc_size_max,
        yspace = args.cc_size_spacing,
        config = config)
