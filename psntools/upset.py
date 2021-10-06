#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    upsetpy.py
#
#    Python implementation of UpSet plots for the analysis of PSNs.
#
#    Copyright (C) 2021 Valentina Sora 
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



# Standard library
import functools
import itertools
# Third-party packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transforms
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
import seaborn.palettes as snspalettes



class UpSetPlot:

    """Class implementing a plotter for UpSet plots [1]_.

    .. [1] Lex, Alexander, et al. "UpSet: visualization of 
       intersecting sets." IEEE transactions on visualization
       and computer graphics 20.12 (2014): 1983-1992.
    """

    def __init__(self,
                 psngroup,
                 item_type = "hubs",
                 sort_by = "degree",
                 ascending = False,
                 min_card = None,
                 max_card = None,
                 min_deg = None,
                 max_deg = None,
                 **kwargs):
        """Initialize an UpSetPlot instance.

        Parameters
        ----------
        psngroup : `PSNGroup` instance
            PSNGroup instance.

        item_type : `str`, accepted values: `"hubs"`, `"edges"`,
                    default: `"hubs"`
            Whether to calculate intersections of hubs or edges.

        sort_by : `str`, accepted values: `"degree"`, `"cardinality"`,
                  default: `"degree"`
            How to sort the sets intersections:
            'degree' means sorting by the intersection
            degree (i.e. how many sets are intersecting),
            'cardinality' means sorting by the number of
            elements in each intersection.

        ascending : `bool`, default: `False`
            Whether to sort in ascending order.

        min_card : `int` or `None`, default: `None`
            Minimum size for an intersection to be considered.

        max_card : `int` or `None`, default: `None`
            Maximum size for an intersection to be considered.

        min_deg : `int` or `None`, default: `None`
            Minimum degree (i.e. minimum number of intersecting
            sets) for an intersection to be considered.

        max_deg : `int` or `None`, default: `None`
            Maximum degree (i.e. maximum number of intersecting
            sets) for an intersection to be considered.
        """

        # Get sets and set names
        self.sets, self.setnames = \
            self._get_sets(psngroup = psngroup,
                           item_type = item_type,
                           **kwargs)

        self._names2ints = \
            {v : k for k, v in enumerate(self.setnames)}

        # Get sets' intersections
        self.inters = \
            self._get_inters(psngroup = psngroup,
                             item_type = item_type,
                             min_card = min_card,
                             max_card = max_card,
                             min_deg = min_deg,
                             max_deg = max_deg,
                             sort_by = sort_by,
                             ascending = ascending,
                             **kwargs)

        # Initialize all the axes that will be created
        # in the plot to None
        self.ax_matrix = None
        self.ax_sizebars = None
        self.ax_setbars = None

        # Get available palettes and available colors
        self.available_palettes = self._get_available_palettes()
        self.available_colors = list(colors.cnames.keys()) 

    
    def _get_available_palettes(self):
        """Return the list of available color
        palettes (`matplotlib` and `seaborn`).
        """
        
        # List all the names of the available Matplotlib palettes
        mpl_palettes = list(plt.colormaps())
        
        # List all the names of the available Seaborn palettes
        sns_palettes = list(snspalettes.SEABORN_PALETTES.keys())
        
        # TODO: come back later to see why I wanted to remove
        # the jet palette so badly

        # mpl_palettes.remove("jet")

        # Return all available palettes
        return mpl_palettes + sns_palettes


    def _get_sets(self, 
                  psngroup,
                  item_type,
                  **kwargs):

        sets = []
        setnames = []

        for psn_label, psn in psngroup.psns.items():
            if item_type == "hubs":
                sets.append(set(psn.get_hubs(**kwargs).keys()))
            elif item_type == "edges":
                sets.append(set(psn.get_edges(**kwargs).keys()))
            
            setnames.append(psn_label)

        return sets, setnames


    def _get_inters(self,
                    psngroup,
                    item_type,
                    min_card,
                    max_card,
                    min_deg,
                    max_deg,
                    sort_by,
                    ascending,
                    **kwargs):
        """Process the sets' intersections for plotting.
        """

        # Minimum size of an intersection is 0
        min_card = min_card if min_card is not None else 0
        
        # Maximum size of an intersection is the sum of the
        # elements in all sets
        sumsets = sum([len(s) for s in self.sets])
        max_card = max_card if max_card is not None else sumsets
        
        # Minimum degree for a meaningful intersection
        # is 2 (two sets intersecting)
        min_deg = min_deg if min_deg is not None else 2
        
        # Maximum degree for an intersection is the maximum
        # number of possible intersecting sets
        max_deg = max_deg if max_deg is not None else len(self.sets)

        # Empty list to store the intersections
        inters = []

        if item_type == "hubs":
            common_items = psngroup.get_common_hubs(**kwargs)
        elif item_type == "edges":
            common_items = psngroup.get_common_edges(**kwargs)

        # For each combination of PSNs
        for inter_label, inter in common_items.items():

            # Get the number of items in the intersection
            inter_card = \
                [len(inter_items) for inter_items in inter.values()][0]

            # Store the intersection only if it is of the correct size
            if inter_card >= min_card and inter_card <= max_card:

                inter_label = tuple(\
                    [self._names2ints[psn_label] for psn_label in \
                     inter_label])

                # We store in the list a tuple for each 
                # intersection, containing the numbers of the
                # intersecting sets as a tuple and the intersection
                # itself as a set  
                inters.append((inter_label, inter_card))

        # Sort the intersections by degree
        if sort_by == "degree":
            # If sorting is ascending
            if ascending:
                # Sort the intersections by degree and return them
                return sorted(inters,
                              key = lambda x: len(x[0]))
            # Otherwise
            else:
                # Reverse the sorting and return the intersections
                return sorted(inters,
                              key = lambda x: len(x[0]),
                              reverse = True)

        # Sort the list by cardinality       
        elif sort_by == "cardinality":
            # If sorting is ascending
            if ascending:
                # Sort the interactions by cardinality and return them
                return sorted(inters,
                              key = lambda x: len(x[1]))
            # Otherwise
            else:
                # Reverse-sort the interactions by cardinality and
                # return them
                return sorted(inters,
                              key = lambda x: len(x[1]),
                              reverse = True)
        
        # Raise an error if some other value was passed 
        else:
            errstr = \
                f"Unrecognized value for 'sort_by': {sort_by}."
            raise ValueError(errstr)
    
    
    def _draw_matrix(self,
                     inter_color,
                     noninter_color,
                     linestyle,
                     yticklabels,
                     custom_fpath):  
        """Draw the matrix showing intersecting
        sets as points connected by lines.
        """


        #---------------------- Data processing ----------------------#


        # Get the combinations generating the intersections
        combos = [inter[0] for inter in self.inters]

        # Generate integer 'names' to be assigned to the combinations
        names = list(range(len(self.inters)))
        
        # Assign an integer 'name' to each combination progressively
        combos2names = dict(zip(combos, names))


        #----------------------- Matrix points -----------------------#


        # Get the x coordinates of the points in the matrix by
        # repeating the 'name' of each combination as many times
        # as the size of the combination itself
        x_interspoints = \
            list(itertools.chain.from_iterable(\
                    [np.repeat(combos2names[key], len(key)).tolist() \
                     for key in combos])) 
        
        # Get the y coordinates of the points in the matrix by chaining
        # all the tuples representing the combinations in a list
        y_interspoints = [i for sub in combos for i in sub]

        # Generate a num.intersections x num.sets grid filled with 1s
        points_grid = np.ones((len(self.inters), len(self.setnames)))
        
        # Zero all the elements corresponding to intersecting sets
        points_grid[x_interspoints,y_interspoints] = 0.0
        
        # Get the coordinates of points representing
        # non-intersecting sets
        x_noninterspoints, y_noninterspoints = \
            np.nonzero(points_grid)

        # Generate the final sets of coordinates
        x_data = x_interspoints + x_noninterspoints.tolist()
        y_data = y_interspoints + y_noninterspoints.tolist()


        #--------------------------- Colors --------------------------#


        # Deal with the user-provided color representation
        is_rgbtuple = isinstance(inter_color, tuple)
        is_string = isinstance(inter_color, str)
        is_hex = False

        # If the color is a HEX code
        if is_string and inter_color.startswith("#"):
            is_hex = True
        
        # If the color is either a RGB tuple or a HEX code
        if is_rgbtuple or is_hex:
            # Color of the points: repeat it for the number of points
            # Color of the vertical lines: repeat it for the number
            #                              of combinations
            facecolor = [inter_color] * len(x_interspoints)
            linecolor = [inter_color] * len(names)

        # If it is a color palette or a named color
        elif is_string and not is_hex:

            # If it is the name of an available palette
            if inter_color in self.available_palettes:
                
                # Get the palette
                palette = sns.color_palette(inter_color, len(names))
                
                # Discretize the palette by getting as many colors
                # as there are combinations
                facecolor = \
                    list(\
                        itertools.chain.from_iterable(\
                            [[palette[i]]*len(key) \
                             for i, key in enumerate(combos)]))
                
                # The linecolor will be the palette
                linecolor = palette

            # If it is the name of a named color
            elif inter_color in self.available_colors:
                # Color of the points: repeat it for the number 
                #                      of points
                # Color of the vertical lines: repeat it for the
                #                              number of combinations
                facecolor = [inter_color] * len(x_interspoints)
                linecolor = [inter_color] * len(names)
            
            # If it is unrecognized
            else:
                # Raise an error
                errstr = \
                    f"{inter_color} is neither a name of a " \
                    f"color palette nor a named color."        
                raise ValueError(errstr)

        # The final list of colors for all points includes both the
        # color of points representing intersecting sets and the
        # color of points representing non-interacting sets
        facecolor += ([noninter_color] * len(x_noninterspoints))


        #--------------------- Matrix generation ---------------------#


        # The matrix is basically rendered as a scatterplot
        self.ax_matrix.scatter(x_data,
                               y_data,
                               s = 10,
                               facecolor = facecolor)


        #----------------------- Vertical lines ----------------------#


        # Get the points representing the extremes of the
        # vertical lines to be drawn
        points_vlines = [[(i, i), (tup[0], tup[-1])] \
                         for i, tup in enumerate(combos)]    
        
        # Plot the lines. Their color should match that of the points
        # they connect
        for points_pair, color in zip(points_vlines, linecolor):
            
            # Get the extremes of the line
            p1, p2 = points_pair
            
            # Plot the line between the two extremes
            self.ax_matrix.plot(p1,
                                p2,
                                color = color,
                                linestyle = linestyle)

        # Dumb patch to align the vertical ax of the matrix
        # with that of the horizontal bar representing the sets,
        # hopefully to be get ridden of at some point
        self.ax_matrix.barh(\
            range(len(self.sets)), np.repeat(1,len(self.sets)),
            facecolor = "None")


        #--------------------------- Spines --------------------------#


        # Hide all spines
        self.ax_matrix.spines["top"].set_visible(False)
        self.ax_matrix.spines["right"].set_visible(False)       
        self.ax_matrix.spines["bottom"].set_visible(False)
        self.ax_matrix.spines["left"].set_visible(False)


        #--------------------------- x-axis --------------------------#


        # Do not plot ticks on the x axis
        self.ax_matrix.set_xticks([])


        #--------------------------- y-axis --------------------------#

        
        # Get the position of the ticks on the y-axis
        yticks = np.arange(len(self.sets))
        # Create the tick labels
        yticklabels = yticklabels if yticklabels else self.setnames
        # Set the ticks
        self.ax_matrix.set_yticks(yticks)
        # Hide the ticks
        self.ax_matrix.tick_params("y", length = 0)
        # Set the tick labels
        self.ax_matrix.set_yticklabels(yticklabels)


        #---------------------------- Font ---------------------------#


        # Create the font properties for the tick labels
        fp_ticklabels = fm.FontProperties(fname = custom_fpath,
                                          size = 7)
            
        # Set the font properties
        for yticklab in self.ax_matrix.get_yticklabels():
            yticklab.set_fontproperties(fp_ticklabels)


    def _draw_sizebars(self,
                       inter_color,
                       top,
                       yspace,
                       custom_fpath):
        """Draw bars representing the size of the intersections.
        """

        # On the x axis we have the numbers of the intersections
        x_data = np.arange(len(self.inters))
        
        # On the y axis we have the size of the intersections
        y_data = [inter[1] for inter in self.inters]

        # Set the bar width
        width = 0.5

        # Deal with the user-provided color representations
        is_rgbtuple = isinstance(inter_color, tuple)
        is_string = isinstance(inter_color, str)
        is_hex = False
        
        # If the string passed is a color HEX code
        if is_string and inter_color.startswith("#"):
            # Set the flag to identify HEX colors to True
            is_hex = True
        
        # If the color is either a RGB tuple or a HEX code
        if is_rgbtuple or is_hex:
            # Repeat it for the number of bars
            facecolor = [inter_color] * len(x_data)

        # If the color is a string (and not a HEX code)
        elif is_string and not is_hex:

            # If the string is among the recognized color
            # palettes   
            if inter_color in self.available_palettes:
                # Set the bars' colors from the color palette
                facecolor = sns.color_palette(\
                                inter_color, \
                                len(x_data))
            
            # If the string is a named color
            elif inter_color in self.available_colors:
                # Sets the bars' color to the named color
                facecolor = [inter_color] * len(x_data)

            # If the string is not recognized, raise an exception
            else:
                errstr = \
                    f"{inter_color} is neither a name of a color " \
                    f"palette nor a named color."
                raise ValueError(errstr)

        # Plot the bars
        barlist = self.ax_sizebars.bar(x = x_data,
                                       height = y_data,
                                       width = width)

        # Color the bars
        for b, c in zip(barlist, facecolor):
            b.set_color(c)

        # Set the upper limit of the y axis
        if top is None:
            top = max(y_data)

        # Set the lower limit of the y axis
        bottom = 0

        # Set the limits for the y axis
        self.ax_sizebars.set_ylim(bottom, top)

        # Set the tick spacing for the y axis
        if yspace is None:
            yspace = \
                self.ax_sizebars.get_yticks()[1] - \
                self.ax_sizebars.get_yticks()[0]

        # plot the ticks for the y ticks
        yticks = np.arange(bottom, top+yspace, yspace)
        self.ax_sizebars.set_yticks(yticks)         

        # Set lower and upper limit for the x axis
        left = min(x_data)
        right = max(x_data)

        # Set the limits for the x axis
        self.ax_sizebars.set_xlim(left-1, right+(width/2))

        # Do not plot ticks for the x axis
        self.ax_sizebars.set_xticks([])
       
        # Hide all spines apart fron the left one
        # (y axis)
        self.ax_sizebars.spines["top"].set_visible(False)
        self.ax_sizebars.spines["right"].set_visible(False)
        self.ax_sizebars.spines["bottom"].set_visible(False)

        # Create the font properties for the tick labels
        fp_ticklabels = fm.FontProperties(fname = custom_fpath,
                                          size = 7)
            
        # Set the font properties
        for yticklab in self.ax_sizebars.get_yticklabels():
            yticklab.set_fontproperties(fp_ticklabels)
        

    def _draw_setbars(self,
                      sets_color,
                      top,
                      yspace,
                      custom_fpath):
        """Draw the horizontal bars representing
        the size of the single sets.
        """
        
        # Set numbers will be on the x axis
        x_data = np.arange(len(self.sets))
        
        # Set sizes will be on the y axis
        y_data = [len(s) for s in self.sets]

        # Invert the x axis so that it goes
        # from right to left
        self.ax_setbars.invert_xaxis()

        # Do not draw ticks or tick labels for the x axis
        self.ax_setbars.set_yticks([])
        self.ax_setbars.set_yticklabels([])

        # Deal with the user-provided color
        # representation
        is_list = isinstance(sets_color, list)
        is_rgbtuple = isinstance(sets_color, tuple)
        is_string = isinstance(sets_color, str)
        is_hex = False

        # If the string passed is a color HEX code       
        if is_string and sets_color.startswith("#"):
            # Set the flag to identify HEX colors to True
            is_hex = True

        # If the color is either a RGB tuple or a HEX code 
        if is_rgbtuple or is_hex:
            # Repeat it for the number of sets
            facecolor = [sets_color] * len(x_data)

        # If the user passed a list of colors (either a list of
        # RGB tuples or a list of named colors)
        elif is_list:
            # Reverse it to match the original sets order
            sets_color.reverse()
            # The list of colors will be used as sets' colors
            facecolor = sets_color

        # If the color is a string (and not a HEX code)
        elif is_string and not is_hex:

            # If the string is among the recognized color
            # palettes
            if sets_color in self.available_palettes:
                # Set the bars' colors from the color palette
                facecolor = sns.color_palette(\
                                sets_color,
                                len(x_data))

            # If the string is a named color         
            elif sets_color in self.available_colors:
                # Sets the bars' color to the named color
                facecolor = [sets_color] * len(x_data)

            # If the string is not recognized, raise an exception
            else:
                errstr = \
                    f"{sets_color} is neither a name of a color " \
                    f"palette nor a named color."
                raise ValueError(errstr)           

        # Plot the bars
        barlist = \
            self.ax_setbars.barh(y = x_data,
                                 width = y_data,
                                 height = 0.8,
                                 align = "center")

        # Set the bar color(s)
        for b, c in zip(barlist,facecolor):
            b.set_color(c)
        
        # Hide all spines apart from the bottom one
        self.ax_setbars.spines["top"].set_visible(False)
        self.ax_setbars.spines["left"].set_visible(False)
        self.ax_setbars.spines["right"].set_visible(False)

        # Generate the font properties
        fp_ticklabels = fm.FontProperties(fname = custom_fpath,
                                          size = 7)
            
        # Set the properties for the tick labels on the x axis
        for xticklab in self.ax_setbars.get_xticklabels():
            xticklab.set_fontproperties(fp_ticklabels)

        # Set the properties for the tick labels on the y axis
        for yticklab in self.ax_setbars.get_yticklabels():
            yticklab.set_fontproperties(fp_ticklabels)


    def _draw_shading(self,
                      shading_color):
        """Shade one every two row of the matrix.
        """
        
        # Set the height of the shaded areas
        height = 0.8

        # For each row every two
        for i in range(0, len(self.sets), 2):
            
            # Create a shaded area
            rect = plt.Rectangle(xy = (-height/2, i-height/2),
                                 width = len(self.inters),
                                 height = height,
                                 facecolor = shading_color,
                                 lw = 0,
                                 zorder = 0)

            # Add the shaded areas to the plot
            self.ax_matrix.add_patch(rect)


    def plot(self,
             outfile = "upset.pdf",
             dpi = 900,
             transparent = True,
             bbox_inches = "tight",
             inter_color = "red",
             sets_color = "Set2",
             noninter_color = "#929292",
             shading_color = "#f1f1f1",
             linestyle_matrix = "-",
             yticklabels = None,
             top_sizebars = None,
             yspace_sizebars = None,
             top_setbars = None,
             yspace_setbars = None,
             custom_fpath = None):
        """Generate an UpSet plot.

        Parameters
        ----------
        outfile : `str`, default: `"upset.pdf"`
            Name of the file where to save
            the plot.

        dpi : `int`, default: 900
            DPI of the final figure.

        transparent : `bool`, default: True
            Draw the plot with a transparent
            background.

        inter_color : `str` or `tuple`
            Color(s) to use for points representing
            intersecting sets.
            If `str`, can be either a HEX code for
            a color, a named color (according to the 
            `matplotlib` convention) or the name of a 
            color palette (`matplotlib` or `seaborn`).
            If `tuple`, it is assumed to be a RGB tuple.

        sets_color : `str` or `tuple` or `list`
            Color(s) to use for bars representing
            sets.
            If `str`, can be either a HEX code for
            a color, a named color (according to the 
            `matplotlib` convention) or the name of a 
            color palette (`matplotlib` or `seaborn`).
            If `tuple`, it is assumed to be a RGB tuple.
            If `list`, it is assumed to be either a list
            of RGB tuples or a list of named colors.

        noninter_color: `str` or `tuple`
            Color to use for points representing
            non-intersecting sets.
            If `str`, it can be either a color HEX code,
            or a named color (according to the 
            `matplotlib` convention).
            If `tuple`, it is assumed to be a RGB tuple.

        shading_color : `str` or `tuple`
            Color of the shaded rows. 
            If `str`, it is assumed to be
            either a named color or a HEX
            code.
            If `tuple`, it is assumed to be
            a RGB tuple.

        linestyle_matrix : `str`, default: '-'
            String representing the line style to use
            for the vertical lines (according to
            the `matplotlib` convention).

        yticklabels : `list` or `None`, default: `None`
            List of labels to use for the ticks of
            the y axis. If `None`, set names will be
            used.

        top_sizebars : `int` or `None`, default: `None`
            Upper limit of the y axis of the size bars. 
            If `None`, it will be the highest interaction 
            size found.

        yspace_sizebars : `int` or `None`, default: `None`
            Spacing for the ticks on the y axis of the
            size bars. If `None`, it will be determined 
            automatically from the data.

        top_setbars : `int` or `None`, default: `None`
            Left (upper) limit of the y axis of the
            sets bars. If `None`, it will be the maximum 
            set size found in the data.

        yspace_setbars : `int` or `None`, default: `None`
            Spacing for the ticks on the y axis of the
            sets bars. If `None`, it will determined 
            automatically from the data.

        custom_fpath : `str` or `None`, default: `None`
            Path to a custom font to be used for
            all texts in the plot.

        Returns
        -------
        None
        """
        
        # Create the figure
        fig = plt.figure()

        # Define the grid for the plots 
        gs = gridspec.GridSpec(nrows = 2,
                               ncols = 3,
                               width_ratios = [0.2, 0.2, 1],
                               height_ratios = [1, 0.8])

        # Set the axes for each plot component
        self.ax_sizebars = plt.subplot(gs[2])
        self.ax_setbars = plt.subplot(gs[3])
        self.ax_matrix = plt.subplot(gs[5],
                                     sharex = self.ax_sizebars) 

        # Draw the bars with the sizes of the intersections
        self._draw_sizebars(inter_color = inter_color,
                            top = top_sizebars,
                            yspace = yspace_sizebars,
                            custom_fpath = custom_fpath)

        # Draw the shaded areas before the matrix
        # in order not to hide it
        self._draw_shading(shading_color = shading_color)
        
        # Draw the matrix
        self._draw_matrix(inter_color = inter_color,
                          noninter_color = noninter_color,
                          linestyle = linestyle_matrix,
                          custom_fpath = custom_fpath,
                          yticklabels = yticklabels)
        
        # Draw the bars with the sizes of the sets
        self._draw_setbars(sets_color = sets_color,
                           top = top_setbars,
                           yspace = yspace_setbars,
                           custom_fpath = custom_fpath)

        # Save the plot
        plt.savefig(outfile,
                    dpi = dpi,
                    transparent = transparent,
                    bbox_inches = bbox_inches)