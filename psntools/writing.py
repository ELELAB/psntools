#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    writing.py
#
#    Utilities to write out results from PSN analyses.
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


# Third-party packages
import pandas as pd
# psntools
import psntools.dataframes as dataframes



# Get the module logger
logger = log.getLogger(__name__)



############################### For PSNs ##############################



def write_psn_csv(psn,
                  outfile,
                  csv_sep = ",",
                  float_fmt = "%2.3f"):
    """Write a 2D `pandas.DataFrame` representation of the PSN
    to a CSV file. Rows and columns represent the nodes of the
    PSN, and each cell stores the value of the edge between two
    nodes. The dataframe is basically a labeled PSN matrix.

    Parameters
    ----------
    psn : `psntools.core.PSN`
        PSN.

    outfile : `str`
        Output file.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    Returns
    -------
    `None`
    """
        
    # Get the dataframe representing the PSN
    df = dataframes.get_psn_df(psn = psn)
    
    # Write the dataframe to the output file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt)


def write_nodes_list(psn, outfile):
    """Write a text file listing all nodes in the PSN as strings.

    Parameters
    ----------
    psn : `psntools.core.PSN`
        PSN.

    outfile : `str`
        Output file.

    Returns
    -------
    `None`
    """

    with open(outfile, "w") as out:
        
        # For all node attributes
        for attrs in psn.graph.nodes.data():
            
            # Write out the nodes as formatted strings
            out.write(psn.NODE_STR_FMT.format(**attrs))
            out.write("\n")


def write_nodes_csv(outfile,
                    df = None,
                    psn = None,
                    metrics = None,
                    csv_sep = ",",
                    float_fmt = "%2.3f"):
    """Write a CSV file containing a data frame with all nodes
    of a PSN. Rows of the data frame will be the PSN nodes, 
    identified by their string representation, and columns
    will contain the the values of node metrics of interest
    for each node.

    Parameters
    ----------
    outfile : `str`
        Output file.

    df : `pandas.DataFrame` or `None`, default: `None`
        Pre-computed data frame containing the values for the
        node metrics of interest (like the one returned by
        `psntools.dataframes.get_nodes_df`).

    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    metrics : `dict` or `None`, default: `None`
        Dictionary of node metrics to be computed, mapping
        to dictionaries with the keyword arguments to be passed
        to the functions computing such metrics.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    Returns
    -------
    `None`
    """

    # If no pre-computed data frame was passed
    if df is None:    
        
        # Build the nodes' data frame
        df = dataframes.get_nodes_df(psn = psn,
                                     metrics = metrics)

    # Write out the data frame as a CSV file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt,
              index = False)

    
def write_edges_csv(outfile,
                    df = None,
                    data = None,
                    psn = None,
                    sort_by = "node",
                    ascending = False,
                    csv_sep = ",",
                    float_fmt = "%2.3f",
                    **kwargs):
    """Write a CSV file listing the edges in the PSN.

    Parameters
    ----------
    outfile : `str`
        Output file.

    df : `pandas.DataFrame` or `None`, default: `None`
        Pre-computed data frame containing the edges
        (like the one returned by 
        `psntools.dataframes.get_edges_df`).

    data : `dict` or `None`, default: `None`
        Pre-computed dictionary containing the edges (like
        the one returned by the 
        `psntools.core.PSN.get_edges` method).

    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    sort_by : `str`, [`"node"`, `"weight"`],
               default: `"weight"`
        Whether to sort the edges by node name or weight.

    ascending : `bool`, default: `False`
        Whether the sorting is ascending. 

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    **kwargs:
        Keyword arguments to be passed to
        `psntools.core.PSN.get_edges`, if neither `df` nor `data`
        have been passed.

    Returns
    -------
    `None`
    """

    # If no pre-computed data frame was passed
    if df is None:

        # Build the data frame
        df = dataframes.get_edges_df(data = data,
                                     psn = psn,
                                     sort_by = sort_by,
                                     ascending = ascending,
                                     **kwargs)
        
    # Write the data frame to the output file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt,
              index = False)


def write_hubs_csv(outfile,
                   df = None,
                   psn = None,
                   sort_by = "degree",
                   ascending = False,
                   csv_sep = ",",
                   **kwargs):
    """Write a CSV file with the hubs found in the PSN.
        
    Parameters
    ----------
    outfile : `str`
        Output file.

    df : `pandas.DataFrame` or `None`, default: `None`
        Pre-computed data frame containing the hubs
        (like the one returned by 
        `psntools.dataframes.get_hubs_df`).

    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    sort_by : `str`, [`"degree"`, `"node"`], 
              default: `"degree"`
        Whether to sort the hubs by degree or by
        node name.

    sort_by : `str`, [`"node"`, `"weight"`],
               default: `"weight"`
        Whether to sort the hubs by node name or node degree.

    ascending : `bool`, default: `False`
        Whether the sorting is ascending.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    **kwargs:
        Keyword arguments to be passed to
        `psntools.core.PSN.get_hubs`, if neither `df` nor `data`
        have been passed.

    Returns
    -------
    `None`
    """

    # If no pre-computed data frame was passed
    if df is None:    
        
        # Build the hubs' data frame
        df = dataframes.get_hubs_df(psn = psn,
                                    sort_by = sort_by,
                                    ascending = ascending,
                                    **kwargs)
        
    # Write the data frame to the output file
    df.to_csv(outfile,
              sep = csv_sep,
              index = False)


def write_connected_components_csv(outfile,
                                   df = None,
                                   data = None,
                                   psn = None,
                                   cc_prefix = "CC_",
                                   node_sep = ".",
                                   csv_sep = ",",
                                   **kwargs):
    """Write a CSV file with the list of connected components found.

    Parameters
    ----------
    outfile : `str`
        Output file.

    df : `pandas.DataFrame` or `None`, default: `None`
        Pre-computed data frame containing the connected components
        (like the one returned by
        `psntools.dataframes.get_connected_components_df`).

    data : `list` or `None`, default: `None`
        List of sets representing the connected components
        (like the one returned by the 
        `psntools.core.PSN.get_connected_components` method).

    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    cc_prefix : `str`, default: `"CC_"`
        Prefix to add to each connected component's name (if it
        is an empty string, the components' names will be
        integers).

    node_sep : `str`, default: `"."`
        Separator for the nodes in each connected component.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    **kwargs
        Arguments to be passed to 
        `psntools.core.PSN.get_connected_components` for
        connected components calculation, if neither `df` nor
        `data` have been passed.

    Returns
    ------
    `None`
    """

    # If no pre-computed data frame was passed
    if df is None:

        # Build the data frame
        df = dataframes.get_connected_components_df(\
                data = data,
                cc_prefix = cc_prefix,
                node_sep = node_sep,
                **kwargs)
        
    # Write the data frame to the output CSV file
    df.to_csv(outfile,
              sep = csv_sep,
              index = False)


def write_shortest_paths_csvs(data = None,
                              psn = None,
                              outfiles_prefix = "path_",
                              sort_by = ("length", "weight"),
                              ascending = (False, False),
                              pair_node_sep = "_",
                              path_node_sep = ".",
                              csv_sep = ",",
                              float_fmt = "%2.3f",
                              **kwargs):
    """For each pair of nodes, write a CSV file containing
    the shortest paths found between them. CSV files will
    be named after the node pairs.
        
    Parameters
    ----------
    data : `dict` or `None`, default: `None`
        Shortest paths, if already calculated.

    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    outfiles_prefix : `str`, default: `"path_"`
        Prefix for the output files.

    sort_by : `tuple`, 
              [`("length", "weight")`, `("weight", "length")`], 
              default: `("length", "weight")`
        Whether to sort the paths primarily by length 
        and secondarily by weight or vice versa.

    ascending : `tuple`, default: `(False, False)`
        Whether the primary and secondary sorting
        will be in ascending order.

    pair_node_sep : `str`, default: `"_"`
        Separator for the nodes in each pair of nodes.

    path_node_sep : `str`, default: `"."`
        Separator for the nodes in each path.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV files.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV files.

    **kwargs
        Arguments to be passed to 
        `psntools.core.PSN.get_shortest_paths`for shortest
        paths calculation, if neither `df` nor `data`
        have been passed.

    Returns
    -------
    `None`
    """

    # If the shortest paths were not passed
    if data is None:
        # Set all the shortest paths
        data = psn.get_shortest_paths(node_fmt = "strings",
                                      **kwargs)
        
    # For each pair of nodes and associated shortest paths
    for pair, sps in data.items():
            
        # Set the file name
        name = f"{pair[0]}{pair_node_sep}{pair[1]}"
        
        # Add prefix and extension to the file name
        outfile = f"{outfiles_prefix}{name}.csv"
            
        # Create a list with the paths as strings,
        # the length and the weight of each path
        sps_str = [(path_node_sep.join(sp), len(sp), w)
                   for sp, w in sps.items()]
            
        # Generate the dataframe containing the shortest paths
        df = pd.DataFrame(sps_str)
            
        # Set the dataframe columns
        df.columns = ["path", "length", "weight"]
            
        # If sorting primarily by length and secondarily by weight
        if sort_by == ("length", "weight"):
            # Get the sorting columns
            sort_cols = [df.columns[1], df.columns[2]]
            
        # If sorting primarily by weight and secondarily by length
        elif sort_by == ("weight", "length"):
            # Get the sorting columns
            sort_cols = [df.columns[2], df.columns[1]]
            
        # Sort the paths
        df = df.sort_values(by = sort_cols,
                            ascending = ascending)
            
        # Write the dataframe to the output file
        df.to_csv(outfile,
                  sep = csv_sep,
                  float_format = float_fmt,
                  index = False)



############################ For PSN groups ###########################



def write_nodes_csv_psngroup(outfile,
                             df = None,
                             psngroup = None,
                             metric = None,
                             csv_sep = ",",
                             float_fmt = "%2.3f"):
    """Write a CSV file where, for a single node metric, values
    for all nodes of the PSNs in the group are reported. Rows of
    the dataframe represent nodes, while each column represents
    a single PSN.

    Parameters
    ----------
    outfile : `str`
        Output file.

    df : `pandas.DataFrame` or `None`, default: `None`
        Nodes' dataframe, if already created.

    psngroup : `psntools.core.PSNGroup` or `None`, default: `None`
        PSN group.

    metric : `dict`
        Label of a node metric to be computed, mapping to a
        dictionary with the keyword arguments to be passed
        to the function computing such metric.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    Returns
    -------
    `None`
    """

    # If the data frame was not passed
    if df is None:
        
        # Get the data frame
        df = dataframes.get_nodes_df_psngroup(psngroup = psngroup,
                                              metric = metric)

    # Write the CSV file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt)



def write_common_hubs_csvs(dfs = None,
                           data = None,
                           psngroup = None,
                           outfiles_prefix = "hubs_",
                           psn_sep = "_",
                           csv_sep = ",",
                           **kwargs):
    """Write the data frames contaning the common hubs for each
    possible subset of PSNs in a PSNgroup to a CSV file.

    Parameters
    ----------
    dfs : `dict` or `None`, default: `None`
        Dictionary of data frames of common hubs,
        if already calculated.

    data : `dict` or `None`, default: `None`
        Pre-computed dictionary of common hubs as obtained with
        `psntools.core.PSNGroup.get_common_hubs`.

    psngroup : `psntools.core.PSNGroup` or `None`, default: `None`
        PSN group.

    outfiles_prefix : `str`, default: `"hubs_"`
        Prefix for the output files.`

    psn_sep : `str`, default: `"_"`
        How to separate the PSN names in the names of the
        output files.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV files.

    **kwargs
        Arguments to be passed to
        `psntools.core.PSNGroup.get_common_hubs` for common
        hubs' calculation, if neither `df` nor `data`
        have been passed.

    Returns
    -------
    `None`
    """

    # If the dictionary of pre-computed data frames was
    # not passed
    if dfs is None:
        # Get the data frame of common hubs
        dfs = \
            dataframes.get_common_hubs_dfs(data = data,
                                           psngroup = psngroup,
                                           psn_sep = psn_sep,
                                           **kwargs)

    # For each name/data frame pair
    for name, df in dfs.items():
        # Join the PSN names
        name = psn_sep.join(name)
        # Add prefix and extension
        outfile = f"{outfiles_prefix}{name}.csv"
        # Write the data frame to the output file
        df.to_csv(outfile, sep = csv_sep)


def write_common_edges_csvs(dfs = None,
                            data = None,
                            psngroup = None,
                            outfiles_prefix = "edges_",
                            psn_sep = "_",
                            node_sep = "_",
                            csv_sep = ",",
                            float_fmt = "%2.3f",
                            **kwargs):
    """Write the data frames of common edges for each possible
    subset of PSNs in a PSNgroup to a CSV file.

    Parameters
    ----------
    dfs : `pandas.DataFrame` or `None`, default: `None`
        Dictionary of data frames of common edges,
        if already calculated.

    data : `dict` or `None`, default: `None`
        Pre-computed dictionary of common edges as obtained with
        `core.PSNGroup.get_common_edges`.

    psngroup : `psntools.core.PSNGroup` or `None`, default: `None`
        PSN group.

    outfiles_prefix : `str`, default: `"edges_"`
        Prefix for the output files.

    psn_sep : `str`, default: `"_"`
        How to separate the PSN names in the names of the
        output files.

    node_sep : `str`, default: `"_"`
        How to separate nodes in each edge.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV files.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV files.

    **kwargs
        Arguments to be passed to
        `psntools.core.PSNGroup.get_common_edges` for common
        edges' calculation, if neither `df` nor `data`
        have been passed.

    Returns
    -------
    `None`
    """

    # If the dictionary of pre-computed data frames was
    # not passed
    if dfs is None:
        
        # Get the data frame of common edges
        dfs = dataframes.get_common_edges_dfs(data = data,
                                              psngroup = psngroup,
                                              psn_sep = psn_sep,
                                              node_sep = node_sep,
                                              **kwargs)

    # For each name/data frame pair
    for name, df in dfs.items():
        
        # Join the PSN names
        name = psn_sep.join(name)
        
        # Add prefix and extension
        outfile = f"{outfiles_prefix}{name}.csv"
        
        # Write the data frame to the output file
        df.to_csv(outfile,
                  sep = csv_sep,
                  float_format = float_fmt)
