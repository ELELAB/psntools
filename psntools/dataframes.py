#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dataframes.py
#
#    Utilities to create data frames containing data from analyses
#    of PSNs or PSN groups.
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
import logging as log
# Third-party packages
import pandas as pd
import networkx as nx
# psntools
from psntools import core



# Get the module logger
logger = log.getLogger(__name__)



################################ For PSN ##############################



def get_psn_df(psn):
    """Get a 2D `pandas.DataFrame` representation of the PSN.
    Rows and columns represent the nodes of the PSN, and 
    each cell stores the value of the edge between two
    nodes. The data frame is basically a labeled PSN matrix.

    Parameters
    ----------
    psn : `psntools.core.PSN`
        PSN.

    Returns
    -------
    df : `pandas.DataFrame`
        Output data frame.
    """
        
    # Convert the graph to a Pandas adjacency matrix (= data frame)
    df = nx.convert_matrix.to_pandas_adjacency(G = psn.graph)
        
    # Generate the labels for the PSN
    df_labels = psn.get_nodes_residues2strings()
    
    # Set index and columns of the data frame 
    df = df.rename(mapper = df_labels, axis = 0)
    df = df.rename(mapper = df_labels, axis = 1)
    
    # Return the data frame
    return df


def get_nodes_df(psn,
                 metrics = None):
    """Create a data frame with nodes of a PSN and node metrics,
    if requested. The node metrics will be computed on the fly or a
    dictionary of dictionaries keyed on the metrics' names and on the
    nodes (they must be `MDAnalysis.core.groups.Residue` instances)
    can be passed if they were already calculated elsewhere.

    Parameters
    ----------
    psn : `psntools.core.PSN`
        PSN.

    metrics : `dict` or `None`, default: `None`
        Dictionary of node metrics to be computed, mapping
        to dictionaries with the keyword arguments to be passed
        to the functions computing such metrics.
        
    Returns
    -------
    df : `pandas.DataFrame`
        Output data frame.
    """
        
    # Get the node attributes
    nodes, attrs = zip(*psn.graph.nodes.data())
        
    # Build a dataframe from such attributes
    df = pd.DataFrame(attrs)
    
    # Set the index column as the column storing
    # the residues' indexes
    df = df.set_index("ix", drop = True)
    
    # Insert a column with the name of the residue as
    # a formatted string
    df["formatted_string"] = \
        pd.Series({n["ix"] : psn.NODE_STR_FMT.format(**n) \
                  for n in attrs})

    # For each metric, node values in the data dictionary
    for m, m_kws in metrics.items():      
        
        # Get the keyword arguments to be passed to the
        # function calculating the metric
        metric_kws = {**metrics[m], **{"node_fmt" : "residues"}}
        
        # Calculate the metric for all nodes
        nv = psn.get_metric(metric = m,
                            kind = "node",
                            metric_kws = metric_kws)
                
        # Add it as a column to the data frame
        df[m] = pd.Series({n.ix : v for n, v in nv.items()})

    # Set the index column as the column storing
    # the residues as formatted strings
    df = df.set_index("formatted_string", drop = True)

    # Return the data frame
    return df


def get_edges_df(data = None,
                 psn = None,
                 sort_by = "node",
                 ascending = False,
                 **kwargs):
    """Create a data frame listing the edges in the PSN.

    Parameters
    ----------
    data : `dict` or `None`, default: `None`
        Pre-computed dictionary containing the edges (like
        the one returned by the `core.PSN.get_edges` method).

    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    sort_by : `str`, [`"node"`, `"weight"`],
               default: `"weight"`
        Whether to sort the edges by node name or weight.

    ascending : `bool`, default: `False`
        Whether the sorting is ascending. 

    **kwargs:
        Keyword arguments to be passed to 
        `psntools.core.PSN.get_edges`, if `data` has not
        been passed.

    Returns
    -------
    `pandas.DataFrame`
        Output data frame.
    """

    # If no pre-computed data were passed
    if data is None:
        # Get the edges
        data = psn.get_edges(node_fmt = "strings", **kwargs)
        
    # Convert the dictionary to a list of flat tuples
    edges = [(n1, n2, v) for (n1, n2), v in data.items()]
    
    # Generate a data frame from the list
    df = pd.DataFrame(edges,
                      columns = ["node1", "node2", "weight"])

    # If sorting by node was requested
    if sort_by == "node":
        
        # Sorting columns are those containing the nodes
        sort_cols = [df.columns[0], df.columns[1]]
        
        # The user's preference of sorting order will
        # be applied to both columns
        ascending = [ascending] * 2
        
    # If sorting by weight was requested
    elif sort_by == "weight":
        
        # Sorting columns are weight and then those
        # containing the nodes
        sort_cols = [df.columns[2], df.columns[0], df.columns[1]]
        
        # The user's preference of sorting will
        # only be applied to the weight column, while
        # the secondary sorting on nodes will happen
        # in ascending lexicographic order
        ascending = [ascending, True, True]

    # Sort the edges
    df = df.sort_values(by = sort_cols,
                        ascending = ascending)

    # Return the data frame
    return df


def get_hubs_df(psn = None,
                sort_by = "degree",
                ascending = False,
                **kwargs):
    """Create a data frame listing the hubs found in the PSN.

    Parameters
    ----------
    psn : `psntools.core.PSN` or `None`, default: `None`
        PSN.

    sort_by : `str`, [`"node"`, `"weight"`],
               default: `"weight"`
        Whether to sort the hubs by node name or node degree.

    ascending : `bool`, default: `False`
        Whether the sorting is ascending. 

    **kwargs:
        Keyword arguments to be passed to
        `psntools.core.PSN.get_hubs`, if `data` has not
        been passed.

    Returns
    -------
    `pandas.DataFrame`
        Output data frame.

    Notes
    -----
    `psn` is always needed since we need to get the
    mapping between MDAnalysis Residue instances and
    the corresponding formatted strings, therefore there
    is no option available to pass directly a dictionary
    of pre-computed hubs.
    """

    # Get the hubs
    data = psn.get_hubs(node_fmt = "residues", **kwargs)

    # Get the mapping between Residue instances
    # and formatted strings
    residues2strings = psn.get_nodes_residues2strings()
        
    # Create a list containing hubs' segment IDs,
    # residue numbers, formatted strings and degrees
    hubs = [(h.segid, h.resnum, residues2strings[h], d)
             for h, d in data.items()]
        
    # Generate a dataframe from the list
    df = pd.DataFrame(hubs)
        
    # If sorting by node name was requested
    if sort_by == "node":
        
        # Sorting columns will be segment ID and
        # residue number
        sort_cols = [df.columns[0], df.columns[1]]
        
        # Both will respect the user's decision
        # about the sorting order
        ascending = [ascending] * 2
        
    # If sorting by degree was requested
    elif sort_by == "degree":
        
        # Sorting columns will be primarily degree,
        # secondarily segment ID and residue number
        sort_cols = [df.columns[3], *df.columns[:2]]
        
        # Only sorting by degree will respect the
        # user's decision about sorting order, while
        # secondary sorting will always be in
        # ascending order
        ascending = [ascending, True, True]     
        
    # Sort the hubs
    df = df.sort_values(by = sort_cols,
                        ascending = ascending)
        
    # Drop the first two columns (only used for sorting)
    df = df.drop(df.columns[:2], axis = 1)
        
    # Set the data frame columns
    df.columns = ["node", "value"]

    # Return the data frame
    return df


def get_connected_components_df(data = None, 
                                psn = None, 
                                cc_prefix = "CC_",
                                node_sep = "."):
    """Get a data frame containing the connected components
    found in a PSN.

    Parameters
    ----------
    data : `list` or `None`, default: `None`
        List of sets representing the connected components
        (as the one returned by the 
        `psntools.core.PSN.get_connected_components` method).

    psn : `core.PSN.PSN` or `None`, default: `None`
        PSN.

    cc_prefix : `str`, default: `"CC_"`
        Prefix to add to each connected component's name (if it
        is an empty string, the components' names will be
        integers).

    node_sep : `str`, default: `"."`
        Separator for the nodes in each connected component.
        
    Returns
    -------
    df : `pandas.DataFrame`
        Output data frame.
    """

    # If no pre-computed data were passed
    if data is None:
        
        # Compute the connected components
        data = psn.get_connected_components(node_fmt = "strings",
                                            **kwargs)

    # Convert each connected component to a string. Start
    # numbering the connected components from 1.
    data = {f"{cc_prefix}{i}" : node_sep.join(cc)
           for i, cc in enumerate(data, 1)}

    # Generate a data frame
    df = pd.DataFrame.from_dict(ccs, orient = "index").reset_index()
        
    # Set the data frame columns
    df.columns = ["cc", "nodes"]

    # Return the data frame
    return df



############################# For PSNGroup ############################



def get_nodes_df_psngroup(psngroup,
                          metric):
    """Get a dataframe where, for a single node metric, values
    for all nodes of the PSNs in the group are reported. Rows of
    the data frame represent nodes, while each column represents
    a single PSN.

    Parameters
    ----------
    psngroup : `psntools.core.PSNGroup`
        PSN group.

    metric : `dict`
        Label of a node metric to be computed, mapping to a
        dictionary with the keyword arguments to be passed
        to the function computing such metric.

    Returns
    -------
    df : `pandas.DataFrame`
        Data frame contaning values for the chosen metric for all
        nodes in all PSNs.
    """

    # Create a dictionary to store the value of the metric for
    # all PSNs in the group
    dfs = {}

    # Create a variable to store the node names in order,
    # because Pandas reindex the dataframe when merging
    df_index = []

    # For each PSN
    for label, psn in psngroup.psns.items():
        
        # Get the nodes' dataframe for the single PSN
        node_df = get_nodes_df(psn = psn,
                               metrics = metric)
        
        # Add the information about chain ID and residue
        # number for each node of the current PSN to the
        # dataframe that will be used for indexing
        df_index.append(node_df[["segid", "resnum"]])
        
        # Add the metric values for the nodes of the current
        # PSN to the dictionary that store all values for
        # all PSNs in the group
        metric_name = list(metric.keys())[0]
        dfs[label] = node_df[metric_name].squeeze()

    # Create the data frame and reindex it
    df = pd.DataFrame(dfs)

    # Rename the index column
    df.index.name = "residues"

    # Return the data frame
    return df


def get_largest_connected_components_df_psngroup(psngroup,
                                                 n_ccs = 5,
                                                 cc_prefix = "CC_"):
    """Get a data frame where the n-th most populated
    connected components for each PSN in the group are
    reported.

    Parameters
    ----------
    psngroup : `psntools.core.PSNGroup`
        PSN group.

    n_ccs : `int` or `None`, default: `5`
        How many of the most populated connected components
        should be reported (components will be sorted by size
        starting from the biggest one).

    cc_prefix : `str`, default: `"CC_"`
        Prefix to add to each connected component's name (if it
        is an empty string, the components' names will be
        integers, starting from 0).

    Returns
    -------
    df : `pandas.DataFrame`
        Data frame contaning the size of the most populated
        connected components in the PSNs of the group.
    """

    # Initialize an empty dictionary to store the data
    data = {}
    
    # For each PSN in the group
    for psn_label, psn in psngroup.psns.items():

        # Get connected components (already sorted by size with the
        # biggest ones first)
        ccs = sorted(psn.get_connected_components(),
                     key = lambda x : len(x),
                     reverse = True)[:n_ccs]

        # Add the data about the size of the most populated connected
        # components in the current PSN
        data[psn_label] = \
            {f"{cc_prefix}{i+1}" : len(cc) for i, cc in enumerate(ccs)}

    # Convert the data into a dataframe
    df = pd.DataFrame(data).T

    # Return the dataframe
    return df


def get_common_hubs_dfs(data = None,
                        psngroup = None,
                        psn_sep = "_",
                        **kwargs):
    """Get the data frames contaning the common hubs for each
    possible subset of PSNs in a PSNgroup.

    Parameters
    ----------
    data : `dict` or `None`, default: `None`
        Pre-computed dictionary of common hubs as obtained with
        `psntools.core.PSNGroup.get_common_hubs`.

    psngroup : `psntools.PSN.PSNGroup` or `None`, default: `None`
        PSN group. Should be passed if `data` is not passed.

    psn_sep : `str`, default: `"_"`
        How to separate the PSN names in the names of the
        output files.

    **kwargs
        Arguments to be passed to
        `psntools.core.PSNGroup.get_common_hubs`
        for common hubs calculation, if `data`
        has not been passed.

    Returns
    -------
    `None`
    """

    # Create a dictionary to store the output data frames
    dfs = {}

    # If the common hubs were not passed
    if data is None:
        
        # Get the common hubs
        data = psngroup.get_common_hubs(**kwargs)

    # For each combination of PSNs
    for combo_label, combo in data.items():
        
        # Get the specific name of the file
        name = psn_sep.join(combo_label)
        
        # Generate the data frame
        df = pd.DataFrame.from_dict(combo, orient = "index")
        
        # Store the data frame under the combo label in the
        # output dictionary
        dfs[combo_label] = df

    # Return the dictionary of data frames
    return dfs


def get_common_edges_dfs(data = None,
                         psngroup = None,
                         psn_sep = "_",
                         node_sep = "_",
                         **kwargs):
    """Get data frames containing the common edges for each 
    possible subset of PSNs in a PSNGroup.

    Parameters
    ----------
    data : `dict` or `None`, default: `None`
        Pre-computed dictionary of common edges as obtained with
        `psntools.core.PSNGroup.get_common_edges`.

    psngroup : `psntools.PSN.PSNGroup` or `None`, default: `None`
        PSN group. Should be passed if `data` is not passed.

    psn_sep : `str`, default: `"_"`
        How to separate the PSN names in the names of the
        output files.

    node_sep : `str`, default: `"_"`
        How to separate nodes in each edge.

    **kwargs
        Arguments to be passed to
        `psntools.core.PSNGroup.get_common_edges`
        for common edges calculation, if `data`
        has not been passed.

    Returns
    -------
    `None`
    """

    # Create a dictionary to store the output data frames
    dfs = {}

    # If the common edges were not passed
    if data is None:
        
        # Get the common edges
        data = psngroup.get_common_edges(**kwargs)

    # For each PSN combination
    for combo_label, combo in data.items():
            
        # Replace tuples representing edges with strings
        combo = {p : {node_sep.join(e) : w for e,w in values.items()}
                 for p, values in combo.items()}
            
        # Get the specific name of the file
        name = psn_sep.join(combo_label)   
            
        # Generate the data frame
        df = pd.DataFrame.from_dict(combo, orient = "index") 

        # Store the data frame under the combo label in the
        # output dictionary
        dfs[combo_label] = df

    # Return the dictionary of data frames
    return dfs
