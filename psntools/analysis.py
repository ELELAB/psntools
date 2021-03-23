#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    analysis.py
#
#    Utilities to extract and manipulate data in PSNs.
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
import pandas as pd
import networkx as nx



################################ For PSN ##############################



def get_psn_df(psn, \
               node_fmt = "residues"):
    """Get a 2D `pandas.DataFrame` representation of the PSN.
    Rows and columns represent the nodes of the PSN, and 
    each cell stores the value of the edge between two
    nodes. The dataframe is basically a labeled PSN matrix.

    Parameters
    ----------
    psn : `psntools.PSN.PSN`
        PSN.

    node_fmt : `str` [`"residues"`, `"strings"`],
          default: `"residues"`
        Whether to represent nodes as 
        `MDAnalysis.core.groups.Residue` instances
        (`"residues"`) or formatted strings (`"strings"`).

    Returns
    -------
    df : `pandas.DataFrame`
        Output dataframe.
    """
        
    # convert the graph to a Pandas adjacency matrix (= dataframe)
    df = nx.convert_matrix.to_pandas_adjacency(G = psn.graph)
        
    # if nodes should be represented as Residue instances
    if node_fmt == "residues":
        # simply return the dataframe
        return df
        
    # if nodes should be represented as strings
    elif node_fmt == "strings":
        # generate the labels for the PSN
        df_labels = psn.get_nodes_residues2strings()
        # set index and columns of the dataframe 
        df = df.rename(mapper = df_labels, axis = 0)
        df = df.rename(mapper = df_labels, axis = 1)
        # return the dataframe
        return df


def get_nodes_df(data = None, \
                 psn = None, \
                 metrics = None):
    """Create a dataframe with nodes of a PSN and node metrics,
    if requested. The node metrics will be computed on the fly or a
    dictionary of dictionaries keyed on the metrics' names and on the
    nodes (they must be `MDAnalysis.core.groups.Residue` instances)
    can be passed if they were already calculated elsewhere.

    Parameters
    ----------
    data : `dict` or `None`, default: `None`
        Dictionary of dictionaries mapping the metrics names to
        internal dictionaries where keys are nodes and values are
        metrics values for the nodes. Nodes must be
        `MDAnalysis.core.groups.Residue` instances.

    psn : `psntools.PSN.PSN` or `None`, default: `None`
        PSN.

    metrics : `dict` or `None`, default: `None`
        Dictionary of node metrics to be computed, mapping
        to dictionaries with the keyword arguments to be passed
        to the functions computing such metrics.
        
    Returns
    -------
    df : `pandas.DataFrame`
        Output dataframe.
    """
        
    # get the node attributes
    attrs = psn.get_nodes_attributes()
        
    # build a dataframe from such attributes
    df = pd.DataFrame(attrs)
    # set the index column as the column storing
    # the residues' indexes
    df = df.set_index("ix", drop = True)
    # insert a column with the name of the residue as
    # a formatted string
    df["formatted_string"] = \
        pd.Series({n["ix"] : psn.NODE_STR_FMT.format(**n) \
                  for n in attrs})      

    # if no data were passed
    if not data:     
        # convert it to a dictionary where each metric
        # maps to an empty dictionary
        data = {m : {} for m in metrics}

    # for each metric, node values in the data dictionary
    for m, nv in data.items():      
        # if there are not values associated to the metric
        # (i.e. the metric needs to be computed)
        if not nv:
            # get the keyword arguments to be passed to the
            # function calculating the metric
            metric_kws = {**metrics[m], **{"node_fmt" : "residues"}}
            # calculate the metric for all nodes
            nv = psn.get_metric(metric = m, \
                                kind = "node", \
                                metric_kws = metric_kws)
                
        # add it as a column to the dataframe
        df[m] = pd.Series({n.ix : v for n, v in nv.items()})

    # set the index column as the column storing
    # the residues as formatted strings
    df = df.set_index("formatted_string", drop = True)

    # return the dataframe
    return df



############################# For PSNGroup ############################



def get_nodes_df_psngroup(psngroup, \
                          metric = None):
    """Get a dataframe where, for a single node metric, values
    for all nodes of the PSNs in the group are reported. Rows of
    the dataframe represent nodes, while each column represents
    a single PSN.

    Parameters
    ----------
    psngroup : `psntools.PSN.PSNGroup`
        PSN group.

    metric : `dict`
        Label of a node metric to be computed, mapping to a
        dictionary with the keyword arguments to be passed
        to the function computing such metric.

    Returns
    -------
    df : `pandas.DataFrame`
        Dataframe contaning values for that metric for all
        nodes in all PSNs.
    """

    # create a dictionary to store the value of the metric for
    # all PSNs in the group
    dfs = {}

    # create a variable to store the node names in order,
    # because Pandas reindex the dataframe when merging
    df_index = []

    # for each PSN
    for label, psn in psngroup.psns.items():
        # get the nodes' dataframe for the single PSN
        node_df = get_nodes_df(psn = psn, \
                               metrics = metric)
        # add the information about chain ID and residue
        # number for each node of the current PSN to the
        # dataframe that will be used for indexing
        df_index.append(node_df[["segid", "resnum"]])
        # add the metric values for the nodes of the current
        # PSN to the dictionary that store all values for
        # all PSNs in the group
        metric_name = list(metric.keys())[0]
        dfs[label] = node_df[metric_name].squeeze()
        
    # create the index to be used in the final dataframe
    index = pd.concat(df_index).drop_duplicates().index

    # create the dataframe and reindex it
    df = pd.DataFrame(dfs).reindex(index)

    # return the dataframe
    return df


def get_connected_components_df_psngroup(psngroup, \
                                         n_ccs = 5, \
                                         cc_prefix = "CC_"):
    """Get a dataframe where, for a single node metric, values
    for all nodes of the PSNs in the group are reported. Rows of
    the dataframe represent nodes, while each column represents
    a single PSN.

    Parameters
    ----------
    psngroup : `psntools.PSN.PSNGroup`
        PSN group.

    n_ccs : `int` or `None`, default: `5`
        How many of the most populated connected components
        should be plotted (components will be sorted by size
        starting from the biggest one before plotting).

    cc_prefix : `str`, default: `"CC_"`
        Prefix to add to each connected component's name (if it
        is an empty string, the components' names will be
        integers).

    Returns
    -------
    df : `pandas.DataFrame`
        Dataframe contaning the size of the most populated
        connected components in the PSNs of the group.
    """

    # initialize an empty dictionary to store the data
    data = {}
    
    # for each PSN in the group
    for psn_label, psn in psngroup.psns.items():

        # get connected components (already sorted by size with the
        # biggest ones first)
        ccs = sorted(psn.get_connected_components(), \
                     key = lambda x : len(x), \
                     reverse = True)[:n_ccs]

        # add the data about the size of the most populated connected
        # components in the current PSN
        data[psn_label] = \
            {f"{cc_prefix}{i+1}" : len(cc) for i, cc in enumerate(ccs)}

    # convert the data into a dataframe
    df = pd.DataFrame(data).T

    # return the dataframe
    return df


