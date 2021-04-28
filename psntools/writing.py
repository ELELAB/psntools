#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    writing.py
#
#    Utilities to write out results from PSN analyses.
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
# psntools
from . import analysis


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
    psn : `psntools.PSN.PSN`
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
        
    # get the dataframe representing the PSN
    df = analysis.get_psn_df(psn = psn, \
                             node_fmt = "strings")
    # write the dataframe to the output file
    df.to_csv(outfile, \
              sep = csv_sep, \
              float_format = float_fmt)


def write_nodes_list(psn, outfile):
    """Write a text file listing all nodes in the PSN as strings.

    Parameters
    ----------
    psn : `psntools.PSN.PSN`
        PSN.

    outfile : `str`
        Output file.

    Returns
    -------
    `None`
    """

    with open(outfile, "w") as out:
        # for all node attributes
        for attrs in psn.get_nodes_attributes():
            # write out the nodes as formatted strings
            out.write(psn.NODE_STR_FMT.format(**attrs))
            out.write("\n")


def write_nodes_csv(psn,
                    outfile,
                    metrics = None,
                    csv_sep = ",",
                    float_fmt = "%2.3f"):
    """Write a CSV file with all nodes of a PSN. Rows of the
    dataframe will be the PSN nodes, identified by their string
    representation. Node metrics can be added to the dataframe,
    and can be computed on the fly or passed as a dictionary
    of dictionaries keyed on the metrics' names and on the nodes
    (nodes must be `MDAnalysis.core.groups.Residue` instances
    in such dictionary).

    Parameters
    ----------
    psn : `psntools.PSN.PSN`
        PSN.

    outfile : `str`
        Output file.
        
    metrics : `list` or `dict` or `None`, default: `None`
        List of node-wise metrics to be included in the CSV file or
        dictionary of dictionaries mapping the metrics names to
        internal dictionaries where keys are nodes and values are
        metrics values for the nodes. Nodes must be
        `MDAnalysis.core.groups.Residue` instances.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    Returns
    -------
    `None`
    """

    # get the nodes' dataframe
    df = get_nodes_df(psn = psn,
                      metrics = metrics)

    # write out the dataframe as a CSV file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt,
              index = False)

    
def write_edges_csv(psn,
                    outfile,
                    metrics = None,
                    sort_by = "node",
                    ascending = False,
                    csv_sep = ",",
                    float_fmt = "%2.3f",
                    **kwargs):
    """Write a CSV file listing edges in the PSN.

    Parameters
    ----------
    psn : `psntools.PSN.PSN`
        PSN.

    outfile : `str`
        Output file.

    metrics : `list` or `dict` or `None`, default: `None`
        List of edge-wise metrics to be included in the CSV file or
        dictionary of dictionaries mapping the metrics names to
        internal dictionaries where keys are edges and values are
        metrics values for the nodes. Edges must have nodes as
        `MDAnalysis.core.groups.Residue` instances.

    sort_by : `str`, [`"node"`, `"weight"`],
               default: `"weight"`
        Whether to sort the edges by node name or weight.

    ascending : `bool`, default: `False`
        Whether the sorting is ascending. 

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    Returns
    -------
    `None`
    """

    # get the edges of the PSN
    edges = psn.get_edges(node_fmt = "strings", **kwargs)
        
    # convert the dictionary to a list of flat tuples
    edges = [(*k, v) for k, v in edges.items()]
    
    # generate a dataframe from the list
    df = pd.DataFrame(edges,
                      columns = ["node1", "node2", "weight"])
        
    # if sorting by node was requested
    if sort_by == "node":
        # sorting columns are those containing the nodes
        sort_cols = [df.columns[0], df.columns[1]]
        # the user's preference of sorting order will
        # be applied to both columns
        ascending = [ascending] * 2
        
    # if sorting by weight was requested
    elif sort_by == "weight":
        # sorting columns are weight and then those
        # containing the nodes
        sort_cols = [df.columns[2], df.columns[0], df.columns[1]]
        # the user's preference of sorting will
        # only be applied to the weight column, while
        # the secondary sorting on nodes will happen
        # in ascending lexicographic order
        ascending = [ascending, True, True]

    # sort the edges
    df = df.sort_values(by = sort_cols,
                        ascending = ascending)
        
    # write the dataframe to the output file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt,
              index = False)


def write_hubs_csv(outfile,
                   hubs = None,
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

    hubs : `dict` or `None`, default: `None`
        Hubs, if already calculated.

    psn : `psntools.PSN.PSN` or `None`, default: `None`
        PSN.

    sort_by : `str`, [`"degree"`, `"node"`], 
              default: `"degree"`
        Whether to sort the hubs by degree or by
        node name.

    ascending : `bool`, default: `False`
        Whether to sort the hubs in ascending
        or descending order.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    **kwargs
        Arguments to be passed to `psn.get_hubs` for hubs
        calculation, if `hubs` has not been passed.

    Returns
    -------
    `None`
    """

    # if the hubs were not passed
    if hubs is None:
        # get the hubs
        hubs = psn.get_hubs(node_fmt = "residues", **kwargs)
        
    # get the mapping between Residue instances
    # and formatted strings
    residues2strings = psn.get_nodes_residues2strings()
        
    # create a list containing hubs' segment IDs,
    # residue numbers, formatted strings and degrees
    hubs = [(h.segid, h.resnum, residues2strings[h], d)
             for h, d in hubs.items()]
        
    # generate a dataframe from the list
    df = pd.DataFrame(hubs)
        
    # if sorting by node name was requested
    if sort_by == "node":
        # sorting columns will be segment ID and
        # residue number
        sort_cols = [df.columns[0], df.columns[1]]
        # both will respect the user's decision
        # about the sorting order
        ascending = [ascending] * 2
        
    # if sorting by degree was requested
    elif sort_by == "degree":
        # sorting columns will be primarily degree,
        # secondarily segment ID and residue number
        sort_cols = [df.columns[3], *df.columns[:2]]
        # only sorting by degree will respect the
        # user's decision about sorting order, while
        # secondary sorting will always be in
        # ascending order
        ascending = [ascending, True, True]     
        
    # sort hubs
    df = df.sort_values(by = sort_cols,
                        ascending = ascending)
        
    # drop the first two columns (only used for sorting)
    df = df.drop(df.columns[:2], axis = 1)
        
    # set the dataframe columns
    df.columns = ["node", "value"]
        
    # write the dataframe to the output file
    df.to_csv(outfile,
              sep = csv_sep,
              index = False)


def write_connected_components_csv(outfile,
                                   connected_components = None,
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

    connected_components : `dict` or `None`, default: `None`
        Connected components, if already calculated.

    psn : `psntools.PSN.PSN` or `None`, default: `None`
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
        Arguments to be passed to `psn.get_connected_components` for
        connected components calculation, if `connected_components`
        has not been passed.

    Returns
    ------
    `None`
    """

    # if the connected components were not passed
    if connected_components is None:
        # get the connected components
        ccs = psn.get_connected_components(node_fmt = "strings",
                                           **kwargs)
        
    # convert each connected component to a string. Start
    # numbering the connected components from 1.
    ccs = {f"{cc_prefix}{i}" : node_sep.join(cc)
           for i, cc in enumerate(ccs, 1)}
        
    # generate a dataframe
    df = pd.DataFrame.from_dict(ccs, orient = "index").reset_index()
        
    # set the dataframe columns
    df.columns = ["cc", "nodes"]
        
    # write the dataframe to the output CSV file
    df.to_csv(outfile,
              sep = csv_sep,
              index = False)


def write_shortest_paths_csvs(shortest_paths = None,
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
    shortest_paths : `dict` or `None`, default: `None`
        Shortest paths, if already calculated.

    psn : `psntools.PSN.PSN` or `None`, default: `None`
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
        Arguments to be passed to `psn.get_shortest_paths`
        for shortest paths calculation, if `shortest_paths`
        has not been passed.

    Returns
    -------
    `None`
    """

    # if the shortest paths were not passed
    if shortest_paths is None:
        # get all the shortest paths
        shortest_paths = psn.get_shortest_paths(node_fmt = "strings",
                                                **kwargs)
        
    # for each pair of nodes and associated shortest paths
    for pair, sps in shortest_paths.items():
            
        # set the file name
        name = f"{pair[0]}{pair_node_sep}{pair[1]}"
        # add prefix and extension to the file name
        outfile = f"{outfiles_prefix}{name}.csv"
            
        # create a list with the paths as strings,
        # the length and the weight of each path
        sps_str = [(path_node_sep.join(sp), len(sp), w)
                   for sp, w in sps.items()]
            
        # generate the dataframe containing the shortest paths
        df = pd.DataFrame(sps_str)
            
        # set the dataframe columns
        df.columns = ["path", "length", "weight"]
            
        # if sorting primarily by length and secondarily by weight
        if sort_by == ("length", "weight"):
            # get the sorting columns
            sort_cols = [df.columns[1], df.columns[2]]
            
        # if sorting primarily by weight and secondarily by length
        elif sort_by == ("weight", "length"):
            # get the sorting columns
            sort_cols = [df.columns[2], df.columns[1]]
            
        # sort the paths
        df = df.sort_values(by = sort_cols,
                            ascending = ascending)
            
        # write the dataframe to the output file
        df.to_csv(outfile,
                  sep = csv_sep,
                  float_format = float_fmt,
                  index = False)



############################ For PSN groups ###########################



def write_common_hubs_csvs(common_hubs = None,
                           psngroup = None,
                           outfiles_prefix = "hubs_",
                           psn_sep = "_",
                           csv_sep = ",",
                           **kwargs):
    """Write the common hubs for each possible subset of PSNs
    in the PSNgroup to a CSV file.

    Parameters
    ----------
    common_hubs : `dict` or `None`, default: `None`
        Common hubs, if already calculated.

    psngroup : `psntools.PSN.PSNGroup` or `None`, default: `None`
        PSN group.

    outfiles_prefix : `str`, default: `"hubs_"`
        Prefix for the output files.`

    psn_sep : `str`, default: `"_"`
        How to separate the PSN names in the names of the
        output files.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV files.

    **kwargs
        Arguments to be passed to `psngroup.get_common_hubs`
        for common hubs calculation, if `common_hubs`
        has not been passed.

    Returns
    -------
    `None`
    """

    # if the common hubs were not passed
    if common_hubs is None:
        # get the common hubs
        common_hubs = psngroup.get_common_hubs(**kwargs)

    # for each combination of PSNs
    for combo_label, combo in common_hubs.items():
        # get the specific name of the file
        name = psn_sep.join(combo_label)
        # add prefix and extension
        outfile = f"{outfiles_prefix}{name}.csv"
        # generate the dataframe
        df = pd.DataFrame.from_dict(combo, orient = "index")
        # write the dataframe to the output file
        df.to_csv(outfile, sep = csv_sep)


def write_common_edges_csvs(common_edges = None,
                            psngroup = None,
                            outfiles_prefix = "edges_",
                            node_sep = "_",
                            psn_sep = "_",
                            csv_sep = ",",
                            float_fmt = "%2.3f",
                            **kwargs):
    """Write the common edges for each possible subset of PSNs
    in the PSNgroup to a CSV file.

    Parameters
    ----------
    common_edges : `dict` or `None`, default: `None`
        Common edges, if already calculated.

    psngroup : `psntools.PSN.PSNGroup` or `None`, default: `None`
        PSN group.

    outfiles_prefix : `str`, default: `"edges_"`
        Prefix for the output files.

    node_sep : `str`, default: `"_"`
        How to separate nodes in each edge.

    psn_sep : `str`, default: `"_"`
        How to separate the PSN names in the names of the
        output files.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV files.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV files.

    **kwargs
        Arguments to be passed to `psngroup.get_common_edges`
        for common edges calculation, if `common_edges`
        has not been passed.

    Returns
    -------
    `None`
    """

    # if the common edges were not passed
    if common_edges is None:
        # get the common edges
        common_edges = psngroup.get_common_edges(**kwargs)

    # for each PSN combination
    for combo_label, combo in common_edges.items():
            
        # replace tuples representing edges with strings
        combo = {p : {node_sep.join(e) : w for e,w in values.items()}
                 for p, values in combo.items()}
            
        # get the specific name of the file
        name = psn_sep.join(combo_label)    
        # add prefix and extension
        outfile = f"{outfiles_prefix}{name}.csv"     
            
        # generate the dataframe
        df = pd.DataFrame.from_dict(combo, orient = "index")         
            
        # write the dataframe to the output file
        df.to_csv(outfile,
                  sep = csv_sep,
                  float_format = float_fmt)


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

    psngroup : `psntools.PSN.PSNGroup` or `None`, default: `None`
        PSN group.

    metric : `str` or `None`, default: `None`
        Label of a single node metric.

    csv_sep : `str`, default: `","`
        Field separator in the output CSV file.

    float_fmt : `str`, default: `"%2.3f"`
        Format for floating point numbers in the output CSV file.

    Returns
    -------
    `None`
    """

    # if the dataframe was not passed
    if df is None:
        # get the dataframe
        df = get_nodes_df_psngroup(metric = metric)

    # write the CSV file
    df.to_csv(outfile,
              sep = csv_sep,
              float_format = float_fmt)