#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-



# Standard library
import os
import subprocess
# Third-party packages
import MDAnalysis as mda
import numpy as np
import pandas as pd
# psntools
import psntools.core as core
import psntools.dataframes as dataframes
import psntools.plotting as plotting
import psntools.upset as upset
import psntools.writing as writing



############################## SETTINGS ###############################



# Set the working directory to the current working directory
wd = os.getcwd()

# Directory containing the configuration files
config_wd = os.path.join(wd, "configuration_files")

# Interpreter used to running sub-processes
interpreter = "python"

# Script to generate the plot showing the size of the largest
# connected component as a function of the interaction strength
# (acPSN)
plot_ccs_imins_script = os.path.join(wd, "plot_ccs_imins.py")


#-------------------------------- PSNs -------------------------------#


# PSN types that have been calculated (i.e. acPSN, hydrogen bonds,
# salt bridges)
psn_types = \
    {"acPSN" : [("0.0_pers50", "acpsn")],
     "h_bonds" : [("sc-sc", "hb"), ("mc-mc", "hb"), ("mc-sc", "hb")],
     "s_bridges" : [("5.0", "sb")]}

# Path where the PSN matrices are stored
matrices_path = \
    "../simulations_analysis/{:s}/md/psn/pyinteraph2/" \
    "{:s}/{:s}/{:s}-graph_all.dat"

# Names of the directories containing data for each system
sysdirs = \
    ["free/2lpc_1-169/replicate1/CHARMM22star_TIP3P",
     "free/2lpc_1-169/replicate3/CHARMM22star_TIP3P",
     "beclin1/2PON_1-156_21-43/replicate1/CHARMM22star_TIP3P",
     "beclin1/2PON_1-156_21-43/replicate2/CHARMM22star_TIP3P",]

# Labels for the systems
labels = \
    ["bclxl_free_replicate1", "bclxl_free_replicate3",
     "bclxl_beclin1_replicate1", "bclxl_beclin1_replicate2"]

# Labels for the systems (UpSet plots)
labels_upset = \
    ["Bcl-xL (1)", "Bcl-xL (2)",
     "Bcl-xL/Beclin-1 (1)", "Bcl-xL/Beclin-1 (2)"]

# Path where the PDBs of the systems are stored
pdbs_path = os.path.join(wd, "reference_structures/{:s}_processed.pdb")

# Create the Universe objects from the PDBs
universes = [mda.Universe(pdbs_path.format(d)) for d in labels]

# Path to the CSV files containing the mapping between nodes in all PSNs
mappings = "node_mappings.csv"
mappings_upset = "node_mappings_upset.csv"

# Node metrics to be computed
metrics = \
    [{"degree" : {}}, 
     {"betweenness_centrality" : {}}, 
     {"closeness_centrality" : {}}]

# Minimum degree for a node to be considered a hub
min_degree = 3

# Minimum weight for an edge to be reported in PSNs weighted on
# the contacts' persistence
min_weight = 20.0

# Range of i_mins to be probed to construct the acPSN
i_mins_range = np.arange(0.0, 5.1, 0.1)


#---------------------------- Data frames ----------------------------#


# Separator to be used for the nodes in data frames
node_sep = "_"

# Separator to be used for the PSNs in data frames
psn_sep = "_"

# Separator to be used in output CSV files
csv_sep = ","

# Format for floating point numbers to be used in the output CSV files
float_fmt = "%2.3f"

# Format for the nodes in data frames
node_fmt = "strings"

# Intersection mode to be used when finding common hubs/edges
inters_mode = "distinct"


#---------------------------- UpSet plots ----------------------------#



# Configuration file for UpSet plot aesthetics for hubs
configfile_upset_hubs = \
    os.path.join(config_wd, "upsetplot_hubs.yaml")

# Configuration file for UpSet plot aesthetics for inter-chain edges
configfile_upset_interchain = \
    os.path.join(config_wd, "upsetplot_interchain_edges.yaml")

# Configuration file for UpSet plot aesthetics for intra-chain edges
configfile_upset_intrachain = \
    os.path.join(config_wd, "upsetplot_intrachain_edges.yaml")

# Configuration file for UpSet plot aesthetics for intra-chain edges
# in acPSN (different from the others, some aesthetics need to be
# further fine-tuned since this figure is used in the main text)
configfile_upset_intrachain_acpsn = \
    os.path.join(config_wd, "upsetplot_intrachain_edges_acpsn.yaml")



############################### ANALYSIS ##############################



# For each PSN type
for t, st_n in psn_types.items():

    # For each PSN subtype
    for st, n in st_n:


        #----------------- PSNs and PSNGroup creation ----------------#


        # Set the current working directory and create it, if
        # it does not exist
        curr_wd = os.path.join(wd, f"{t}/{st}")
        os.makedirs(curr_wd, exist_ok = True)

        # For all systems, get the PSN matrices corresponding to the
        # current PSN type and subtype
        matrices = [matrices_path.format(d, t, st, n) for d in sysdirs]

        # Create an empty list to store the PSNs
        psns = []

        # Create an empty list to store the data frames
        # containing the size of the largest connected
        # component at each i_min value for each PSN
        num_nodes_largest_cc_dfs = []
        
        # For each matrix/universe/label/label for UpSet plots
        for m, u, label, label_upset in \
            zip(matrices, universes, labels, labels_upset):
            
            # Create the PSN
            psn = core.PSN(matrix = m, universe = u)
            
            # If the PSN is not an acPSN
            if t != "acPSN":

                # Select the edges based on the minimum
                # weight (persistence) defined
                psn.select_edges(min_weight = min_weight)
            
            # If the PSN is an atomic-contacts-based PSN
            else:

                # Create an empty list to store the size (number
                # of nodes) of the largest connected component
                # found at each i_min
                num_nodes_largest_cc = []
                
                # For each i_min
                for i_min in i_mins_range:
                    
                    # Create a copy of the PSN so that the original
                    # is not modified (can't use copy functions
                    # because MDAnalysis instances apparently do
                    # not support them)
                    psn_copy = core.PSN(matrix = m, universe = u)

                    # On the copy of the PSN, select only those
                    # edges whose weight is higher than the i_min
                    psn_copy.select_edges(min_weight = i_min)
                    
                    # Calculate the connected components of this
                    # filtered PSN
                    ccs = \
                        psn_copy.get_connected_components(\
                            ascending = False)

                    # Append the i_min and the size of
                    # the largest connected component found
                    # to the list
                    num_nodes_largest_cc.append(\
                        ("{:1.1f}".format(i_min), len(ccs[0])))

                # Create a data frame from the list
                df = pd.DataFrame(\
                        num_nodes_largest_cc,
                        columns = ["i_min", label_upset])

                # Set the i_min column as index
                df = df.set_index("i_min")

                # Store the data frame in the list of data frames
                num_nodes_largest_cc_dfs.append(df)

                # Filter by interaction strength according to the
                # values found looking at the connected components'
                # distributions
                if label == "bclxl_free_replicate1":
                    psn.select_edges(min_weight = 3.1)
                elif label == "bclxl_free_replicate3":
                    psn.select_edges(min_weight = 3.1)
                elif label == "bclxl_beclin1_replicate1":
                    psn.select_edges(min_weight = 2.7)
                elif label == "bclxl_beclin1_replicate2":
                    psn.select_edges(min_weight = 2.4)

            # Save the PSN
            psns.append(psn)
            
            # Get a data frame containing the edges of the PSN
            edges_df = dataframes.get_edges_df(psn = psn,
                                               sort_by = "node",
                                               ascending = True)
            
            # Write out the data frame
            outfile = os.path.join(curr_wd, f"edges_{label}.csv")
            writing.write_edges_csv(df = edges_df,
                                    outfile = outfile)

        # If there are data frames about the size of the largest
        # connected component at different i_mins (acPSN)
        if num_nodes_largest_cc_dfs:

            # Concatenate the data frames          
            final_df = pd.concat(num_nodes_largest_cc_dfs,
                                 axis = 1)

            # Save the data frame to a CSV file
            outfile = os.path.join(curr_wd, "num_nodes_largest_cc.csv")
            final_df.to_csv(outfile, sep = ",")

            # Set the colors to be used in the plot
            colors_ccs_imins = \
                ["#5e5ce6", "#40c8e0", "#ff9d0a", "#ffd60a"]

            # Set the output file that will contain the plot
            outfile_plot = os.path.join(curr_wd, "ccs_imins.pdf")

            # Interaction thresholds used for the different acPSNs
            # (will be highlighted in the plot with vertical lines)
            imins_vlines = [3.1, 3.1, 2.7, 2.4]

            # Generate the plot showing the behavior of the
            # size of the largest connected component as a
            # function of the interaction strength threshold
            subprocess.run([interpreter,
                            plot_ccs_imins_script,
                            "-i", outfile,
                            "-o", outfile_plot,
                            "-c", *colors_ccs_imins,
                            "--imins-vlines", *imins_vlines,
                            "--imin-min", "0.0",
                            "--imin-max", "5.0",
                            "--imin-spacing", "1.0",
                            "--cc-size-min", "0",
                            "--cc-size-max", "150",
                            "--cc-size-spacing", "50"])

        # Create the PSNGroup from the PSNs
        pg = core.PSNGroup(psns = psns,
                           labels = labels,
                           mappings = mappings)

        # Create the PSNGroup from the PSNs (for the UpSet plots)
        pg_upset = core.PSNGroup(psns = reversed(psns),
                                 labels = reversed(labels_upset),
                                 mappings = mappings_upset)

        # Create the PSNGroup from the PSNs (for the per-node metrics)
        pg_metric = core.PSNGroup(psns = psns,
                                  labels = labels_upset,
                                  mappings = mappings_upset)


        #------------------------ Node metrics -----------------------#


        # For each metric
        for metric in metrics:

            # Set the output names for the data frame containing
            # the metric values for all PSNs in the PSN group and
            # for the corresponding heatmap
            outfile = \
                os.path.join(curr_wd, 
                             f"{t}_{st}_{list(metric.keys())[0]}.csv")
            outfile_plot = \
                os.path.join(curr_wd, 
                             f"{t}_{st}_{list(metric.keys())[0]}.pdf")

            # Get the path to the correct configuration file
            # for plotting
            configfile = \
                os.path.join(config_wd,
                             f"heatmap_nodes_{list(metric.keys())[0]}.yaml")

            # Generate the data frame containing the value of the
            # metric for all PSNs
            df = dataframes.get_nodes_df_psngroup(psngroup = pg_metric,
                                                  metric = metric)

            # Save the data frame containing the metric
            # for all PSNs
            writing.write_nodes_csv_psngroup(df = df,
                                             outfile = outfile,
                                             csv_sep = csv_sep,
                                             float_fmt = float_fmt)

            # Generate a heatmap from the data frame
            plotting.plot_heatmap_nodes(df = df,
                                        outfile = outfile_plot,
                                        configfile = configfile)

            # Only for the degree of side chain - side chain
            # hydrogen bonds
            if t == "h_bonds" and st == "sc-sc" \
                and list(metric.keys())[0] == "degree":

                # Create a list of nodes of interest
                selected_nodes = \
                    ["A-22TYR", "A-71GLN", "A-75THR",
                     "A-78THR", "A-96ASN", "A-137HIS"]

                # Create a list of the labels to be used for such
                # nodes of interest
                node_labels = \
                    ["Y22", "Q111", "T115",
                     "T118", "N136", "H177"]
                
                # Set the path to the output file that will
                # contain the plot
                outfile_plot_subset = os.path.join(\
                    curr_wd,
                    f"{t}_{st}_{list(metric.keys())[0]}_selected.pdf")

                # Set the configuration file to be used for plotting
                configfile_subset = os.path.join(\
                    config_wd,
                    f"heatmap_nodes_{list(metric.keys())[0]}_selected.yaml")

                # Generate a heatmap from the data frame
                plotting.plot_heatmap_nodes(\
                    df = df,
                    outfile = outfile_plot_subset,
                    configfile = configfile_subset,
                    selected_nodes = selected_nodes,
                    node_labels = node_labels)


        #------------------------ Common hubs ------------------------#


        # Get the data frames of common hubs
        common_hubs_dfs = \
            dataframes.get_common_hubs_dfs(\
                psngroup = pg,
                psn_sep = psn_sep,
                min_degree = min_degree,
                inters_mode = inters_mode,
                node_fmt = node_fmt)

        # Filter out the empty data frames
        common_hubs_dfs = \
            {label : df for label, df \
             in common_hubs_dfs.items() if not df.empty}

        # If there are non-empty dataframes
        if common_hubs_dfs:

            # Create the directory where the results of the analysis
            # of common hubs will be stored and enter it
            common_hubs_dir = os.path.join(curr_wd, "common_hubs")
            os.makedirs(common_hubs_dir, exist_ok = True)
            os.chdir(common_hubs_dir)

            # Write out CSV files for all common hubs
            writing.write_common_hubs_csvs(\
                dfs = common_hubs_dfs,
                outfiles_prefix = f"hubs_{t}_{st}_",
                csv_sep = csv_sep,
                float_fmt = float_fmt)

            # Set the output file that will contain the UpSet plot
            p_hubs_out = \
                os.path.join(common_hubs_dir, 
                             f"hubs_{t}_{st}_upset.pdf")

            # Generate the UpSet plot
            plotting.plot_upsetplot(\
                psngroup = pg_upset,
                item_type = "hubs",
                outfile = p_hubs_out,
                configfile = configfile_upset_hubs)


        #------------------ Common inter-chain edges -----------------#


        # Get the data frames of common interchain edges
        common_edges_interchain_dfs = \
            dataframes.get_common_edges_dfs(\
                psngroup = pg,
                psn_sep = psn_sep,
                node_sep = node_sep,
                mode = "interchain",
                inters_mode = inters_mode,
                node_fmt = node_fmt)

        # Filter out the empty data frames
        common_edges_interchain_dfs = \
            {label : df for label, df \
             in common_edges_interchain_dfs.items() if not df.empty}

        # If there are non-empty data frames
        if common_edges_interchain_dfs:

            # Create the directory where the results of the analysis 
            # of common inter-chain edges will be stored and enter it
            common_edges_inter_dir = \
                os.path.join(curr_wd, "common_edges_interchain")
            os.makedirs(common_edges_inter_dir, exist_ok = True)
            os.chdir(common_edges_inter_dir)

            # Write out CSV files for all common interchain edges
            writing.write_common_edges_csvs(\
                dfs = common_edges_interchain_dfs,
                outfiles_prefix = f"edges_inter_{t}_{st}_",
                csv_sep = csv_sep,
                float_fmt = float_fmt)

            # Set the output file that will contain the UpSet plot
            p_edges_inter_out = \
                os.path.join(common_edges_inter_dir, 
                             f"edges_inter_{t}_{st}_upset.pdf")

            # Generate the UpSet plot
            plotting.plot_upsetplot(\
                psngroup = pg_upset,
                item_type = "edges",
                outfile = p_edges_inter_out,
                configfile = configfile_upset_interchain,
                mode = "interchain")


        #------------------ Common intra-chain edges -----------------#


        # Get the data frames of common intrachain edges
        common_edges_intrachain_dfs = \
            dataframes.get_common_edges_dfs(\
                psngroup = pg,
                psn_sep = psn_sep,
                node_sep = node_sep,
                mode = "intrachain",
                node_fmt = node_fmt)

        # Filter out empty data frames        
        common_edges_intrachain_dfs = \
            {label : df for label, df \
             in common_edges_intrachain_dfs.items() if not df.empty}

        # If there are non-empty data frames
        if common_edges_intrachain_dfs:

            # Create the directory where the results of the analysis 
            # of common intra-chain edges will be stored and enter it
            common_edges_intra_dir = \
                os.path.join(curr_wd, "common_edges_intrachain")
            os.makedirs(common_edges_intra_dir, exist_ok = True)
            os.chdir(common_edges_intra_dir)

            # Write out CSV files for all common intrachain edges
            writing.write_common_edges_csvs(\
                dfs = common_edges_intrachain_dfs,
                outfiles_prefix = f"edges_intra_{t}_{st}_",
                csv_sep = csv_sep,
                float_fmt = float_fmt)
            
            # Set the output file that will contain the UpSet plot
            p_edges_intra_out = \
                os.path.join(common_edges_intra_dir, 
                             f"edges_intra_{t}_{st}_upset.pdf")

            # If it is the acPSN
            if t == "acPSN" and st == "0.0_pers50":
                
                # Generate the UpSet plot of the intra-chain edges
                # of acPSN
                plotting.plot_upsetplot(\
                    psngroup = pg_upset,
                    item_type = "edges",
                    outfile = p_edges_intra_out,
                    configfile = configfile_upset_intrachain_acpsn,
                    mode = "intrachain")

            # Otherwise
            else:
                
                # Generate the UpSet plot
                plotting.plot_upsetplot(\
                    psngroup = pg_upset,
                    item_type = "edges",
                    outfile = p_edges_intra_out,
                    configfile = configfile_upset_intrachain,
                    mode = "intrachain")


        # Go back to the working directory
        os.chdir(wd)
