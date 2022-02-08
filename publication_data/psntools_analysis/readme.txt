# In order to run the psntools analyses, you need these directories/files in place:
#
# - reference_structures -> Directory where the PDB reference structures
#			    used to generate the Universe objects that
#			    will be passed to the PSN constructor are
#			    stored.
#
# - configuration files -> Directory where the configuration files used
# 			   by psntools to generate the plots are stored.
#			   You will need to specify your own font (or
#			   leave it '!!null') in each of the configuration
#			   files, since we cannot redistribute the font
#			   that we used for the figures (Nexa Light) and
#			   the plotting utilities will try to find that font
#			   if the analysis is run without modifying the
#			   configuration files first.
#
# - node_mappings.csv -> File containining the one-to-one node mapping for
#			 the creation of the PSNgroup from the PSNs of the
#			 ensembles.
#
# - node_mappings_upset.csv -> File containing the one-to-one node mapping for
#			       the creation of the PSNgroup from the PSNs of the
#			       ensembles used for gnerating UpSet plots (different
#			       labeling of the ensembles with respect to
#			       node_mappings.csv).
#
# - run_analysis.py -> Python script responsible for running the analyses
#		       automatically.
#
# - plot_ccs_imins.py -> Python script responsible for plotting the size
#			 of the largest connected component of each acPSN
#			 as a function of the interaction strength's
#			 cut-off used (it is called by run_analysis.py).

# You need the latest version of psntools (as per 01-Feb-2022, commit 1b47833)
# installed to run the analyses.

# Run the analyses

python run_analysis.py

# The run will create the following directories (and subdirectories) and files:
#
# - acPSN -> Directory containing the analyses for acPSNs.
# - h_bonds -> Directory containing the analyses for hydrogen bonds networks.
# - s_bridges -> Directory containing the analyses for salt bridges networks.
# - ccs_imins.pdf -> Plot of the size of the most populated connected component
#		     in acPSNs as a function of the interaction strength cut-off
#		     (used for panel A of figure 3).

# Other directories contain data used for the creation of the figures used
# in the publication:

# - ai_figures -> Adobe Illustrator and corresponding TIF files of the
#		  final figures.
# - pymol_figures -> PyMol renderings of panels A and B of figure 2.
# - pymol_sessions -> PyMol sessions used to generate panels A and B of figure 2.
 
