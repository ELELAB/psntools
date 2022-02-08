#!/bin/bash

# Trajectory, topology, reference structure

traj="../../../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/Mol_An/center_traj.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"

# Options

ff_masses=charmm27
acpsn_perco=0.5
acpsn_co=4.5
acpsn_proxco=1
acpsn_imin=3.0
acpsn_graph=acpsn-graph.dat
acpsn_nf_file=normalization_factors.ini

# Run pyinteraph

pyinteraph -s $gro -t $traj -r $pdb --ff-masses $ff_masses -a --acpsn-co $acpsn_co --acpsn-perco $acpsn_perco --acpsn-proxco $acpsn_proxco --acpsn-imin $acpsn_imin --acpsn-graph $acpsn_graph --acpsn-nf-file $acpsn_nf_file -v
