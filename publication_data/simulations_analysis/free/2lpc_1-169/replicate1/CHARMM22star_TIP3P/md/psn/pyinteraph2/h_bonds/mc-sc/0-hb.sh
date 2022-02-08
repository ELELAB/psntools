#!/bin/bash

# Trajectory, topology, reference structure

traj="../../../../../../../../../../simulations/free/2lpc_1-169/replicate1/CHARMM22star_TIP3P/md/Mol_An/center_traj.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"

# Options

ff_masses=charmm27
hb_perco=0
hb_class=mc-sc
hb_graph=hb-graph.dat
hb_ad_file=hydrogen_bonds.ini

pyinteraph -s $gro -t $traj -r $pdb --ff-masses $ff_masses -y --hb-perco $hb_perco --hb-class $hb_class --hb-graph $hb_graph --hb-ad-file $hb_ad_file -v

~                                                                               
