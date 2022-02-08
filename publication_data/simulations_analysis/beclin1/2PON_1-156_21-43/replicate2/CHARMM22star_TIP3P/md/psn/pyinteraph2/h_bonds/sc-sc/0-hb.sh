#!/bin/bash

# Trajectory, topology, reference structure

traj="../../../../../../../../../../simulations/beclin1/2PON_1-156_21-43/replicate2/CHARMM22star_TIP3P/md/Mol_An/traj_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"

# Options

ff_masses=charmm27
hb_perco=0
hb_class=sc-sc
hb_graph=hb-graph.dat
hb_ad_file=hydrogen_bonds.ini

pyinteraph -s $gro -t $traj -r $pdb --ff-masses $ff_masses -y --hb-perco $hb_perco --hb-class $hb_class --hb-graph $hb_graph --hb-ad-file $hb_ad_file -v

~                                                                               
