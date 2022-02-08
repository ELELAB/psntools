#!/bin/bash

# Trajectory, topology, reference structure

traj="../../../../../../../../../../simulations/beclin1/2PON_1-156_21-43/replicate2/CHARMM22star_TIP3P/md/Mol_An/traj_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"

# Options

ff_masses=charmm27
sb_perco=0
sb_co=5.5
sb_mode=different_charge
sb_graph=sb-graph.dat
sb_cg_file=charged_groups.ini

# Run pyinteraph

pyinteraph -s $gro -t $traj -r $pdb --ff-masses $ff_masses -b --sb-co $sb_co --sb-perco $sb_perco --sb-mode $sb_mode --sb-graph $sb_graph --sb-cg-file $sb_cg_file -v
