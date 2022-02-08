#!/bin/bash

# Load GROMACS

. /usr/local/gromacs-2021.2/bin/GMXRC

# Trajectory, topology, index file

xtc="../../../../../../../../simulations/beclin1/2PON_1-156_21-43/replicate1/CHARMM22star_TIP3P/md/Mol_An/traj_centered.xtc"
tpr="../../../../../../../../simulations/beclin1/2PON_1-156_21-43/replicate1/CHARMM22star_TIP3P/md/md.tpr"
ndx="protein.ndx"

# Create the index file containing only the 'Protein' group

gmx make_ndx -f $tpr -o $ndx << eof
keep 1
q
eof

# Extract the first structure from the trajectory

gmx trjconv -f $xtc -s $tpr -n $ndx  -b 0 -e 1 -o model0.pdb
