#!/bin/bash

# Load GROMACS

. /usr/local/gromacs-2021.2/bin/GMXRC

# Trajectory, topology, index file

xtc="../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/Mol_An/center_traj.xtc"
tpr="../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/md.tpr"
ndx="protein.ndx"

# Create the index file containing only the 'Protein' group

gmx make_ndx -f $tpr -o $ndx << eof
keep 1
q
eof

# Extract the first structure from the trajectory

gmx trjconv -f $xtc -s $tpr -n $ndx  -b 0 -e 1 -o model0.pdb
