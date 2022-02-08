#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    renumber_chain.py
#
#    Renumber a chain in a PDB file.
#
#    Copyright (C) 2020 Valentina Sora 
#                       <sora.valentina1@gmail.com>
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



# Standard library
import argparse
# Third-party packages
from Bio.PDB import PDBParser, PDBIO, Select



# Create the argument parser
argparser = argparse.ArgumentParser()

# Add the arguments
i_help = "Input PDB file"
argparser.add_argument("-i",
                       dest = "in_pdb_file",
                       type = str,
                       help = i_help)

o_help = "Output PDB file"
argparser.add_argument("-o",
                       dest = "out_pdb_file",
                       type = str,
                       help = o_help)

chain_help = "Chain to be renumbered"
argparser.add_argument("--chain",
                       dest = "chain",
                       type = str,
                       help = chain_help)

new_start_help = "New starting number"
argparser.add_argument("--new-start",
                       dest = "new_start",
                       type = int,
                       help = new_start_help)

# Get the arguments
args = argparser.parse_args()
in_pdb_file = args.in_pdb_file
out_pdb_file = args.out_pdb_file
chain = args.chain
new_start = args.new_start

# Create a PDB parser
parser = PDBParser()

# Parse the structure
name = in_pdb_file.replace(".pdb", "")
structure = parser.get_structure(name, in_pdb_file)

# Renumber the chain (for all models)
for mod in structure:
    for ch in mod:
        if ch.id == chain:
            i = new_start
            for res in ch:
                resid = res.id
                res.id = (resid[0], i, resid[2])
                i += 1

# Save the processed structure
w = PDBIO()
w.set_structure(structure)
w.save(out_pdb_file)