#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    assign_chain.py
#
#    Assign a chain ID to a range of atoms in a PDB file.
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



# Create the argument parser
argparser = argparse.ArgumentParser()


# Add the arguments
i_help = "Input PDB file."
argparser.add_argument("-i",
                       dest = "in_pdb_file",
                       type = str,
                       help = i_help)

o_help = "Output PDB file."
argparser.add_argument("-o",
                       dest = "out_pdb_file",
                       type = str,
                       help = o_help)

chain_help = "Chain ID to be assigned."
argparser.add_argument("--chain",
                       dest = "chain",
                       type = str,
                       help = chain_help)

start_help = "Where the chain starts (atom number)."
argparser.add_argument("--start",
                       dest = "start",
                       type = int,
                       help = start_help)

end_help = "Where the chain ends (atom number)."
argparser.add_argument("--end",
                       dest = "end",
                       type = int,
                       help = end_help)

# Get the arguments
args = argparser.parse_args()
in_pdb_file = args.in_pdb_file
out_pdb_file = args.out_pdb_file
chain = args.chain
start = args.start
end = args.end

# Range of serial numbers of the atoms which will
# be assigned the chain ID (add a 1 to the
# because it is passed as closed, but it
# is open in Python 'range')
range_serials = list(range(start, end+1))

# Open the input and output PDB files
with open(in_pdb_file, "r") as f, \
     open(out_pdb_file, "w") as out:

    # For each line in the input file
    for line in f:
        
        # If it is an ATOM record
        if line.startswith("ATOM"):
            
            # Strip the newline character and get rid of
            # all empty strings when splitting the line
            split_line = \
                [i for i in line.rstrip("\n").split(" ") \
                 if i != ""]
            # Get the atom serial number
            atom_serial = int(split_line[1])
            # Create a copy of the line to be modified
            tmp_line = list(line)
            # If the serial number is inside the range of interest
            if atom_serial in range_serials:
                # Change the chain ID field
                tmp_line[21] = chain
            # Write out the modified ATOM record
            out.write("".join(tmp_line))
        
        # If the line contains other data
        else:
            # Just write the line as it is
            out.write(line)