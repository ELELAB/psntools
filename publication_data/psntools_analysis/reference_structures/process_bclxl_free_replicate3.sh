#!/bin/bash

in_pdb=bclxl_free_replicate3.pdb
out_pdb=bclxl_free_replicate3_processed.pdb

python3.7 assign_chain.py -i $in_pdb -o $out_pdb --chain A --start 1 --end 2643
