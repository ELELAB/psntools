#!/bin/bash

in_pdb=bclxl_beclin1_replicate2.pdb
out_pdb=bclxl_beclin1_replicate2_processed.pdb

python3.7 assign_chain.py -i $in_pdb -o tmp1 --chain A --start 352 --end 2801
python3.7 assign_chain.py -i tmp1 -o tmp2 --chain B --start 1 --end 351
python3.7 renumber_chain.py -i tmp2 -o $out_pdb --chain B --new-start 106

rm tmp1 tmp2
