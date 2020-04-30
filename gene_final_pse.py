###############################################################################
"""
This script generates a pse from protein and an object for each pdb that
contributes to the incorporation of waters (just waters added) 

IMPRTANT: To be run with "pymol -cqr prep_PDBs.py working_dir path name_of_model"
returns 1 file:
    model_HW.pse, model plus added waters 
Homolwat 6 of 6
__author__: Eduardo Mayol

"""
###############################################################################


from pymol import cmd
import os
import sys
import glob

req_num = sys.argv[3]
path = sys.argv[4]
query_name = sys.argv[5]

path_req = path + req_num + "/"
path_pse = path_req + "HW/pse/"

# load pdbs in pymol and save as pse
cmd.do("cd " + path_pse)
cmd.do("load protein.pdb")

# load wats from protein
if os.path.isfile(path_pse + "protein_wats.pdb"):
    cmd.do("load protein_wats.pdb")

# load wats from dowser
if os.path.isfile(path_pse + "dowser_wats.pdb"):
    cmd.do("load dowser_wats.pdb")

# load other waters
folder = glob.glob(path_pse + "*wats.pdb")

for pdb in sorted(folder):
    pdbcode = pdb[len(path_pse) :]
    if pdbcode[:4] != "prot" and pdbcode[:4] != "dows":
        cmd.do("load " + pdbcode)

cmd.do("cd " + path_req)
cmd.do("save " + query_name + "_HW.pse")

cmd.do("quit")
