"""
In order to run HomolWat it is required to generate many files
this script will create some folders in a specific location
for that you shouls indicate where is going to run
This script will:
- prepare folders and model to run Homolwat
- from the model it will prepare protein.pdb ligand.pdb and wat.pdb
    protein.pdb will be translated to prot_dw.pdb where "unusual" residue names
will be replaced for stardar names.
Homolwat 1 of 6
__author__: Eduardo Mayol
"""

aa_nonStd = {
    "GLH": "GLU",
    "HIE": "HIS",
    "HID": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
    "CYP": "CYS",
    "ASH": "ASP",
    "MSE": "MET",
}

aa_std = [
    "CYS",
    "ASP",
    "SER",
    "GLN",
    "LYS",
    "ILE",
    "PRO",
    "THR",
    "PHE",
    "ASN",
    "GLY",
    "HIS",
    "LEU",
    "ARG",
    "TRP",
    "ALA",
    "VAL",
    "GLU",
    "TYR",
    "MET",
]

import sys
import os
import shutil
import numpy as np
import shutil
from EM_functions import *


for i in range(6):
    print(i, sys.argv[i])

# This program is meant to work in a sever
# however is adapted to run in local computers


# asign paths: defined at run_HomolWat.sh
pdb_in = sys.argv[1]  # name of the pdb to work with -> MODEL
path_prev = sys.argv[2]  # where the folders will be created
path = sys.argv[3]  # name of the sulfoder
path_pdbs = sys.argv[4]  # fordel where the MODEL is located
folder_name = sys.argv[5]  # name of the folder where all is done
path_scripts = sys.argv[6]  # where the scripts are located
print(folder_name)


### open MODEL and get name of file
f_in = open(path_pdbs + pdb_in, "r")
f_r = f_in.readlines()
f_in.close()
path_parts = pdb_in.split("/")
filename = path_parts[-1][:-4]


# make folder new request
path_newReq = path + folder_name
os.mkdir(path_newReq)
# and subfolders
os.mkdir(path_newReq + "/waters")
os.mkdir(path_newReq + "/REC")
path_HW = path_newReq + "/HW"
os.mkdir(path_HW)
path_dow = path_newReq + "/dowser"
os.mkdir(path_dow)

# You can chose if the model is in active or inactive state
# inactive by default
# to prioritize Active sctructures you should change INACT to ACT for this.
shutil.copy(path_scripts + "Pdbs_resol_INACT", path_newReq + "/Pdbs_resol_L")
# in "Pdbs_resol_INACT" the GPCR pdbs are sorted by family and resolution

# Get fasta of gpcr MODEL

pdb2fasta_save(f_r, path_HW)

# separate ATOMS and HETATM and save header and cryst1
pdb_info_lines = []
wat_ions = ["HOH", "SOL", "WAT", "NA", "CL"]
protein, ligand, wat_ion, dowprot = [], [], [], []
cbs_coords = []
for line in f_r:
    if line.startswith("HEADER"):
        pdb_info_lines.append(line)
    elif line.startswith("CRYST1"):
        pdb_info_lines.append(line)
    elif line.startswith("TITLE"):
        pdb_info_lines.append(line)
    elif line.startswith("ATOM") or line.startswith("TER"):
        protein.append(line)
        resname = line[17:20].strip()
        if resname in list(aa_nonStd.keys()):
            std_rsdname = aa_nonStd[resname]
            new_line = "ATOM  " + line[6:17] + std_rsdname + line[20:]
            dowprot.append(new_line)
        else:
            dowprot.append(line)
        tipo = line[12:16].strip()
        if tipo == "CB":
            coxA = line[30:38].strip()
            coyA = line[38:46].strip()
            cozA = line[47:54].strip()
            coord_i = np.array([coxA, coyA, cozA], float)
            cbs_coords.append(coord_i)
    elif line.startswith("HETATM"):
        molec = line[17:20].strip()
        if molec in wat_ions:
            wat_ion.append(line)
        elif molec in list(aa_nonStd.keys()):
            protein.append(line)
            std_rsd = aa_nonStd[molec]
            new_line = "ATOM  " + line[6:17] + std_rsd + line[20:]
            dowprot.append(new_line)
            tipo = line[12:16].strip()
            if tipo == "CB":
                coxA = line[30:38].strip()
                coyA = line[38:46].strip()
                cozA = line[47:54].strip()
                coord_i = np.array([coxA, coyA, cozA], float)
                cbs_coords.append(coord_i)
        else:
            ligand.append(line)

# save data in separated files
f_out_prot = open(path_newReq + "/protein.pdb", "w")
for l1 in protein:
    f_out_prot.write(l1)
f_out_prot.close()


### save dowser_prot in case it is used
f_out_dow = open(path_newReq + "/dowser/dow_protein.pdb", "w")
for l2 in dowprot:
    f_out_dow.write(l2)
f_out_dow.close()


# save header, title and cryst1
f_out_header = open(path_HW + "/info_header.pdb", "w")
for l3 in pdb_info_lines:
    f_out_header.write(l3)
f_out_header.close()

########################


int_mol, ligand_lines = [], []
if ligand != []:
    # vamos a guardar solo el ligando interno
    # 1o hacemos lista de los ligandos del pdb
    molecs = []
    for l in ligand:
        molec = l[17:20].strip()
        num = l[22:27].strip()
        info = [molec, num]
        if info not in molecs:
            molecs.append(info)
    # 2o sacamos el COM de cada uno
    for i, mol in enumerate(molecs):
        coords = []
        # print i, mol
        for l2 in ligand:
            molec = l2[17:20].strip()
            num = l2[22:27].strip()
            if molec == mol[0] and num == mol[1]:
                # print l2
                coxA = l2[30:38].strip()
                coyA = l2[38:46].strip()
                cozA = l2[47:54].strip()
                coord_i = np.array([coxA, coyA, cozA], float)
                coords.append(coord_i)
        com_mol = center_of_molec(coords)
        molecs[i].append(com_mol)
    # 3o sacamos el CV de la molecula
    for j, mol in enumerate(molecs):
        cv = check_CircVar(cbs_coords, mol[2], 20)
        molecs[j].append(cv)
    # for x in molecs:
    # print x
    # 4o ordenamos la lista de molecs por cv y nos quedamos el ultimo (max CV)
    sorted_mols = sorted(molecs, key=lambda cv: cv[3])
    # print "ordenados"
    int_mol_num = [sorted_mols[-1][0], sorted_mols[-1][1]]
    # print int_mol_num
    # 5o se guardan las lineas del ligando en una lista
    if len(molecs) > 1:
        f_out_other_ligands = open(path_HW + "/other_ligs.pdb", "w")
    for line in ligand:
        molec = line[17:20].strip()
        num = line[22:27].strip()
        if molec == int_mol_num[0] and num == int_mol_num[1]:
            # print line
            ligand_lines.append(line)
        else:
            f_out_other_ligands.write(line)
    # 6o guardamos el ligando en un archivo
    f_out_ligand_het = open(path_HW + "/ligand_het.pdb", "w")
    f_out_ligand_atom = open(path_HW + "/ligand_atom.pdb", "w")
    for l3 in ligand_lines:
        f_out_ligand_het.write(l3)
        atom_l3 = "ATOM  " + l3[6:]
        f_out_ligand_atom.write(atom_l3)
    f_out_ligand_het.close()
    f_out_ligand_atom.close()

# save waters and ions from model
if wat_ion != []:
    f_out_wat_ions = open(path_HW + "/wations.pdb", "w")
    for l4 in wat_ion:
        f_out_wat_ions.write(l4)
    f_out_wat_ions.close()


# RUN BLASTP to get the order of receptors
# extracts the list of sorted receptores (por orden filogenetico)

fasta_in = path_HW + "/HWprot.fasta"
path2save = path_newReq + "/"
path_wd = path_prev + "scripts/"
run_blastp(fasta_in, path_wd, path2save)

