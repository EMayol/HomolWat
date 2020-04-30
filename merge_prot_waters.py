"""
merge output of HomoWat
order: protein + (MOL) + (other ligands) + (NA) + WATS
returns 2 files:
    model_HW.pdb, model plus added waters 
    model_info_HW.csv, information about waters incorporated
Homolwat 5 of 6
__author__: Eduardo Mayol
"""

import os
import sys
import glob
import shutil
import numpy as np
from EM_functions import Euc_dist


def file_exist(path):
    if os.path.isfile(path):
        p = True
    else:
        p = False
    return p

def get_coords(atom_line):
    coxA = atom_line[30:38].strip()
    coyA = atom_line[38:46].strip()
    cozA = atom_line[47:54].strip()
    coord_A = np.array([coxA, coyA, cozA], float)
    return coord_A

# path = "/people/common/add_water/"
path = sys.argv[1]
query_name = sys.argv[2]

f_in_query = open(path + "protein-wat.pdb", "r")
f_r_query = f_in_query.readlines()
f_in_query.close()

# chek if dowser has been run
if file_exist(path[:-3] + "waters/wats_HWDp.pdb"):
    f_in_watsH = open(path[:-3] + "waters/wats_HWDp.pdb", "r")
    f_r_watsH = f_in_watsH.readlines()
    f_in_watsH.close()
else:
    f_r_watsH = []


protein, mol, sodium, waters, wat_pdb_info = [], [], [], [], []
# save prot and wats WITH Hs
for l1 in f_r_query:
    if l1.startswith("ATOM") or l1.startswith("TER"):
        protein.append(l1)
    elif l1.startswith("HETATM"):
        tipo = l1[17:20].strip()
        if tipo != "HOH":
            if tipo == "NA":
                sodium.append(l1)
            else:
                mol.append(l1)
        else:
            wat_pdb_info.append(l1)

# Dowser waters are renamed and the other are kept from protein-wat

if f_r_watsH != []:
    tmp_wats = []
    fin_line_O = "  1.00 99.99           O "
    for w_d in f_r_watsH:
        wcoords = w_d[30:54].strip()
        at_name = w_d[12:16].strip()
        if at_name == "O":
            for w_p in wat_pdb_info:
                w_p_coords = w_p[30:54].strip()
                if w_p_coords == wcoords:
                    tmp_wats.append(w_p)
        elif at_name == "OW":
            dw_coords = get_coords(w_d)
            new_line = w_d[:13] + "O " + w_d[15:54] + fin_line_O + "dowser\n"
            # avoid clashes with MOL and NA
            if mol != []:
                distances = []
                for mol_at in mol:
                    mol_coords = get_coords(mol_at)
                    dist = Euc_dist(dw_coords, mol_coords)
                    distances.append(dist)
                if min(distances) > 2.4:
                    if sodium != []:
                        for na_line in sodium:
                            na_coord = get_coords(na_line)
                            dist_na_wat = Euc_dist(dw_coords, na_coord)
                            if dist_na_wat > 2.1:
                                tmp_wats.append(new_line)
                    else:
                        tmp_wats.append(new_line)
            else:
                if sodium != []:
                    for na_line in sodium:
                        na_coord = get_coords(na_line)
                        dist_na_wat = Euc_dist(dw_coords, na_coord)
                        if dist_na_wat > 2.1:
                            tmp_wats.append(new_line)
                else:
                    tmp_wats.append(new_line)
    wat_pdb_info = tmp_wats

print("info num aguas:", len(wat_pdb_info))

###other ligands
other_ligands = []
if file_exist(path + "other_ligs.pdb"):
    f_in_oth = open(path + "other_ligs.pdb")
    f_r_oth = f_in_oth.readlines()
    for line in f_r_oth:
        other_ligands.append(line)

### Delete waters clashing with other ligands heavy atoms

new_wats_renum = []
delete_wats = []
for wat in wat_pdb_info:
    at_name = wat[12:16].strip()
    wat_num = wat[22:26].strip()
    wat_info = wat[22:54].strip()
    coor_wat = [
        float(wat[29:38].strip()),
        float(wat[38:46].strip()),
        float(wat[46:54].strip()),
    ]
    if at_name == "O":
        for atm in other_ligands:
            atm_name = atm[12:16].strip()
            if not atm_name.startswith("H"):
                coor_oth = [
                    float(atm[29:38].strip()),
                    float(atm[38:46].strip()),
                    float(atm[46:54].strip()),
                ]
                if Euc_dist(coor_wat, coor_oth) < 2.4:
                    delete_wats.append(wat_info)

uniq_delete = set(delete_wats)
if uniq_delete != []:
    for wat in wat_pdb_info:
        wat_info = wat[22:54].strip()
        if wat_info not in uniq_delete:
            new_wats_renum.append(wat)
print("without clashes:", len(new_wats_renum))

uniq_waters = []
for w1 in new_wats_renum:
    if w1 not in uniq_waters:
        uniq_waters.append(w1)

print("without repeats:", len(uniq_waters))

# RENUM WATERS and Generate file wats4colors for NGL

cont = len(sodium)

wats4colors, final_wats, water_numbers = [], [], []
for w in uniq_waters:
    old_watnum = w[22:26].strip()
    new_line = w
    num_a = str(4000 + cont)
    num_b = str(400 + cont)
    water_nums = (num_b, old_watnum)
    water_numbers.append(water_nums)
    new_line = w[:7] + num_a + w[11:21] + "  " + num_b + w[26:]
    cont += 1
    # watcol = (num_b, w[-5:])
    if w[-4:-1].replace(" ", "") == "O":
        watcol = (num_b, w[-5:-1])
    elif w[-7:-1].replace(" ", "") == "dowser":
        watcol = (num_b, w[-7:-1])
    else:
        wend = w[79:]
        elemwat = wend.split()
        pdbfam = elemwat[0].strip()
        elempdb = pdbfam.split("-")
        family = elempdb[0]
        pdbcode = elempdb[1]
        code = pdbcode + "_" + family
        watcol = (num_b, code)
    wats4colors.append(watcol)
    water = str(new_line)
    final_wats.append(water)

water_info = open(path + "wats4colors.txt", "w")
for elem in wats4colors:
    water_info.write(str(elem[0]) + ":" + str(elem[1]) + "\n")


# save header, title, cryst1 if exist
head, title, cryst1 = False, False, False
header_out, title_out, cryst1_out = [], [], []
if file_exist(path + "info_header.pdb"):
    f_in_head = open(path + "info_header.pdb", "r")
    f_r_head = f_in_head.readlines()
    f_in_head.close()
    if f_r_head != []:
        for l_head in f_r_head:
            if l_head.startswith("HEADER"):
                header_out.append(l_head)
                head = True
            elif l_head.startswith("TITLE"):
                title_out.append(l_head)
                title = True
            elif l_head.startswith("CRYST1"):
                cryst1_out.append(l_head)
                cryst1 = True

if head == False:
    new_head = "HEADER    " + query_name + " waters incorpored using HomolWat\n"
    header_out.append(new_head)

if cryst1 == False:
    new_cryst = (
        "CRYST1    1.0      1.0      1.0    90     90     90    P 1           1 \n"
    )
    cryst1_out.append(new_cryst)

# join all lists
system = (
    header_out
    + title_out
    + cryst1_out
    + protein
    + mol
    + other_ligands
    + sodium
    + final_wats
)

# save system merged
f_out_model = open(path[:-3] + query_name + "_HW.pdb", "w")
for line in system:
    f_out_model.write(line)
f_out_model.write("END")
f_out_model.close()

### END save solvated MODEL:

## INIT MODEL CSV

## gene table info output | open files with info needed
receptores, pdbs_used, wats_used = [], [], []
for i, wat_out in enumerate(final_wats):
    wat_num = wat_out[22:26].strip()
    oldwatnum = water_numbers[i][1]
    print(wat_num, oldwatnum, water_numbers[i][0])
    watend = wat_out[79:]
    if wat_out[-4:-1].replace(" ", "") == "O":
        info_w = (wat_num, oldwatnum, 'from model', '', '')
    elif wat_out[-7:-1].replace(" ", "") == "dowser":
        info_w = (wat_num, '', 'from dowser', '', '')
    else:
        elemwat = watend.split()
        local_rmsd = elemwat[1].strip()
        pdbfam = elemwat[0].strip()
        elempdb = pdbfam.split("-")
        family = elempdb[0]
        pdbcode = elempdb[1]
        info_w = (wat_num, oldwatnum, pdbcode[:4], pdbcode[-1], local_rmsd)
        
        if family not in receptores:
            receptores.append(family)
        if pdbcode not in pdbs_used:
            pdbs_used.append(pdbcode)
    wats_used.append(info_w)


print(receptores)
print(pdbs_used)
print("waters:", len(wats_used))

# open files with info about RMSD and Blast score

f_in_rmsds = open(path[:-3] + "info_rmsds.txt", "r")
f_r_rmsds = f_in_rmsds.readlines()
f_in_rmsds.close()

f_in_blast = open(path[:-3] + "rec_sort_list.csv", "r")
f_r_blast = f_in_blast.readlines()
f_in_blast.close()


# save info blast and get pdbs
print("blast-info:")
blast_vals = []
for recept in receptores:
    for line_b in f_r_blast:
        if line_b != f_r_blast[0]:
            elem_b = line_b.split(",")
            uniprot = elem_b[0]
            blastp = elem_b[1]
            evalue = elem_b[4]
            pdb_list = elem_b[5]
            pdbs = pdb_list.split()
            if recept == uniprot:
                info_blast = (uniprot, blastp, evalue, pdbs)
                print(info_blast)
                blast_vals.append(info_blast)
            # info_pdbs = (uniprot, pdbs)
            # pdbs_blast.append(info_pdbs)

# save info rmsds
rmsd_values = []
print("Globalrmsd-info:")
for pdbid in pdbs_used:
    pdb_in = pdbid[:4]
    chain_in = pdbid[-1]
    # print(pdb_in, chain_in)
    for line_r in f_r_rmsds:
        elem_r = line_r.split()
        pdbcode = elem_r[0][:4]
        chain = elem_r[0][-1]
        rmsd = elem_r[1]
        # print("compare to:", pdbcode, chain, rmsd)
        if pdb_in == pdbcode and chain_in == chain:
            info_rmsd = (pdbcode, chain, rmsd)
            print(info_rmsd)
            rmsd_values.append(info_rmsd)

# merge info in table
table_out1 = []
for line_o in rmsd_values:
    pdb_out = line_o[0]
    chain_out = line_o[1]
    globrmsd_out = line_o[2]
    for line2 in blast_vals:
        unip_out = line2[0]
        blastp_out = line2[1]
        evalue_out = line2[2]
        pdb_list_blast = line2[3]
        for pdbid in pdb_list_blast:
            if pdb_out == pdbid:
                info_out = (
                    pdb_out,
                    chain_out,
                    unip_out,
                    globrmsd_out,
                    blastp_out,
                    evalue_out,
                )
                print(info_out)
                if info_out not in table_out1:
                    table_out1.append(info_out)

# for each water look at local rmsd
# info_w = (wat_num, pdbcode[:4], pdbcode[-1], wat_coords)

table_out2 = []
print("para tabla final")
for wat in wats_used:
    print(wat)
    wnum_new = wat[0]
    wnum_old = wat[1]
    wpdb = wat[2]
    wchain = wat[3]
    wlocalrmsd = wat[4]
    if wpdb.startswith('from'):
        info_out2 = (wnum_new, wpdb, '','', '', '', '', '')
        table_out2.append(info_out2)
    else:
        for row in table_out1:
            if wpdb == row[0] and wchain == row[1]:
                info_out2 = (wnum_new, wpdb, wchain, row[2], wlocalrmsd, row[3], row[4], row[5])
                table_out2.append(info_out2)



f_out_info = open(path[:-3] + query_name + "_info_HW.csv", "w")
f_out_info.write("WaterNum,PDBid,Chain,Receptor,Local RMSD,Global RMSD,Blastp score,e-value\n")
for element in sorted(table_out2):
    line_out = (
        element[0]
        + ","
        + element[1]
        + ","
        + element[2]
        + ","
        + element[3]
        + ","
        + element[4]
        + ","
        + element[5]
        + ","
        + element[6]
        + ","
        + element[7]
        + "\n"
    )
    f_out_info.write(line_out)
f_out_info.close()

### END MODEL INFO CSV

### INIT PSE

### copy files to generate pse
os.mkdir(path + "pse")
path_prot = path[:-3]
path_pse = path + "pse/"
shutil.copyfile(path_prot + 'protein.pdb', path_pse + 'protein.pdb')
# folder_wats = glob.glob(path_waters + "*.pdb")

watnums = []
for p1 in pdbs_used:
    wat_pdbch  = []
    for w2 in final_wats:
        watnum = w2[22:26].strip()
        wend2 = w2[79:].strip()
        if wend2 != '' and wend2 != '\n' and wend2 != 'dowser':
            elemwat2 = wend2.split()
            pdbfam2 = elemwat2[0].strip()
            elempdb2 = pdbfam2.split("-")
            family2 = elempdb2[0]
            pdbcode2 = elempdb2[1]
            if pdbcode2 == p1:
                wat_pdbch.append(w2)
                watnums.append(watnum)
    f_out= open(path_pse + p1 + "_wats.pdb", "w")
    for w3 in wat_pdbch:
        f_out.write(w3)
    f_out.close()

wats_prot, wats_dowser = [], []
for w4 in final_wats:
    if w4[-4:-1].replace(" ", "") == "O":
        wats_prot.append(w4)
    elif w4[-7:-1].replace(" ", "") == "dowser":
        wats_dowser.append(w4)

print("watsprot", len(wats_prot))
print("watsdowser", len(wats_dowser))

if wats_prot != []:
    f_out_p = open(path_pse + "protein_wats.pdb", "w")
    for w5 in wats_prot:
        f_out_p.write(w5)
    f_out_p.close()

if wats_dowser != []:
    f_out_d = open(path_pse + "dowser_wats.pdb", "w")
    for w6 in wats_dowser:
        f_out_d.write(w6)
    f_out_d.close()
