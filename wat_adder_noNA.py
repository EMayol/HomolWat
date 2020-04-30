"""
script to add waters WITHOUTinternal NA
usage ==> python wat_adder_noNA.py req_num
Homolwat 4A of 6
__author__: Eduardo Mayol
"""

import math
import os
import sys
import numpy as np

np.seterr(invalid="warn")
# distances by default
d_wr = 2.3
d_ww = 2.4
d_wn = 2.1
d_nr = 1.8
d_nn = 5.0
path = os.getcwd()

aa3 = [
    "ILE",
    "VAL",
    "LEU",
    "PHE",
    "TYR",
    "TRP",
    "MET",
    "PRO",
    "GLY",
    "ALA",
    "CYS",
    "SER",
    "THR",
    "ASN",
    "GLN",
    "ASP",
    "GLU",
    "LYS",
    "ARG",
    "HIS",
    "MSE",
    "ACE",
    "NME",
    "GLH",
    "HIE",
    "HID",
    "HIP",
    "CYX",
    "CYP",
    "ASH",
]


def Euc_dist(p, q):
    """ distance between 2 3D points"""
    dist = math.sqrt((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2)
    return dist


def check_CircVar(CB_coords, water_line):
    """ Compute Circular Variance from CBs of protein
    of water molecule at a radium of 10A"""
    rcv = 10.0
    # info_CVs = []
    coxB = water_line[30:38].strip()
    coyB = water_line[38:46].strip()
    cozB = water_line[47:54].strip()
    # w_num = water_line[23:26].strip()
    coord_w = np.array([coxB, coyB, cozB], float)
    sumvec = 0
    sumcount = 0
    for cb in CB_coords:
        coord_cb = np.array(cb, float)
        vec = coord_w - coord_cb
        modvec = np.linalg.norm(vec)
        if modvec < rcv:
            sumvec += vec / modvec
            sumcount += 1
    cv = 1 - (np.linalg.norm(sumvec)) / sumcount
    return cv


def center_of_mass(coord_list):
    """ Compute Center of a group of coordinates """
    x_list, y_list, z_list = [], [], []
    for coord in coord_list:
        x_list.append(float(coord[0]))
        y_list.append(float(coord[1]))
        z_list.append(float(coord[2]))
    x_temp = sum(x_list) / len(coord_list)
    y_temp = sum(y_list) / len(coord_list)
    z_temp = sum(z_list) / len(coord_list)
    cf = [("%.3f" % x_temp), ("%.3f" % y_temp), ("%.3f" % z_temp)]
    coord_fin = []
    for el in cf:
        coord_fin.append(float(el))
    return coord_fin


def file_exist(path):
    if os.path.isfile(path):
        p = True
    else:
        p = False
    return p


def add_h2o_to_model(h2o_i, pdb_model):
    """
    Main funtion to incorporate waters
    Checks if it is possible to add water to model
    Avoid clashes with resiudes and previos added waters
    returns True or False whether it fits inthe model
    """
    N_atom_list = ["N", "NZ", "ND1", "NE2", "NE", "NH1", "NH2", "NE1", "ND2"]
    O_atom_list = ["O", "OD1", "OD2", "OH", "OE1", "OE2", "OG1", "OG"]
    S_atom_list = ["SG"]  # SD de MET NO forma pdH
    hbond_acc_don = N_atom_list + O_atom_list + S_atom_list
    global d_wr, d_wn
    # tipo_w = h2o_i[17:20].strip()
    todas_dists = []
    polar_dists = []
    coxA = h2o_i[30:38].strip()
    coyA = h2o_i[38:46].strip()
    cozA = h2o_i[47:54].strip()
    coord_i = (float(coxA), float(coyA), float(cozA))
    x = 0
    for line in pdb_model:
        if line.startswith("ATOM"):
            tipo = line[12:16].strip()
            molec = line[17:20].strip()
            if not tipo.startswith("H"):
                x += 1
                coxB = line[30:38].strip()
                coyB = line[38:46].strip()
                cozB = line[47:54].strip()
                coord_j = (float(coxB), float(coyB), float(cozB))
                dist = Euc_dist(coord_i, coord_j)
                todas_dists.append(dist)
                # Avoid clash with heavy atoms
                if dist < d_wr:
                    return False
                # Filter polar atom dist
                elif tipo in hbond_acc_don:
                    if dist < 4:
                        polar_dists.append(dist)
        elif line.startswith("HETATM"):
            tipo = line[12:16].strip()
            molec = line[17:20].strip()
            if molec == "NA":
                coxB = line[30:38].strip()
                coyB = line[38:46].strip()
                cozB = line[47:54].strip()
                coord_j = (float(coxB), float(coyB), float(cozB))
                dist = Euc_dist(coord_i, coord_j)
                todas_dists.append(dist)
                if dist < d_wn:
                    return False
            elif molec == "HOH":
                coxB = line[30:38].strip()
                coyB = line[38:46].strip()
                cozB = line[47:54].strip()
                coord_j = (float(coxB), float(coyB), float(cozB))
                dist = Euc_dist(coord_i, coord_j)
                todas_dists.append(dist)
                # Avoid clash with previous waters
                if dist < d_ww:
                    return False
                # Filter polar atom dist
                elif dist < 4:
                    polar_dists.append(dist)
    # Apply filter polar dist
    if len(polar_dists) == 0:
        return False
    else:
        return True


def add_na_to_model(na_list, pdb_model):
    """
    first we need to compute de COM of the prot
    then place a sodium if is near D2.50, just 1
    """
    counter = 0
    global d_nr
    c_list = []
    for atom in pdb_model:
        tipo = atom[12:16].strip()
        if tipo == "CA":
            cox = atom[30:38].strip()
            coy = atom[38:46].strip()
            coz = atom[47:54].strip()
            coord = (float(cox), float(coy), float(coz))
            c_list.append(coord)
    COM_coord = center_of_mass(c_list)
    inner_sodiums = []
    for sod in na_list:
        coxZ = sod[30:38].strip()
        coyZ = sod[38:46].strip()
        cozZ = sod[47:54].strip()
        coord_Z = (float(coxZ), float(coyZ), float(cozZ))
        distance_to_COM = Euc_dist(COM_coord, coord_Z)
        if distance_to_COM < 10:
            inner_sodiums.append(sod)
    distances_n_r = []
    if inner_sodiums != []:
        for sodium in inner_sodiums:
            coxG = sodium[30:38].strip()
            coyG = sodium[38:46].strip()
            cozG = sodium[47:54].strip()
            coordG = (float(coxG), float(coyG), float(cozG))
            for atom2 in pdb_model:
                if atom2.startswith("ATOM"):
                    coxF = atom2[30:38].strip()
                    coyF = atom2[38:46].strip()
                    cozF = atom2[47:54].strip()
                    coordF = (float(coxF), float(coyF), float(cozF))
                    dist_at_na = Euc_dist(coordG, coordF)
                    if dist_at_na > d_nr:
                        counter += 1
                        return sodium
                    else:
                        distances_n_r.append(dist_at_na)
            min_dist = min(distances_n_r)
            if counter == 0:
                print(min_dist)
                print("ANY SODIUM ADDED!!!!!!!!!")
                return []
    else:
        print("ANY SODIUM ADDED!!!!!!!!!")
        return []


############
# EMPIEZA EL PROTOCOLO
##########

req_num = sys.argv[1]
path = sys.argv[2]
path_req = path + req_num + "/"


path_out = path_req + "HW/"
path_wats = path_req + "waters/"

filog_tree = open(path_req + "filter_rec_list", "r")
ord_list = filog_tree.readlines()

ref = ord_list[0]
print("REFERENCE:", ref)
# print d_ww, "esta es la dist wat-wat minima!!!!"


##########################
# add SODIUM just if it exists in MODEL
##########################

sod_list = []

file_wat_ions = path_out + "wations.pdb"
# first we get the Na from user model
if os.path.isfile(file_wat_ions):
    wations = open(file_wat_ions, "r")
    wations_r = wations.readlines()
    for atom in wations_r:
        molec = atom[17:20].strip()
        if molec == "NA":
            sod_list.append(atom)
            print(atom)


print(len(sod_list), "SODIOS EN LA LISTA")

# open and save MODEL
model = open(path_req + "protein.pdb", "r")
r_model = model.readlines()
prot_atoms = r_model
# lines_in_model = len(r_model))

#### INCLUDE LIGAND!
# include ligand if exists in MODEL

lines_in_ligand = 0
if file_exist(path_out + "ligand_atom.pdb"):
    f_lig = open(path_out + "ligand_atom.pdb", "r")
    f_lir_r = f_lig.readlines()
    for line in f_lir_r:
        r_model.append(line)
    lines_in_ligand = len(f_lir_r)  # num lines ligand to delete at the end

## 1st include sodium
if sod_list != []:
    na_added = add_na_to_model(sod_list, r_model)
    print(na_added)
    if na_added != []:
        r_model.append(na_added)

# then include internal waters
# compute Circular Variance from CBs of protein

CB_atoms_list = []
for line in r_model:
    if line.startswith("ATOM"):
        tipo = line[12:16].strip()
        if tipo == "CB":
            coxA = line[30:38].strip()
            coyA = line[38:46].strip()
            cozA = line[47:54].strip()
            coord_i = np.array([coxA, coyA, cozA], float)
            CB_atoms_list.append(coord_i)

#########################
### XTAL WATERS: first add waters from MODEL
#########################

xtal_wats = []
if os.path.isfile(file_wat_ions):
    wations = open(file_wat_ions, "r")
    wations_r = wations.readlines()
    for atom in wations_r:
        molec = atom[17:20].strip()
        if molec == "HOH" and atom[16] != "B":
            xtal_wats.append(atom)
            ################

refined_wats = []

if xtal_wats != []:
    for wat in xtal_wats:
        r_model.append(wat)  # add to model
        refined_wats.append(wat)  # add to wat list
        #     #if wat_list == []:
        #     wat_list.append(wat)


w_descartes_1 = []

## Add waters from crystals following blast order

for receptor in ord_list:
    rec_name = receptor.strip()
    path_rec_wats = path_wats + rec_name + "_wats.pdb"
    if file_exist(path_rec_wats):
        wat_pdb = open(path_rec_wats, "r")
        wat_r = wat_pdb.readlines()
        for wat in wat_r:
            # resi = wat[23:26]
            # end_line_h2o = wat[79:-1]
            local_rmsd = float(wat[-6:].strip())
            # #print 'HOH', resi, end_line_h2o, resolution
            if local_rmsd <= 2:
                # if 1 == 1: #test
                if add_h2o_to_model(wat, r_model):
                    cv = check_CircVar(CB_atoms_list, wat)
                    if cv > 0.6:
                        r_model.append(wat)  # add to model
                        refined_wats.append(wat)  # add to wat list
                else:
                    w_descartes_1.append(wat)
    else:
        wat_r = []

#  print("From UNIPROT:", receptor[:-1], "there are", len(wat_r), "waters")

# print(len(refined_wats), "TOTAL of water molecules included")

### write to pdb

model_out = open(path_out + "protein-wat.pdb", "w")
for atom_line in r_model:
    # model_out.write(atom_line)
    if atom_line.startswith("ATOM"):
        molec = atom_line[17:20].strip()
        if molec in aa3:
            model_out.write(atom_line)

## append LIGAND(s) in HETATM form

###SI EXISTE LIGANDO!!!!!!!!!!!!!!!!!

if file_exist(path_out + "ligand_het.pdb"):
    f_lig_in = open(path_out + "ligand_het.pdb")
    f_r_lig = f_lig_in.readlines()
    for line in f_r_lig:
        model_out.write(line)

for na_line in r_model:
    if na_line.startswith("HETATM"):
        tipo = na_line[12:16].strip()
        if tipo == "NA":
            new_l = na_line
            new_l = na_line[:7] + "4000" + na_line[11:21] + "  " + "400" + na_line[26:]
            sod = str(new_l)
            model_out.write(sod)

for w in refined_wats:
    model_out.write(w)

model_out.close()

