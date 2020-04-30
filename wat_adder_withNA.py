"""
script to add waters WITH internal NA
usage ==> python wat_adder_withNA.py req_num
Homolwat 4B of 6
__author__: Eduardo Mayol
"""
# import MySQLdb as mdb
import math
import os

from datetime import datetime
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
    global d_wr, d_wn, d_ww
    # tipo_w = h2o_i[17:20].strip()
    todas_dists = []
    polar_dists = []

    coxA = h2o_i[30:38].strip()
    coyA = h2o_i[38:46].strip()
    cozA = h2o_i[47:54].strip()
    coord_wat = (float(coxA), float(coyA), float(cozA))

    # with the h2o coord, chwck if fits in the model
    for line in pdb_model:
        if line.startswith("ATOM"):
            tipo = line[12:16].strip()
            molec = line[17:20].strip()
            if not tipo.startswith("H"):  # avoid Hs in residues
                coxB = line[30:38].strip()
                coyB = line[38:46].strip()
                cozB = line[47:54].strip()
                coord_rsd = (float(coxB), float(coyB), float(cozB))
                dist = Euc_dist(coord_wat, coord_rsd)
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
                coxC = line[30:38].strip()
                coyC = line[38:46].strip()
                cozC = line[47:54].strip()
                coord_na = (float(coxC), float(coyC), float(cozC))
                dist = Euc_dist(coord_wat, coord_na)
                todas_dists.append(dist)
                # avoid clash with NA
                # print dist
                if dist < d_wn:
                    return False
                elif dist < 4:
                    polar_dists.append(dist)
            elif molec == "HOH":
                coxD = line[30:38].strip()
                coyD = line[38:46].strip()
                cozD = line[47:54].strip()
                coord_wat2 = (float(coxD), float(coyD), float(cozD))
                dist = Euc_dist(coord_wat, coord_wat2)
                todas_dists.append(dist)
                # Avoid clash with previous added waters
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


def na_near_d250(na_line, pdb_model):
    """ first we need to compute de COM of the prot
        then place a sodium if is near D250
        returns True or False
    """
    # get CAs from model
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
    # get_coord from NA
    coxZ = na_line[30:38].strip()
    coyZ = na_line[38:46].strip()
    cozZ = na_line[47:54].strip()
    coord_Z = (float(coxZ), float(coyZ), float(cozZ))
    distance_to_COM = Euc_dist(COM_coord, coord_Z)
    if distance_to_COM < 10:
        return True
    else:
        return False


def add_na_to_model(na2check, pdb_model):
    """
    Compute if Na ion fits in model
    returns True or False
    """
    distances_n_r = []
    coxG = na2check[30:38].strip()
    coyG = na2check[38:46].strip()
    cozG = na2check[47:54].strip()
    coordG = (float(coxG), float(coyG), float(cozG))
    for atom2 in pdb_model:
        if atom2.startswith("ATOM"):
            coxF = atom2[30:38].strip()
            coyF = atom2[38:46].strip()
            cozF = atom2[47:54].strip()
            coordF = (float(coxF), float(coyF), float(cozF))
            dist_at_na = Euc_dist(coordG, coordF)
            if dist_at_na > d_nr:
                distances_n_r.append(dist_at_na)
            else:
                return False  # clash
    min_dist = min(distances_n_r)
    if min_dist > 4:
        return False  # too far from rsd
    else:
        return True  # Na fits


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

# se abre el archivo de la proteina y se guarda en r_model
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

# get lsit of CBs coords to compute CV

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


##########################
######### SODIUM ##########
##########################

# check if there's sodium in model AND if it's near D250

sodios_xtal = []
file_wat_ions = path_out + "wations.pdb"
# first we get the Na from user model
if os.path.isfile(file_wat_ions):
    wations = open(file_wat_ions, "r")
    wations_r = wations.readlines()
    for atom in wations_r:
        molec = atom[17:20].strip()
        if molec == "NA":
            sodios_xtal.append(atom)
            r_model.append(atom)
            print(atom)

sod_xtal_d250 = False
# And check if its near D250
if sodios_xtal != []:
    dentro = 0
    for naxtl in sodios_xtal:
        if na_near_d250(naxtl, r_model):
            dentro += 1
    if dentro > 0:
        sod_xtal_d250 = True


# if any we get from other xtals
waterssodium = []
if sod_xtal_d250 == False:
    sod_list = []
    na_file = open(path_wats + "sodiums.pdb", "r")
    na_r = na_file.readlines()
    for na in na_r:
        sod_list.append(na)
    na_file.close()

    contador = 0
    for na_ion in sod_list:
        # print "contador", contador
        if contador == 0:
            if na_near_d250(na_ion, r_model):
                if add_na_to_model(na_ion, r_model):
                    r_model.append(na_ion)
                    contador += 1
                    end_line_na = na_ion[79:-1]
                    # print end_line_na
                    elem1 = end_line_na.split()
                    pdb_code_na = elem1[0]
                    print("pdb sodio:", pdb_code_na)
                    # save water near sodiom from the same xtal
                    f_in_watna = open(path_wats + "wat_near_sodiums.pdb", "r")
                    f_r_watna = f_in_watna.readlines()
                    f_in_watna.close()

                    wats_na = []
                    for watline in f_r_watna:
                        wat_end_line = watline[79:-1]
                        elem2 = wat_end_line.split()
                        pdb_code_wat = elem2[0]
                        if pdb_code_wat == pdb_code_na:
                            # print "pdb sodio:", pdb_code_na
                            # print watline
                            if add_h2o_to_model(watline, r_model):
                                waterssodium.append(watline)
        else:
            break

else:
    print("Na included from model")

print("waters around NA:", len(waterssodium))
for x in waterssodium:
    print(x)


refined_wats = []


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
            cv = check_CircVar(CB_atoms_list, atom)
            if cv > 0.6:
                xtal_wats.append(atom)
        ################
#### MODIF!!!
##SE PILLAN TODAS LAS AGUAS INTERNAS DEL CRISTAL


if xtal_wats != []:
    for wat in xtal_wats:
        r_model.append(wat)  # add to model
        refined_wats.append(wat)  # add to wat list


# se incluyen las aguas q estancerca del NA
if waterssodium != []:
    for wat_na in waterssodium:
        r_model.append(wat_na)  # add to model
        refined_wats.append(wat_na)  # add to wat list


w_descartes_1 = []

## de la lista de receptores se pillan las aguas de la carpeta waters
for receptor in ord_list:
    rec_name = receptor.strip()
    path_rec_wats = path_wats + rec_name + "_wats.pdb"
    if file_exist(path_rec_wats):
        wat_pdb = open(path_rec_wats, "r")
        wat_r = wat_pdb.readlines()
        for wat in wat_r:
            local_rmsd = float(wat[-6:].strip())
            # if 1 == 1:  # test
            if local_rmsd <= 2:
                if add_h2o_to_model(wat, r_model):
                    cv = check_CircVar(CB_atoms_list, wat)
                    if cv > 0.6:
                        r_model.append(wat)  # add to model
                        refined_wats.append(wat)  # add to wat list
                else:
                    w_descartes_1.append(wat)
    else:
        wat_r = []


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

cont = 0

# RENUMBER NA ion (residue number)
for na_line in r_model:
    if na_line.startswith("HETATM"):

        tipo = na_line[12:16].strip()
        if tipo == "NA":
            num_a = str(4000 + cont)
            num_b = str(400 + cont)
            # new_l = na_line
            new_l = na_line[:7] + num_a + na_line[11:21] + "  " + num_b + na_line[26:]
            sod = str(new_l)
            model_out.write(sod)
            cont += 1
# Add waters to model
for w in refined_wats:
    model_out.write(w)

model_out.close()

