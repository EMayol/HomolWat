###############################################################################
"""
This script groups water molecules from crystals in families
according to resolution and Bfactor
It also saves Na ion near D2.50 position and water molecules around it
in case it is necessary to incorporate it the final model
wether wat_adder_withNA.py or wat_adder_noNA.py is used
Homolwat 3 of 6
__author__: Eduardo Mayol
"""
###############################################################################

import os
import sys
import glob
import numpy as np
import math
from EM_functions import Euc_dist


def get_coords(atom_line):
    coxA = atom_line[30:38].strip()
    coyA = atom_line[38:46].strip()
    cozA = atom_line[47:54].strip()
    coord_A = np.array([coxA, coyA, cozA], float)
    return coord_A


def get_wats_filtered(letter, pdb_id, recep):
    wats_all, wats_good, sodiums = [], [], []
    f_in = open(path_rw + pdb_id + letter + recep + ".pdb", "r")
    f_r = f_in.readlines()
    xtal_id = recep + "-" + pdb_id
    for line in f_r:
        if line.startswith("HETATM"):
            res_n = line[17:21].strip()
            Bfactor = float((line[60:66].strip()))
            new_line = ""
            p1 = line[:79]
            new_line = p1 + xtal_id + letter[:2] + "\n"
            if res_n == "HOH":  # or res_n == "NA":
                wats_all.append(new_line)
                if Bfactor < 45.0:
                    wats_good.append(new_line)
            elif res_n == "NA":
                if Bfactor < 45.0:
                    sodiums.append(new_line)
    wat_near_na = []
    if sodiums != []:
        coorNA = get_coords(sodiums[0])
        for wat in wats_good:
            coord_wat = get_coords(wat)
            if Euc_dist(coorNA, coord_wat) < 3:
                wat_near_na.append(wat)
    return wats_all, wats_good, sodiums, wat_near_na


def file_exist(path):
    if os.path.isfile(path):
        p = True
    else:
        p = False
    return p


def wat_Not_repeated(w_coord, listofwat_coords):
    distances = []
    for element in listofwat_coords:
        dist = Euc_dist(element, w_coord)
        distances.append(dist)
    if min(distances) > 0.1:
        return True
    else:
        return False


# Save paths and open files
req_num = sys.argv[1]
path = sys.argv[2]
path_req = path + req_num + "/"

path_script = path[:-8] + "scripts/"
path_out = path_req + "waters/"
path_rw = path_req + "REC/"

folder = glob.glob(path_rw + "*.pdb")

Recs_resol = []
f_in = open(path_req + "Pdbs_resol_L", "r")
f_r = f_in.readlines()
f_in.close()

f_in2 = open(path_req + "filter_rec_list", "r")
f_order = f_in2.readlines()
f_in2.close()

for rec in f_order:
    print(rec)
    for line in f_r:
        elem = line.split("_")
        recep = elem[0]
        if recep == rec.strip():
            pdbs = elem[1].strip()
            if len(pdbs) >= 1:
                Recs_resol.append([recep, pdbs.strip(" ")])

# Incoporates info from origin of water and local RMSD
Resol_pdbs = {}
f_resol_in = open(path_req + "info_rmsds.txt", "r")
f_resol = f_resol_in.readlines()
f_resol_in.close()
for line in f_resol:
    elem = line.split()
    pdbcode = elem[0]
    rmsd_pdb = elem[1]
    Resol_pdbs[pdbcode] = rmsd_pdb

# Save Na ion and waters around
sodiums = []
w_prior_na = []
for receptor in Recs_resol:
    rec = receptor[0]
    pdbs = receptor[1].split(" ")
    f_out = open(path_out + rec + "_wats.pdb", "w")
    for elem in pdbs:
        wats_good, wats_good_coords = [], []
        for repli in ["A", "B", "C", "D"]:
            repli_let = "_" + repli + "_"
            pdbID = elem + repli_let + rec
            if file_exist(path_rw + pdbID + ".pdb"):
                # open file of rmsd values for each water and save in dictionary
                Rmsd_water = {}
                water_rmsds_file = open(
                    path_rw + "info_rmsds" + pdbID + "_wats.txt", "r"
                )
                water_rmsds = water_rmsds_file.readlines()
                water_rmsds_file.close()
                for line_rms in water_rmsds:
                    column = line_rms.split()
                    water_num = column[2]
                    rmsd = column[3].strip()
                    Rmsd_water[water_num] = rmsd
                # get waters, sodiums and waters near sodiums
                wat_ions = get_wats_filtered(repli_let, elem, rec)
                w_all, w_good, sods, w_near_sod = (
                    wat_ions[0],
                    wat_ions[1],
                    wat_ions[2],
                    wat_ions[3],
                )
                for w2 in w_good:
                    # get info of water and origin
                    watnum = w2[22:27].strip()
                    infoXtalPdb = w2[79:-1]
                    water_chain = w2[21:22].strip()
                    if water_chain == repli:
                        # openfile of refined wat and save to list
                        water_ref_file = open(
                            path_rw + "water_" + watnum + "_" + pdbID + ".pdb", "r"
                        )
                        water_line = water_ref_file.readlines()
                        water_ref_file.close()
                        new_wat_line = (
                            water_line[0][:-1]
                            + infoXtalPdb
                            + " "
                            + Rmsd_water[watnum]
                            + "\n"
                        )
                        # Filter to avoid repeated waters, same posit at 0.1A
                        water_coords = get_coords(new_wat_line)
                        if wats_good == []:
                            wats_good.append(new_wat_line)
                            wats_good_coords.append(water_coords)
                        else:
                            if wat_Not_repeated(water_coords, wats_good_coords):
                                wats_good.append(new_wat_line)
                                wats_good_coords.append(water_coords)
                for na1 in sods:
                    sodiums.append(
                        na1[:-1] + "  " + Resol_pdbs[elem + repli_let[:-1]] + "\n"
                    )
                for w3 in w_near_sod:
                    w_prior_na.append(
                        w3[:-1] + "  " + Resol_pdbs[elem + repli_let[:-1]] + "\n"
                    )
        # sort by Bfactor
        temp_w, sorted_w = [], []
        for wat in wats_good:
            p1 = wat[0:60]
            Bfactor = float(wat[60:66].strip())
            p2 = wat[66:]
            t_w = (p1, Bfactor, p2)
            temp_w.append(t_w)
            # print wat
        sorted_w = sorted(temp_w, key=lambda bf: bf[1])
        for b in sorted_w:
            if b[1] < 10.00:
                l = b[0] + "  " + str(("%.2f" % b[1])) + b[2]
            else:
                l = b[0] + " " + str(("%.2f" % b[1])) + b[2]
            f_out.write(l)
    # prnt ln(wats_all), len(wats_good)
    f_out.close()
## sort sodiums by Bfactor
f_out_na = open(path_out + "sodiums.pdb", "w")
temp_na, sorted_na = [], []
# print "total sodiums:", len(sodiums)
for na in sodiums:
    f_out_na.write(na)
f_out_na.close()

f_out_watna = open(path_out + "wat_near_sodiums.pdb", "w")
for watna in w_prior_na:
    f_out_watna.write(watna)
f_out_watna.close()
