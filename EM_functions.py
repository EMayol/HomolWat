"""
functions from Eduardo Mayol
to use in HomolWat
"""

import numpy as np
import math
from subprocess import check_output

## insert path were blast+ is installed
path_blast = "/path/to/blast/" + "bin/"

aa = {
    "ALA": "A",
    "ARG": "R",
    "LEU": "L",
    "MET": "M",
    "LYS": "K",
    "GLN": "Q",
    "GLU": "E",
    "GLH": "E",
    "ILE": "I",
    "TRP": "W",
    "SER": "S",
    "TYR": "Y",
    "PHE": "F",
    "VAL": "V",
    "HIS": "H",
    "HIE": "H",
    "HID": "H",
    "HIP": "H",
    "ASN": "N",
    "THR": "T",
    "CYS": "C",
    "CYX": "C",
    "CYP": "C",
    "ASP": "D",
    "ASH": "D",
    "GLY": "G",
    "PRO": "P",
    "ACE": ".",
    "NME": ".",
    "WAT": ".",
    "SOL": ".",
    "HOH": ".",
    "MOL": "*",
}


def pdb2fasta_save(pdb_file, path2save):
    # open file to save and wite header
    fout = open(path2save + "prot.fasta", "w")
    header = "> " + path2save + "prot" + "\n"
    fout.write(header)
    # convert pdb to fasta and save
    reso = 0
    count = -1
    for line in pdb_file:
        if line.startswith("ATOM"):
            tipo = line[12:16].strip()
            if tipo == "CA":
                resn = str(line[17:20])
                resi = int(line[23:26])
                if resi != reso:
                    reso = resi
                    count += 1
                    if count == 72:
                        fout.write("\n")
                        count = 0
                fout.write(aa[resn])
    fout.close()


def pdb2fasta_chains(pdb_file, pdbname, chain, path2save):
    # open file to save and wite header
    fout = open(path2save + "prot_" + pdbname + "_ch" + chain + ".fasta", "w")
    header = "> " + path2save + "prot_" + pdbname + "_ch" + chain + "\n"
    fout.write(header)
    # convert pdb to fasta and save
    reso = 0
    count = -1
    for line in pdb_file:
        if line.startswith("ATOM"):
            tipo = line[12:16].strip()
            ch = line[21:22].strip()
            if chain == ch:
                if tipo == "CA":
                    resn = str(line[17:20])
                    resi = int(line[23:26])
                    if resi != reso:
                        reso = resi
                        count += 1
                        if count == 72:
                            fout.write("\n")
                            count = 0
                    fout.write(aa[resn])
    fout.close()


def center_of_molec(coord_list):
    x_list, y_list, z_list = [], [], []
    for coord in coord_list:
        # print coord
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


def Euc_dist(p, q):
    dist = math.sqrt((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2)
    return dist


def check_CircVar(CB_coords, at_coords, rcv):
    # info_CVs = []
    coxB = at_coords[0]
    coyB = at_coords[1]
    cozB = at_coords[2]
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


def run_blastp(fasta_in, path_scripts, path_out):
    ### ojo con esto (modif y borrar)
    # path_db = "/var/www/HW/hw/scripts/"
    path_db = path_scripts
    f_in1 = open(path_db + "Pdbs_resol_INACT", "r")
    f_r1 = f_in1.readlines()
    unipsWithwat = []
    info_unips = {}
    for line1 in f_r1:
        elem0 = line1.split("_")
        unip = elem0[0]
        pdbs = elem0[1].strip()
        if len(pdbs) > 0:
            info_unips[unip] = pdbs
            unipsWithwat.append(unip)
    blast1 = check_output(
        path_blast
        + "blastp -db "
        + path_db  # path were files ref_gpcr.phr, ref_gpcr.pin and ref_gpcr.psq are located
        + "ref_gpcr -query "
        + fasta_in
        + ' -outfmt "10 sseqid score pident qcovs evalue" -evalue 100 ',
        shell=True,
    ).split("\n")
    rec_sorted_list = []
    tmp_rec_list = []
    for rec in blast1:
        if len(rec) > 1:
            print(rec)
            # print type(rec), len(rec)
            elem = rec.split(",")
            unip = elem[0]
            #  score = elem[1].strip()
            # print unip, pident
            if unip not in tmp_rec_list:
                rec_sorted_list.append(rec)
                tmp_rec_list.append(unip)
    # save output
    f_out = open(path_out + "rec_sort_list", "w")
    f_out2 = open(path_out + "rec_sort_list.csv", "w")
    f_filt_out = open(path_out + "filter_rec_list", "w")
    f_out2.write("Uniprot,blastp score,pident,coverage,evalue,PDB codes\n")
    recep_score = []
    for rec2 in rec_sorted_list:
        element = rec2.split(",")
        unip, score, pident, coverage, evalue = (
            element[0],
            element[1],
            element[2],
            element[3],
            element[4].strip(),
        )
        info = (unip, score, pident, coverage, evalue)
        recep_score.append(info)

    # sort blast output
    blast_recs = []
    for receptor in recep_score:
        blast_recs.append(receptor[0])
        if receptor[0] in unipsWithwat:
            info_out = receptor[0] + " " + str(receptor[1])
            info_out2 = (
                receptor[0]
                + ","
                + str(receptor[1])
                + ","
                + str(receptor[2])
                + ","
                + str(receptor[3])
                + ","
                + str(receptor[4])
                + ","
                + str(info_unips[receptor[0]])
            )
            f_filt_out.write(receptor[0] + "\n")
            f_out.write(info_out + "\n")
            f_out2.write(info_out2 + "\n")
    # in case there's no output write score 0
    for uniprot in unipsWithwat:
        if uniprot not in blast_recs:
            info_out_b = uniprot + " 0"
            info_out2_b = uniprot + ",0,0,0,>100," + str(info_unips[uniprot])
            f_filt_out.write(uniprot + "\n")
            f_out.write(info_out_b + "\n")
            f_out2.write(info_out2_b + "\n")


def run_blastp_chains(fasta_in, path_scripts, path_out):
    path_db = path_scripts
    f_in1 = open(path_db + "Pdbs_resol_INACT", "r")
    f_r1 = f_in1.readlines()
    unipsWithwat = []
    info_unips = {}
    for line1 in f_r1:
        elem0 = line1.split("_")
        unip = elem0[0]
        pdbs = elem0[1].strip()
        if len(pdbs) > 0:
            info_unips[unip] = pdbs
            unipsWithwat.append(unip)
    # run blast
    blast2 = check_output(
        path_blast
        + "blastp -db "
        + path_db  # path were files ref_gpcr.phr, ref_gpcr.pin and ref_gpcr.psq are located
        + "ref_gpcr -query "
        + fasta_in
        + ' -outfmt "10 sseqid evalue" -evalue 100 ',
        shell=True,
    ).split("\n")
    return blast2
