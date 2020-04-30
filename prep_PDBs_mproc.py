###############################################################################
"""

This script is superpose all gpcr crystals with waters to a reference MODEL (Global Alignment)
Then, for each water molecule, a local alignmet is done to refine the position
Alignment are done using PyMol

This is the most cpu demanding process
Meant to use 8 CPUS, adapt it for your resources!

It also generates many files that will be removed at the end of HomolWat process
(up to 110MB in ~5000 files, in folder /REC/ )

IMPRTANT: To be run with "pymol -cqr prep_PDBs.py file.pdb"

Internal water molecules from crystals are taken from HomolWatDB
Homolwat 2 of 6
__author__: Eduardo Mayol
"""
###############################################################################

# import MySQLdb as mdb
from pymol import cmd
import os, sys
import glob
from multiprocessing import Pool
from datetime import datetime


def chunklist(pdb_paths, parts):
    """divide la lista en N listas de similar length"""
    avg = len(pdb_paths) / float(parts)
    out_lists = []
    last = 0.0
    while last < len(pdb_paths):
        out_lists.append(pdb_paths[int(last) : int(last + avg)])
        last += avg
    return out_lists


def superimpose(listapdbs):
    """ make global and local alignments for each structure in listapdbs"""
    path_req = listapdbs.pop()
    path_recWat = listapdbs.pop()
    path_out = path_req + "/REC/"
    # Save information of Global RMSD
    f_out_rmsd = open(path_req + "info_rmsds.txt", "a+")

    allres = "ALA+ARG+ASP+GLU+THR+TYR+HIS+LEU+ILE+VAL+PRO+GLY+TRP+LYS+PHE+ASN+GLN+MET+CYS+SER+ACE+NH2+DI7+DI8"
    wat_ions = "HOH+NA"
    watPDBs = []
    wat_pdbs_file = open(path_req + "Pdbs_resol_L", "r")
    f_rw = wat_pdbs_file.readlines()
    for line in f_rw:
        recep = line.split("_")
        if recep[1].strip() != []:
            pdbs = recep[1].split(" ")
            for pdb in pdbs:
                if pdb != "" and pdb != "\n":
                    watPDBs.append(pdb.strip())
    allowed_uniprots = []
    filog_tree = open(path_req + "filter_rec_list", "r")
    ord_list = filog_tree.readlines()
    for unip in ord_list:
        allowed_uniprots.append(unip.strip())

    # start PyMol process
    cmd.do("cd " + path_req)
    cmd.do("load protein.pdb")

    for pdb in sorted(listapdbs):
        # FILTER CRYSTALS WITH INTERNAL WATERS (Pdbs_resol_L)
        pdbname = pdb[len(path_recWat) : len(path_recWat) + 4]
        unipcode = pdb[len(path_recWat) + 7 : -4]
        pdbID = pdb[len(path_recWat) : -4]
        if pdbname in watPDBs and unipcode in allowed_uniprots:
            pdbname = pdb[len(path_recWat) : len(path_recWat) + 4]
            pdbID = pdb[len(path_recWat) : -4]
            waters = []
            f_in = open(pdb, "r")
            f_r = f_in.readlines()
            for line in f_r:
                if line.startswith("HETATM"):
                    res_n = line[17:21].strip()
                    if res_n == "HOH":  # or res_n == "NA":
                        waters.append(line)
            # Load crystal PDB; do Super AND Align; keep best
            cmd.do("load " + pdb)
            cmd.do("select prot, resn " + allres + " and " + pdbID)
            cmd.do("select close, (br. prot around 5) and " + pdbID)
            cmd.do(
                "select far, not (byres (close or prot) around 5) and not (close or prot) and "
                + pdbID
            )
            cmd.do("remove far")
            rms_data_super = cmd.super(pdbID + "////CA", "protein////CA")
            rms_data_align = cmd.align(pdbID, "protein")
            rms_super = rms_data_super[0]
            rms_align = rms_data_align[0]
            pdb_chain = pdb[len(path_recWat) : len(path_recWat) + 6]
            chain = pdb_chain[-1]
            if rms_super > rms_align:
                rmsd_val = "%.3f" % rms_align
                info = pdb_chain + " " + str(rmsd_val) + " ALIGN\n"
            else:
                rmsd_val = "%.3f" % rms_super
                info = pdb_chain + " " + str(rmsd_val) + " " + " SUPER\n"
            f_out_rmsd.write(info)
            if rms_align > rms_super:
                rms_data_super = cmd.super(pdbID, "protein")
            cmd.do("select not resn " + allres + "+" + wat_ions + " and " + pdbID)
            cmd.do("cd " + path_out)
            cmd.do("save " + pdbID + ".pdb" + ", " + pdbID)
            f_out_rmsd_wats = open(path_out + "info_rmsds" + pdbID + "_wats.txt", "w")
            # Start Local alignments and save RMSD value
            for wat in waters:
                watnum = wat[22:27].strip()
                chain_wat = wat[21:22].strip()
                if chain_wat == chain:
                    cmd.do(
                        "select wat_env, byres "
                        + pdbID
                        + " within 10 of (resi "
                        + watnum
                        + " and resn HOH and "
                        + pdbID
                        + ")"
                    )
                    cmd.create("wat_out", "wat_env")
                    cmd.do("save " + pdbID + "_wat_" + watnum + ".pdb, wat_out")
                    cmd.do("load " + pdbID + "_wat_" + watnum + ".pdb")
                    cmd.do(
                        "select wat2ref, resi "
                        + watnum
                        + " and "
                        + pdbID
                        + "_wat_"
                        + watnum
                    )
                    cmd.do("select template, byres wat2ref around 10 and protein")
                    rmsd_wat_super = cmd.super(
                        pdbID + "_wat_" + watnum + " and name CA",
                        "template and name CA",
                    )
                    rmsd_wat = rmsd_wat_super[0]
                    value_rmsd = "%.3f" % rmsd_wat
                    info_wat = pdbID + " HOH " + watnum + " " + str(value_rmsd) + "\n"
                    f_out_rmsd_wats.write(info_wat)
                    cmd.create("watrefined", "wat2ref")
                    cmd.do("save water_" + watnum + "_" + pdbID + ".pdb, watrefined")


## Determinar variables fijas, paths, etc
req_num = sys.argv[3]  # name of specific folder of the process
path = sys.argv[4]  # path where the above folder is located

path_req = path + req_num + "/"
path_recWat = (
    path[:-8] + "REC_WAT/"
)  # path where crystals with internal waters are located


startTime = datetime.now()


# Specify number of CPUs used --> if this value changes
# you have to modify  the lines below in accordance
num_cpus = 8

folder_RW = glob.glob(path_recWat + "*.pdb")

# print len(folder_RW), "total!"

# if there are >10 pdbs in folder_RW, divide list in N (num_cpus)
# path_req is included in list as pool just accept 1 argument
# inside superimpose function it is removed from list with .pop


if len(folder_RW) > 10:
    div_pathsRW = chunklist(folder_RW, num_cpus)
    for partial_list in div_pathsRW:
        partial_list.append(path_recWat)
        partial_list.append(path_req)
    pathsRW1, pathsRW2, pathsRW3, pathsRW4 = (
        div_pathsRW[0],
        div_pathsRW[1],
        div_pathsRW[2],
        div_pathsRW[3],
    )
    pathsRW5, pathsRW6, pathsRW7, pathsRW8 = (
        div_pathsRW[4],
        div_pathsRW[5],
        div_pathsRW[6],
        div_pathsRW[7],
    )
    # execute in multiprocess
    pool = Pool(processes=num_cpus)
    pool.map(
        superimpose,
        [
            pathsRW1,
            pathsRW2,
            pathsRW3,
            pathsRW4,
            pathsRW5,
            pathsRW6,
            pathsRW7,
            pathsRW8,
        ],
    )
else:
    folder_RW.append(path_recWat)
    folder_RW.append(path_req)
    superimpose(folder_RW)

endTime = datetime.now()
print("execution time:", endTime - startTime)

# Finish PyMol process
cmd.do("quit")

