#!/bin/bash -l 
start=`date +%s`

# MODIFY PATHS 
path='/home/edu/HomolWat/'
path_scripts='/home/edu/scripts/'
path_rec=${path}ref_gpcr/
path_wats=${path}REC_WAT/
path_jobs=${path}validation/
path_pymol='/opt/soft/pymol-2/bin/pymol' 

prot=$1
pdb_folder=${prot:0:4}_HW 
folder_out=$path_jobs$pdb_folder
# prepare folders and MODEL  1/6
python ${path_scripts}prep_all.py $prot $path $path_jobs $path_rec $pdb_folder $path_scripts
echo 'running pymol part...'
# Superpose crystals to model 2/6
$path_pymol -cqr ${path_scripts}prep_PDBs_mproc.py $pdb_folder $path_jobs $path_wats $path_scripts> pymol_info.txt
# Sort and group waters 3/6
python ${path_scripts}gene_wats_file_refined_v4.py $pdb_folder $path_jobs $path_scripts

# Chose one of the script to add (or not internal Na ion)
# Add internal water molecules 4/6
python ${path_scripts}wat_adder_withNA.py $pdb_folder $path_jobs $path_scripts
#python ${path_scripts}wat_adder_noNA.py $pdb_folder $path_jobs $path_scripts

# merge MODEL and waters, renumber waters and prepare output 5/6
python ${path_script_EM}merge_prot_waters.py ${folder_out}/HW/ ${pdb_folder}
# prepare pymol sesion with output 6/6
$path_pymol -cqr ${path_scripts}gene_final_pse.py $pdb_folder $path_jobs $pdb_folder
end=`date +%s`
runtime=$((end-start))
echo "execution time:" $runtime

# Delete files 
rm ${folder_out}/REC/*.pdb

echo "output is located at "${folder_out} "AS" ${pdb_folder}_HW.pdb

