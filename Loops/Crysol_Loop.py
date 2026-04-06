#Crysol Loop

import subprocess
import os
import pandas

#fix range to include all files specified
def range1(start, end):
    return range(start, end+1)

##Parameters
path_to_pbds = "/home/malab/Desktop/PBD_ID/Structures" #path to folder with pdbs
pdb_file = range1(0, 952) #how many pdb files you want to sim and their ids
smax = 1 #maximum q value
output_directory = "/home/malab/Desktop/BioEn-master/examples/scattering/files/raw_crysol_data/raw_abs" #path to output folder

#Setting up terminal- FOR LINUX RUN ONLY
env = os.environ.copy()
env['ATSAS'] = '/home/malab/ATSAS-3.2.1-1'
env['PATH'] = f"{env['ATSAS']}/bin:{env['PATH']}"

#Run Counter
counter = 0

##Crysol Loop
for i in pdb_file:
    pdb_name = f"mm016_{i:03}.pdb" #Name of individual pdb file
    pdb_path = os.path.join(path_to_pbds, pdb_name)  # Create full path to pdb
    #crysol_output = os.path.join(output_directory, pdb_name[:-4])  # Create the output to be 'path' + file name without the 'pdb'
    
    result = subprocess.run(["crysol", pdb_path, "--smax={}".format(smax)],
    env=env, cwd=output_directory, capture_output=True, text=True) #run crysol for each pdb file with specific smax and save to output directory

    ##Crysol outputs 4 files- for which file(s) you wish to keep, put a comment. Anything without comment will get deleted
    int_delete = f"mm016_{i:03}.int"
    #abs_delete = f"mm016_{i:03}.abs"
    log_delete = f"mm016_{i:03}.log"
    alm_delete = f"mm016_{i:03}.alm"


    #Place comment here as well following above rules
    os.remove(os.path.join(output_directory, int_delete))
    #os.remove(os.path.join(output_directory, abs_delete))
    os.remove(os.path.join(output_directory, log_delete))
    os.remove(os.path.join(output_directory, alm_delete))


    if result.returncode != 0: #error
        print(f" crysol failed for {result}")
        print(result.stderr)
        continue
    
    counter += 1

print(counter)







