#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 11:16:06 2025

@author: malab
"""

import os
import subprocess


def range1(start, end):
    return range(start, end+1)


# Parameters
path_to_pepsi = '/home/malab/Desktop/Pepsi-SAXS-Linux'
pdb_range = range1(0, 640)
output_directory = '/home/malab/iBME/pepsi_output'
path_to_pdbs = '/home/malab/Desktop/PBD_ID/Structures'

# Counter
counter = 0

print("Current working directory: {0}".format(os.getcwd()))
os.chdir('/home/malab/Desktop/Pepsi-SAXS-Linux')
print("Current working directory: {0}".format(os.getcwd()))

# Pepsi Loop
for i in pdb_range:
    pdb_name = f"mm016_{i:03}.pdb"  # Name of individual pdb file
    pdb_path = os.path.join(path_to_pdbs, pdb_name)  # Create full path to pdb

    result = subprocess.run(["{}/Pepsi-SAXS".format(path_to_pepsi), pdb_path], 
                            cwd=output_directory, capture_output=True, text=True)

    # Delete log file (optional)
    log_pepsi = f"mm016_{i:03}.log"
    os.remove(os.path.join(output_directory, log_pepsi))
    
    #Counter
    counter += 1
    
    #Error
    if result.returncode != 0: #error
        print(f" crysol failed for {result}")
        print(result.stderr)
        continue

print(counter)
