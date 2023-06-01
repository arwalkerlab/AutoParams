#!/usr/bin/env python
import shutil
import argparse
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Application.utilities.utilities import *
from Application.utilities.main_job import *

if __name__ == "__main__":
    # Initialization
    startdir = os.path.abspath(os.getcwd())
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file",dest="file",help="/path/to/3D/structure",required=True)
    parser.add_argument("-c","--charge",dest="charge",help="integer charge of structure",default=0)
    parser.add_argument("-s","--spin",dest="spin",help="integer spin multiplicity of structure",default=1)
    parser.add_argument("--method",dest="method",help="Method to use for calculations",default="b3lyp")    
    parser.add_argument("--basis",dest="basis",help="Basis set to use for calculations",default="6-31g**")
    parser.add_argument("--headconnect",dest="headconnect",help="Atom name for head connection to previous residue",default="0")
    parser.add_argument("--tailconnect",dest="tailconnect",help="Atom name for tail connection to next residue",default="0")
    parser.add_argument("--optimize-first",dest="optimize",help="Runs unrestrained geometry optimization before parametrization.",action="store_false")
    parser.add_argument("--restype",dest="restype",help="structure type (DNA, RNA, Protein)",default="RNA")
    parser.add_argument("-v","--verbose",dest="verbose",help="Enable console logging.",action="store_true")

    # Parse and process CL arguments
    args = parser.parse_args()
    args.charge = int(args.charge)
    args.spin = int(args.spin)
    connections=[args.headconnect,args.tailconnect]

    # Generate job folder
    job_id = random_job_identifier()
    if args.verbose:
        print("Job ID:",job_id)

    # Set up MainJob structure
    current_job = MainJob(job_id,
                          int(args.charge),
                          int(args.spin),
                          args.restype,
                          args.optimize,
                          True,  ## Forcibly skipping DB checks in batch mode.
                          args.method,
                          args.basis,
                          connections)
    # Copy PDB file into job folder and enter job folder
    shutil.copy(os.path.abspath(args.file),current_job._job_folder)
    if args.verbose:
        print("File name:",os.path.basename(args.file))

    current_job.file_list["Original PDB"] = os.path.join(current_job._job_folder,os.path.basename(args.file))
    current_job.file_list["FRCMOD"] = current_job.file_list["Original PDB"].replace(".pdb",".frcmod")
    current_job.file_list["MOL2"] = current_job.file_list["Original PDB"].replace(".pdb",".mol2")
    current_job.file_list["PRMTOP"] = current_job.file_list["Original PDB"].replace(".pdb",".prmtop")
    current_job.file_list["INPCRD"] = current_job.file_list["Original PDB"].replace(".pdb",".inpcrd")
    current_job.file_list["ChemDraw"] = current_job.file_list["Original PDB"].replace(".pdb",".png")
    current_job.file_list["Working PDB"] = current_job.file_list["Original PDB"]
    current_job.file_list["SMILES"] = current_job.file_list["Original PDB"].replace(".pdb",".png").replace(".pdb",".smi")
    current_job.file_list["LeapLog"] = os.path.join(current_job._job_folder, "leap.log")
    
    # Optimize structure (checks optimization flag internally)
    if not current_job.Optimize():
        if args.verbose:
            print("Unable to optimize structure:", os.path.basename(args.file))
        sys.exit()

    # Generate RESP charges.
    if not current_job.RESPCharges():
        if args.verbose:
            print("Unable to generate RESP charges:", os.path.basename(args.file))
        sys.exit()

    # Generate MOL2 and FRCMOD
    if not current_job.Parametrize():
        if args.verbose:
            print("Unable to generate parameters:", os.path.basename(args.file))
        sys.exit()

    # Test resulting parameters.
    if not current_job.TestParams():
        if args.verbose:
            print("Unable to generate MD inputs as test files:", os.path.basename(args.file))
        sys.exit()

    # Copy parameters to original submission folder.
    shutil.copy(current_job.file_list["FRCMOD"],startdir)
    shutil.copy(current_job.file_list["MOL2"],startdir)

    # Generate Database Folder
    if not current_job.AddResultsToDB():
        if args.verbose:
            print("Unable to add parameters to database:", os.path.basename(args.file))
        sys.exit()

