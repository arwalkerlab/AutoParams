#!/usr/bin/env python 
import psiresp
import qcelemental as qcel
import parmed as pm
import psi4
import sys
sys.path.append("/app/")
from utilities.defaults import AVAILABLE_PSI4_MEMORY
psi4.set_memory(AVAILABLE_PSI4_MEMORY)


def PDB_to_XYZ_data(pdbfile,charge,mult):
    xyzdata = f"{charge} {mult}"
    mol = pm.load_file(pdbfile)
    for atom in mol.atoms:
        xyzdata += f"\n{atom.element_name:<3} {atom.xx:>10.3f} {atom.xy:>10.3f} {atom.xz:>10.3f}"
    return xyzdata

def GetCapAtomIndices(pdbfile,cap_atoms):
    indices = []
    mol = pm.load_file(pdbfile)
    for atom in mol.atoms:
        if atom.name in cap_atoms:
            indices.append(atom.idx)
    return indices

def RunJob(pdbfile,charge,mult,workdir,cap_atoms):
    print("Running Job...")
    qc_molecule = qcel.models.Molecule.from_data(PDB_to_XYZ_data(pdbfile,charge,mult))
    print("qc_molecule generated")
    psi_molecule = psiresp.Molecule(qcmol=qc_molecule,charge=charge,multiplicity=mult)
    print("psi_molecule_generated")
    psi_molecule.add_conformer(qcmol=qc_molecule)
    cap_indices = []
    if cap_atoms != []:
        print("Adding cap atom charge restraints")
        constraints = psiresp.ChargeConstraintOptions()
        print("added contraints to system")
        cap_indices = GetCapAtomIndices(pdbfile,cap_atoms)
        constraints.add_charge_sum_constraint_for_molecule(psi_molecule, charge=0, indices=cap_indices)
        print("added charge constraint to specified atoms.")

    psi_job = psiresp.Job(molecules=[psi_molecule],working_directory=workdir)
    print("Running psi4 job...")
    psi_job.run() ## MUST HAVE BOTH INSTANCES OF psi_job.run() TO FUNCTION CORRECTLY.
    print("Obtaining RESP charges from psi4 job...")
    respcharges = psi_job.run()
    with open("resp.out","w") as f:
        for i,ch in enumerate(respcharges[0]):
            if i in cap_atoms:
                continue
            f.write(f"{ch}\n")
    print("Job complete!")

def PrintToConsole(args):
    print("#" * 60)
    print("PDBFile:        ",args.pdbfile)
    print("Charge:         ",args.charge)
    print("Multiplicity:   ",args.multiplicity)
    print("Work Directory: ",args.workdir)
    print("Cap Atoms:      ",args.cap_atoms)
    print("#" * 60)



if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-p",dest="pdbfile",required=True)
    parser.add_argument("-c",dest="charge",default=0)
    parser.add_argument("-m",dest="multiplicity",default=1)
    parser.add_argument("-d",dest="workdir",default="./")
    parser.add_argument("-z",dest="cap_atoms",default=[],nargs="+")
    args = parser.parse_args()
    
    PrintToConsole(args)

    RunJob(args.pdbfile,args.charge,args.multiplicity,args.workdir,args.cap_atoms)
