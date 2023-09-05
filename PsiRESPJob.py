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
    qc_molecule = qcel.models.Molecule.from_data(PDB_to_XYZ_data(pdbfile,charge,mult))
    psi_molecule = psiresp.Molecule(qcmol=qc_molecule,charge=charge,multiplicity=mult)
    psi_molecule.add_conformer(qcmol=qc_molecule)
    if cap_atoms != []:
        constraints = psiresp.ChargeConstraintOptions()
        constraints.add_charge_sum_constraint_for_molecule(psi_molecule, charge=0, indices=GetCapAtomIndices(pdbfile,cap_atoms))

    psi_job = psiresp.Job(molecules=[psi_molecule],working_directory=workdir)
    psi_job.run() ## MUST HAVE BOTH INSTANCES OF psi_job.run() TO FUNCTION CORRECTLY.
    respcharges = psi_job.run()
    with open("resp.out","w") as f:
        for ch in respcharges[0]:
            f.write(f"{ch}\n")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-p",dest="pdbfile",required=True)
    parser.add_argument("-c",dest="charge",default=0)
    parser.add_argument("-m",dest="multiplicity",default=1)
    parser.add_argument("-d",dest="workdir",default="./")
    parser.add_argument("-z",dest="cap_atoms",default=[],nargs="+")
    args = parser.parse_args()

    RunJob(args.pdbfile,args.charge,args.multiplicity,args.workdir,args.cap_atoms)
