#!/usr/bin/env python
import psi4
import parmed as pm
from ..defaults import *

def OptimizePDB(pdbfile,charge=0,mult=1,method="scf/6-31g**"):
    psi4.set_memory(AVAILABLE_PSI4_MEMORY)
    xyzdata = f"{charge} {mult}"
    mol = pm.load_file(pdbfile)
    mol.write_pdb(pdbfile+".ORIG")
    for atom in mol.atoms:
        xyzdata += f"\n{atom.element_name:<3} {atom.xx:>10.3f} {atom.xy:>10.3f} {atom.xz:>10.3f}"
    _psi4_molecule = psi4.geometry(xyzdata)
    print("In the OptimizePDB function now")
    psi4.optimize(method, molecule=_psi4_molecule)
    for i,atom in enumerate(mol.atoms):
        atom.xx = _psi4_molecule.x(i)
        atom.xy = _psi4_molecule.y(i)
        atom.xz = _psi4_molecule.z(i)
    mol.write_pdb(pdbfile)
    del mol
    return 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--pdbfile",dest="pdbfile",required=True)
    parser.add_argument("-c","--charge",dest="charge",default=0)
    parser.add_argument("-s","--spin",dest="spin",default=1)
    parser.add_argument("-m","--method",dest="method",default="scf/6-31g**")
    args = parser.parse_args()
    OptimizePDB(args.pdbfile,args.charge,args.spin,args.method)