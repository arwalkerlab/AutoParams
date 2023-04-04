### from .optimizer import *
import psi4
import parmed as pm

def OptimizePDB(pdbfile,charge=0,mult=1,method="scf/6-31g**"):
    xyzdata = f"{charge} {mult}"
    mol = pm.load_file(pdbfile)
    mol.write_pdb(pdbfile+".ORIG")
    for atom in mol.atoms:
        xyzdata += f"\n{atom.element_name:<3} {atom.xx:>10.3f} {atom.xy:>10.3f} {atom.xz:>10.3f}"
    _psi4_molecule = psi4.geometry(xyzdata)
    psi4.optimize(method, molecule=_psi4_molecule)
    for i,atom in enumerate(mol.atoms):
        atom.xx = mol.x(i)
        atom.xy = mol.y(i)
        atom.xz = mol.z(i)
    mol.write_pdb(pdbfile)
    

