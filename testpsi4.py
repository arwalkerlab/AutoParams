### from .optimizer import *
import psi4
import parmed as pm
pdbfile = "test.pdb"

def PDBtoXYZ(pdbfile,charge=0,mult=1):
    xyzfile = f"{charge} {mult}"
    mol = pm.load_file(pdbfile)
    # xyzfile+= f"{len(mol.atoms)}\n"
    for atom in mol.atoms:
        xyzfile += f"\n{atom.element_name:<3} {atom.xx:>10.3f} {atom.xy:>10.3f} {atom.xz:>10.3f}"
    return xyzfile

xyzfile = PDBtoXYZ(pdbfile,0,1)
# with open("tmp.xyz","w") as f:
#     f.write(xyzfile)



psi4.set_memory('1024 MB')
_psi4_molecule = psi4.geometry(xyzfile)

psi4.set_options({'reference': 'rhf'})
psi4.optimize('scf/6-31g**', molecule=_psi4_molecule)

