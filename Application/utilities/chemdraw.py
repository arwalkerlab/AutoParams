import glob
from rdkit import Chem
from rdkit.Chem import Draw
from .defaults import *
from .utilities import *

def PDBtoChemDraw(pdb,output):
    mol = Chem.MolFromPDBFile(pdb)
    smiles = Chem.MolToSmiles(mol)
    smi_file = pdb.replace(".pdb",".smi")
    with open(smi_file,"w") as f:
        f.write(smiles)
        f.write("    ")
        f.write(pdb)
    Draw.MolToFile(Chem.MolFromSmiles(smiles), output, size=(1268,720))
    return smiles
