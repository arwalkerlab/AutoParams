from .ParameterUtilities import *
from .AtomTyper import *
from .classes import *
import parmed


def GetPDBAtoms(pdbfile):
    tmp = parmed.load_file(pdbfile)
    return [atom.name for atom in tmp.atoms]

def GetAtomTypes(pdbfile:str,template:str):
    base_dict = nucleobase_name_type_dicts[template]
    atomtypes = []
    for line in open(pdbfile).readlines():
        if any(["ATOM" in line, "HETATM" in line]):
            atomname = line.split()[2]
            if atomname in base_dict.keys():
                atomtypes.append(base_dict[atomname])
            else:
                atomtypes.append(GenerateAtomTypeFromPDB(pdbfile,atomname))
    return atomtypes
        
def Generate_Mol2_Bonds(pdbfile):
    pdb = pmd.load_file(pdbfile)
    mol2_bonds = []
    bond_orders = {}
    for key in [atom.name for atom in pdb.atoms]:
        bond_orders[key] = 0
    max_bond_count = {"C":4,"N":3,"O":2,"S":5,"P":6}
    found_bonds = []
    for main_atom in pdb.atoms:
        for bond in main_atom.bonds:
            a1 = bond.atom1
            a2 = bond.atom2
            if sorted([a1.idx+1,a2.idx+1]) in found_bonds:
                continue
            found_bonds.append(sorted([a1.idx+1,a2.idx+1]))
            bond_order = 0
            if any([a1.element_name == "H",a2.element_name == "H"]):
                # All bonds with hydrogens are order 1
                bond_order = 1
            elif any([all([bond.atom2.element_name == "O",len(bond.atom2.bonds)==2]),all([bond.atom1.element_name == "O",len(bond.atom1.bonds)==2])]):
                # All oxygens with two bonds MUST have bond order of 1 for each bond
                bond_order = 1
            elif any([all([a1.element_name == "O", len(a1.bonds)==1,a2.element_name=="C", len(a2.bonds)==3]),all([a2.element_name == "O", len(a2.bonds)==1,a1.element_name=="C", len(a1.bonds)==3])]):
                # All bonds with carbonyl oxygens have a bond order of 2.
                # requires:  sp2 carbon in bond, oxygen with only one bonded atom (how to deal with deprotonated?)
                bond_order = 2
            elif any([all([a2.element_name == "C", len(a2.bonds)==4]),all([a1.element_name == "C", len(a1.bonds)==4])]):
                # All bonds to sp3 carbons must have bond order of 1
                bond_order = 1
            elif any([all([a1.element_name == "C",len(a1.bonds)==3]),all([a2.element_name == "C",len(a2.bonds)==3])]):
                # sp2 carbons have 3 bonded atoms, one of the bonds is bond-order 2
                avail_bonds_1 = max_bond_count[a1.element_name] - bond_orders[a1.name]
                avail_bonds_2 = max_bond_count[a2.element_name] - bond_orders[a2.name]
                if all([avail_bonds_1 > 1,avail_bonds_2 > 1]):
                    bond_order = 2
                else:
                    bond_order = 1
            elif any([a1.element_name == "P",a2.element_name == "P"]):
                avail_bonds_1 = max_bond_count[a1.element_name] - bond_orders[a1.name]
                avail_bonds_2 = max_bond_count[a2.element_name] - bond_orders[a2.name]
                if all([avail_bonds_1 > 1,avail_bonds_2 > 1,any([a1.element_name=="O",a2.element_name=="O"])]): ## Deals with phosphate groups
                    bond_order = 2
                else:
                    bond_order = 1
            else:
                report_failure(f"Unable to assign bond order. {bond}")

            if bond_order > 0:
                atom_list = sorted([a1.idx+1,a2.idx+1])
                mol2bond = [atom_list[0],atom_list[1],bond_order]
                mol2_bonds.append(mol2bond)
                bond_orders[bond.atom1.name] += bond_order
                bond_orders[bond.atom2.name] += bond_order
            else:
                report_failure(f"Unable to assign bond order for {bond}.  Check structure of {pdbfile}")
    # Resort Mol2Bonds to be ordered numerically by first and second atoms
    recheck = True
    while recheck:
        recheck = False
        for i,bond in enumerate(mol2_bonds[:-1]):
            if bond[0] == mol2_bonds[i+1][0]:
                if bond[1] > mol2_bonds[i+1][1]:
                    temp_2 = bond[1]
                    temp_ord = bond[2]

                    bond[1] = mol2_bonds[i+1][1]
                    mol2_bonds[i+1][1] = temp_2
                    bond[2] = mol2_bonds[i+1][2]
                    mol2_bonds[i+1][2] = temp_ord
                    recheck = True
    return mol2_bonds

def PDB_to_Mol2(pdbfile:str,respcharges:list,template:str,charge:int):
    # Initialize Mol2File
    mol2filename = pdbfile.replace(".pdb",".mol2")
    Mol2obj = Mol2File(mol2filename,charge)
    # Generate Atoms,X,Y,Z,Types,Resname,RESPCharges.
    resname = GetResName(pdbfile)
    atomlist = GetPDBAtoms(pdbfile,resname)
    atomtypes = GetAtomTypes(pdbfile,template)
    charges = respcharges
    n_atoms = len(atomlist)
    if not all([ len(atomtypes) == n_atoms,len(charges) == n_atoms]):
        report_failure("Inconsistent number of atomnames,charges, and/or atomtypes.")
    for i in range(n_atoms):
        atomname = atomlist[i][1]
        x = atomlist[i][2]
        y = atomlist[i][3]
        z = atomlist[i][4]
        atomtype = atomtypes[i]
        respcharge = charges[i]
        Mol2obj.AddAtom(atomname,x,y,z,atomtype,resname,respcharge)
    for bond in Generate_Mol2_Bonds(pdbfile):
        Mol2obj.AddBond(bond[0],bond[1],bond[2])
    Mol2obj.WriteFile()
