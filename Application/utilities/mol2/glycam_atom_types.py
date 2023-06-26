import parmed

def RewriteMol2Glycam(mol2file,pdbfile):
    tmpmol = Molecule(pdbfile=pdbfile)
    GlycamAtomTypes(tmpmol)
    mol2lines = open(mol2file).readlines()
    n_atoms = int(mol2lines[2].strip().split()[0])
    first_chunk = []
    for i,line in enumerate(mol2lines):
        first_chunk.append(line)
        if "@<TRIPOS>ATOM" in line:
            break
    last_chunk = []
    atom_chunk = []
    for line in mol2lines[i+1:i+n_atoms+1]:
        atom_chunk.append(line)
    for line in mol2lines[i+n_atoms+1:]:
        last_chunk.append(line)
    for atom in tmpmol.atoms:
        if atom.atomtype != "XX":
            tmp_line = atom_chunk[atom.index-1][:50] + atom.atomtype + atom_chunk[atom.index-1][52:]
            atom_chunk[atom.index-1]=tmp_line
    with open(mol2file,"w") as f:
        f.write("".join(first_chunk))
        f.write("".join(atom_chunk))
        f.write("".join(last_chunk))




class Molecule():
    def __init__(self,pdbfile):
        self.pdb = parmed.load_file(pdbfile)
        self.atoms = []
        self.bonds = []
        for atom in self.pdb.atoms:
            tmpatom = Atom()
            tmpatom.FromParmedAtom(atom)
            self.atoms.append(tmpatom)
        for bond in self.pdb.bonds:
            a1 = bond.atom1.idx
            a2 = bond.atom2.idx
            self.bonds.append([self.atoms[a1],self.atoms[a2]])

class Atom():
    def __init__(self):
        self.element="X"
        self.atomname="XX"
        self.atomtype="XX"
        self.x = 0.00
        self.y = 0.00
        self.z = 0.00
        self.n_bonds = 0
        self.partial_charge = 0.0
        self.index = 0
    def FromParmedAtom(self,atom):
        self.__init__()
        self.element = atom.element_name
        self.atomname = atom.name
        self.x = float(atom.xx)
        self.y = float(atom.xy)
        self.z = float(atom.xz)
        self.n_bonds = len(atom.bonds)
        self.index = atom.idx+1
        self.resname = atom.residue.name
    def SetAtomType(self,atomtype):
        self.atomtype = atomtype
    def SetPartialCharge(self,charge):
        self.partial_charge = charge
    def Mol2String(self):
        return f"{self.index:>7} {self.atomname:<4} {self.x:>14.4f} {self.y:>10.4f} {self.z:>10.4f} {self.atomtype:<2} {1:>10} {self.resname:<4} {self.partial_charge:>14.6f}\n"
    def PDBString(self):
        pass
    

##### START OF GLYCAM ATOMTYPING #####
def Glycam_Hydrogens(mol):
    for bond in mol.bonds:
        if bond[0].element == "H":
            H_atom = bond[0]
            other_element = bond[1]
        elif bond[1].element == "H":
            H_atom = bond[1]
            other_element = bond[0]
        else:
            continue
        if H_atom.atomtype != "XX":
            continue
        if other_element.element == "N":
            # "H bonded to nitrogen"
            H_atom.atomtype="H "
        elif other_element.element == "O":
            # "H in hydroxyl group"
            H_atom.atomtype="Ho"
        elif other_element.element == "C":
            if other_element.n_bonds == 4:
                #sp3 carbon
                H_atom.atomtype="Hc" #H bonded to aliphatic C, no EWD
            elif other_element.n_bonds == 3:
                H_atom.atomtype="Ha" #H bonded to C in alkenes
                
def Glycam_Sulfurs(mol):
    for atom in mol.atoms:
        if any([atom.element != "S",atom.atomtype != "XX"]):
            continue
        bonded_elements=[]
        for bond in mol.bonds:
            if bond[0].index == atom.index:
                bonded_elements.append(bond[1].element)
            elif bond[1].index == atom.index:
                bonded_elements.append(bond[0].element)
        def count_in_array(array,element):
            i=0
            for val in array:
                if val == element:
                    i+=1
            return i
        if count_in_array(bonded_elements,"O") == 3:
            # Sulfur in Sulfates
            atom.atomtype="S "
        elif count_in_array(bonded_elements,"C") == 2:
            # sulfure in carbohydrate linkages (thioether)
            atom.atomtype="Sm"
            
def Glycam_Oxygens(mol):
    for atom in mol.atoms:
        if any([atom.element != "O",atom.atomtype != "XX"]):
            continue
        bonded_elements=[]
        hold_bonds = []
        for bond in mol.bonds:
            if bond[0].index == atom.index:
                bonded_elements.append(bond[1].element)
                hold_bonds.append(bond)
            elif bond[1].index == atom.index:
                bonded_elements.append(bond[0].element)
                hold_bonds.append(bond)
        def count_in_array(array,element):
            i=0
            for val in array:
                if val == element:
                    i+=1
            return i
        if count_in_array(bonded_elements,"C") == 2:
            # oxygen in ether linkage
            atom.atomtype = "Os"
        elif count_in_array(bonded_elements,"H") == 1:
            #hydroxyl group oxygen
            atom.atomtype = "Oh"
        elif atom.n_bonds == 1:
            if "C" in bonded_elements:
                # get C index to get all other atoms:
                # carbonyl/carboxylate
                for bond in hold_bonds:
                    if bond[0].index == atom.index:
                        c_index = bond[1].index
                    else:
                        c_index = bond[0].index
                for bond in mol.bonds:
                    if any([bond[0].index == c_index,bond[1].index == c_index]):
                        if "O" == bond[0].element:
                            index = bond[0].index
                        elif "O" == bond[1].element:
                            index = bond[1].index
                        else:
                            continue
                        if index < atom.index:
                            atom.atomtype = "O2"
                        else:
                            atom.atomtype = "O "
                        break
            elif "P" in bonded_elements:
                # phosphate group oxygen
                atom.atomtype = "O2"
            elif "S" in bonded_elements:
                # phosphate group oxygen
                atom.atomtype = "O2"
            
def Glycam_Nitrogens(mol):
    for atom in mol.atoms:
        if any([atom.element != "N",atom.atomtype != "XX"]):
            continue
        if atom.n_bonds==1:
            # nitrile, not handled by glycam
            continue
        if atom.n_bonds==2:
            # sp2, such as found in aromatic heterocycles
            atom.atomtype="Ng"
            continue
        if atom.n_bonds==4:
            # four-connections = sp3 charged Nitrogen
            atom.atomtype = "N3"
            continue
        bonded_elements=[]
        for bond in mol.bonds:
            if bond[0].index == atom.index:
                bonded_elements.append(bond[1].element)
            elif bond[1].index == atom.index:
                bonded_elements.append(bond[0].element)
        def count_in_array(array,element):
            i=0
            for val in array:
                if val == element:
                    i+=1
            return i
        if any([count_in_array(bonded_elements,"O") == 1,count_in_array(bonded_elements,"C") == 2]):
            # Nitrogen bonded to oxygen in sulfates
            # OR
            # Nitrogen as secondary amine/amide
            atom.atomtype="Ng"
            continue
        else:
            # sp3 nitrogen
            atom.atomtype="NT"

def Glycam_Carbons(mol):
    for atom in mol.atoms:
        if any([atom.element != "C",atom.atomtype != "XX"]):
            continue
        if atom.n_bonds==1:
            # weird and broken carbon, not handled by glycam
            continue
        if atom.n_bonds==2:
            # sp, such as in nitriles OR between two double-bonds.  Either way, not handled by glycam
            continue
        bonded_elements=[]
        for bond in mol.bonds:
            if bond[0].index == atom.index:
                bonded_elements.append(bond[1].element)
            elif bond[1].index == atom.index:
                bonded_elements.append(bond[0].element)
        def count_in_array(array,element):
            i=0
            for val in array:
                if val == element:
                    i+=1
            return i
        if atom.n_bonds==3:
            #sp2 carbons, check for carbonyl, then regular sp2 carbon
            if count_in_array(bonded_elements,"O")>0:
                atom.atomtype="C "
                continue
            atom.atomtype="Ck"
            continue
        if atom.n_bonds != 4:
            # error handling for strange connections, break out from this carbon and leave it for later.
            continue
        # four-connections = sp3 carbon
        if sorted(list(set(["H","C","N"]))) == sorted(list(set(bonded_elements))):
            #protein alpha-carbon
            atom.atomtype = "CX"
            continue
        atom.atomtype = "Cg"

def Glycam_Phosphorus(mol):
    for atom in mol.atoms:
        if all([atom.element == "P",atom.atomtype != "XX",atom.n_bonds==4]):
            atom.atomtype = "P "

def GlycamAtomTypes(mol):
    Glycam_Hydrogens(mol)
    Glycam_Sulfurs(mol)
    Glycam_Oxygens(mol)
    Glycam_Nitrogens(mol)
    Glycam_Carbons(mol)
    Glycam_Phosphorus(mol)
##### END OF GLYCAM ATOMTYPING #####

