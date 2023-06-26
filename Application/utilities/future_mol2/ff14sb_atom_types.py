##### START OF FF14SB ATOMTYPING #####
def ff14SB_Hydrogens(mol):
    for atom in mol.atoms:
        if any([atom.element!= "H",atom.atomtype != "XX"]):
            continue
        bond = atom.bonds[0]
        if bond[0].element == "H":
            H_atom = bond[0]
            other_element = bond[1]
        elif bond[1].element == "H":
            H_atom = bond[1]
            other_element = bond[0]
        if other_element.element == "N":
            # "H bonded to nitrogen"
            H_atom.atomtype = "H "
            continue
        if other_element.element == "O":
            # "H in hydroxyl group"
            H_atom.atomtype = "HO"
            continue
        if other_element.element == "S":
            # H bonded to sulfur
            H_atom.atomtype = "HS"
            continue
        if other_element.element == "C":
            if other_element.n_bonds == 4:
                #sp3 carbon
                H_atom.atomtype="HC" #H bonded to aliphatic C, no EWD
                continue
            if other_element.n_bonds == 3:
                H_atom.atomtype="HA" #H bonded to C in alkenes
                continue
            if other_element.n_bonds == 2:
                H_atom.atomtype="HZ" #H bonded to sp C
                continue

def ff14SB_Sulfurs(mol):
    for atom in mol.atoms:
        if any([atom.element != "S",atom.atomtype != "XX"]):
            continue
        bonded_elements=[]
        for bond in mol.bonds:
            if bond[0].index == atom.index:
                bonded_elements.append(bond[1].element)
            elif bond[1].index == atom.index:
                bonded_elements.append(bond[0].element)
        if "S" in bonded_elements:
            # Sulfur in disulfide linkage
            atom.atomtype="S "
            continue
        if "H" in bonded_elements:
            atom.atomtype="SH"
            continue
            
def ff14SB_Oxygens(mol):
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
            # oxygen in ether/ester linkage
            atom.atomtype = "OS"
        elif count_in_array(bonded_elements,"H") == 1:
            #hydroxyl group oxygen
            atom.atomtype = "OH"
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
                # sulfate group oxygen
                atom.atomtype = "O2"
            
def ff14SB_Nitrogens(mol):
    for atom in mol.atoms:
        if any([atom.element != "N",atom.atomtype != "XX"]):
            continue
        if atom.n_bonds==1:
            # nitrile
            atom.atomtype="NY"
            continue
        if atom.n_bonds==2:
            # sp2, such as found in aromatic heterocycles
            atom.atomtype="NB"
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
        if "H" in bonded_elements:
            atom.atomtype="NA"
            continue
        if any([count_in_array(bonded_elements,"O") == 1,count_in_array(bonded_elements,"C") == 2]):
            # Nitrogen bonded to oxygen in sulfates
            # OR
            # Nitrogen as secondary amine/amide
            atom.atomtype="N "
            continue
        # sp3 nitrogen
        atom.atomtype="NT"

def ff14SB_Carbons(mol):
    for atom in mol.atoms:
        if any([atom.element != "C",atom.atomtype != "XX"]):
            continue
        if atom.n_bonds==1:
            # weird and broken carbon
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
        if atom.n_bonds==2:
            # sp, such as in nitriles OR between two double-bonds.
            if "N" in bonded_elements:
                atom.atomtype="CY"
                continue
            atom.atomtype="CZ"
            continue
        if atom.n_bonds==3:
            #sp2 carbons, check for carbonyl, then regular sp2 carbon
            if count_in_array(bonded_elements,"O")>0:
                atom.atomtype="C "
                continue
            atom.atomtype="CA"
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

def ff14SB_Phosphorus(mol):
    for atom in mol.atoms:
        if all([atom.element == "P",atom.atomtype != "XX",atom.n_bonds==4]):
            atom.atomtype = "P "

def FF14SBAtomTypes(mol):
    ff14SB_Hydrogens(mol)
    ff14SB_Sulfurs(mol)
    ff14SB_Oxygens(mol)
    ff14SB_Nitrogens(mol)
    ff14SB_Carbons(mol)
    ff14SB_Phosphorus(mol)
##### END OF FF14SB ATOMTYPING #####

