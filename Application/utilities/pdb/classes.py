import parmed

class Atom():
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.resp_charge = None
        self.atom_name = None
        self.atom_type = None
        pass
    def set_position(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    def set_resp_charge(self,respcharge):
        self.resp_charge = respcharge
    def set_atom_name(self,atomname):
        self.atom_name = atomname
    def set_atom_type(self,atomtype):
        self.atom_type = atomtype

class Bond():
    def __init__(self):
        self.atom1 = ""
        self.atom2 = ""
        self.order = 0
    
class PDBGeometry():
    def __init__(self,source_file):
        self.ParseBonds(source_file)
        self.ParseAngles()
        self.ParseDihedrals()
        
    def ParseBonds(self,pdbfile):
        tmp = parmed.load_file(pdbfile)
        for bond in tmp.bonds:
            self.bonds.append(sorted([bond.atom1.idx,bond.atom2.idx]))
        self.bonds = sorted(self.bonds,key=lambda x: x[0])
    
    def ParseAngles(self):
        self.angles = []
        for i,[a1,a2] in enumerate(self.bonds):
            for [a3,a4] in self.bonds[i+1:]:
                if a4 == a1:
                    if [a3,a4,a2] not in self.angles:
                        self.angles.append([a3,a4,a2])
                elif a3 == a2:
                    if [a1,a2,a4] not in self.angles:
                        self.angles.append([a1,a2,a4])
                elif a3 == a1:
                    if [a2,a1,a3] not in self.angles:
                        self.angles.append([a2,a1,a3])
                elif a4 == a2:
                    if [a1,a2,a3] not in self.angles:
                        self.angles.append([a1,a2,a3])
        self.angles = sorted(self.angles,key=lambda x: x[1])

    def ParseDihedrals(self,):
        self.dihedrals = []
        for i,[a1,a2,a3] in enumerate(self.angles):
            for [b1,b2,b3] in self.angles[i+1:]:
                if all([a2 == b1, a3==b2, [a1,a2,a3,b3] not in self.dihedrals]):
                    self.dihedrals.append([a1,a2,a3,b3])
                elif all([a2 == b3, a1==b2, [b1,a1,a2,a3] not in self.dihedrals]):
                    self.dihedrals.append([b1,a1,a2,a3])