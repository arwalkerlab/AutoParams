class Mol2File():
    def __init__(self,filename,charge):
        self.filename = filename
        self.atoms = []
        self.bonds = []
        self.headtail = ""
        self.resconnect = ""
        self.charge = charge
        self.resname = ""
    def AddAtom(self,atomname,x,y,z,atomtype,resname,respcharge):
        self.atoms.append(f"{len(self.atoms)+1:>7} {atomname:<4} {float(x):>14.4f} {float(y):>10.4f} {float(z):>10.4f} {atomtype:<2} {1:>10} {resname:<4} {float(respcharge):>14.6f}\n")
    def AddBond(self,atom1,atom2,order):
        self.bonds.append(f"{len(self.bonds)+1:>6} {atom1:>5} {atom2:>5} {order:>1}\n")
    def WriteFile(self):
        with open(self.filename,"w") as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{self.resname:>3}\n")
            f.write(f"{len(self.atoms):>5} {len(self.bonds):>5}     1     0     0\n")
            f.write(f"{len(self.atoms):>5} {len(self.bonds):>5}     1     0     0\n")
            f.write(f"{self.charge}\n")
            f.write("\n\n")
            f.write("@<TRIPOS>ATOM\n")
            for atom in self.atoms:
                f.write(str(atom))
            f.write("@<TRIPOS>BOND\n")
            for bond in self.bonds:
                f.write(str(bond))
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write(f"     1   {self.resname:>3}     1 TEMP\n")
            if self.headtail != "":
                f.write("@<TRIPOS>HEADTAIL\n")
                f.write(f"{self.headtail}\n")
            if self.resconnect != "":
                f.write("@<TRIPOS>RESIDUECONNECT\n")
                f.write(f"{self.resconnect}\n")
