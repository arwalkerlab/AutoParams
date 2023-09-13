import os
this_folder = os.path.dirname(os.path.abspath(__file__))

MASSES_LIB_FILE = this_folder+"/masses.txt"
BONDS_LIB_FILE = this_folder+"/bonds.txt"
ANGLES_LIB_FILE = this_folder+"/angles.txt"
DIHEDRALS_LIB_FILE = this_folder+"/dihedrals.txt"
TORSIONS_LIB_FILE = this_folder+"/torsions.txt"
NONBONDED_LIB_FILE = this_folder+"/nonbonded.txt"

def CleanBondsLibraryFile():
    new_file_dict = []
    for line in open(BONDS_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:5]
        [a1,a2] = key.split("-")
        a1=a1.strip()
        a2=a2.strip()
        [a1,a2] = sorted([a1,a2])
        newline = [a1,a2,line[5:].replace("\n","")]
        if newline not in new_file_dict:
            new_file_dict.append(newline)
    new_file_dict = sorted(new_file_dict,key=lambda x: x[0])
    lines = []
    for newline in new_file_dict:
        lines.append(f"{newline[0]:<2}-{newline[1]:<2}{newline[2]}")
    with open(BONDS_LIB_FILE,"w") as f:
        f.write("\n".join(x for x in lines))
        f.write("\n")

def CleanAnglesLibraryFile():
    new_file_dict = []
    for line in open(ANGLES_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:8]
        [a1,a2,a3] = key.split("-")
        a1=a1.strip()
        a2=a2.strip()
        a3=a3.strip()
        [a1,a3] = sorted([a1,a3])
        newline = [a1,a2,a3,line[8:].replace("\n","")]
        if newline not in new_file_dict:
            new_file_dict.append(newline)
    new_file_dict = sorted(new_file_dict,key=lambda x: x[1])
    lines = []
    for newline in new_file_dict:
        lines.append(f"{newline[0]:<2}-{newline[1]:<2}-{newline[2]:<2}{newline[3]}")
    with open(ANGLES_LIB_FILE,"w") as f:
        f.write("\n".join(x for x in lines))
        f.write("\n")

def CleanDihedralsLibraryFile():
    new_file_dict = []
    for line in open(DIHEDRALS_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:11]
        [a1,a2,a3,a4] = key.split("-")
        a1=a1.strip()
        a2=a2.strip()
        a3=a3.strip()
        if [a3,a2] == sorted([a2,a3]):
            [a1,a2,a3,a4] = [a4,a3,a2,a1]
        newline = [a1,a2,a3,a4,line[11:].replace("\n","")]
        if newline not in new_file_dict:
            new_file_dict.append(newline)
    new_file_dict = sorted(new_file_dict,key=lambda x: [x[1],x[2]])
    lines = []
    for newline in new_file_dict:
        lines.append(f"{newline[0]:<2}-{newline[1]:<2}-{newline[2]:<2}-{newline[3]:<2}{newline[4]}")
    with open(DIHEDRALS_LIB_FILE,"w") as f:
        f.write("\n".join(x for x in lines))
        f.write("\n")

def CleanTorsionsLibraryFile():
    new_file_dict = []
    for line in open(TORSIONS_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:11]
        [a1,a2,a3,a4] = key.split("-")
        a1=a1.strip()
        a2=a2.strip()
        a3=a3.strip()
        if [a3,a2] == sorted([a2,a3]):
            [a1,a2,a3,a4] = [a4,a3,a2,a1]
        newline = [a1,a2,a3,a4,line[11:].replace("\n","")]
        if newline not in new_file_dict:
            new_file_dict.append(newline)
    new_file_dict = sorted(new_file_dict,key=lambda x: [x[1],x[2]])
    lines = []
    for newline in new_file_dict:
        lines.append(f"{newline[0]:<2}-{newline[1]:<2}-{newline[2]:<2}-{newline[3]:<2}{newline[4]}")
    with open(TORSIONS_LIB_FILE,"w") as f:
        f.write("\n".join(x for x in lines))
        f.write("\n")

def GenerateParameterDictionaries():
    # Initialize dictionaries of bonds, angles, dihedrals, and torsions that may be missing.
    known_masses_dict    = {}
    known_bonds_dict     = {}
    known_angles_dict    = {}
    known_dihedrals_dict = {}
    known_torsions_dict  = {}
    known_nonbonded_dict = {}
    # Process known masses from local masses.txt file
    for line in open(MASSES_LIB_FILE).readlines():
        if line.strip()=="":
            continue
        key = line[:2]
        known_masses_dict[key] = line

    # Process known bonds from local bonds.txt file
    for line in open(BONDS_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:5]
        known_bonds_dict[key] = line

    # Process known angles from local angles.txt file.
    for line in open(ANGLES_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:8]
        known_angles_dict[key] = line

    # Process known dihedrals from local dihedrals.txt file.
    for line in open(DIHEDRALS_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:11]
        # because there can be multiple lines for a given dihedral, 
        # we will include them all as a single entry in the dictionary.
        if key not in known_dihedrals_dict.keys():
            known_dihedrals_dict[key] = line
        else:
            known_dihedrals_dict[key] += "\n"
            known_dihedrals_dict[key] += line

    # Process known improper torsions from local torsions.txt file
    for line in open(TORSIONS_LIB_FILE).readlines():
        if line.strip() == "":
            continue
        key = line[:11]
        known_torsions_dict[key] = line

    # Process known nonbonded from local nonbonded.txt file
    for line in open(NONBONDED_LIB_FILE).readlines():
        if line.strip()=="":
            continue
        key = line[:2]
        known_nonbonded_dict[key] = line

    return {"MASSES":known_masses_dict,
            "BONDS":known_bonds_dict,
            "ANGLES":known_angles_dict,
            "DIHEDRALS":known_dihedrals_dict,
            "TORSIONS":known_torsions_dict,
            "NONBONDED":known_nonbonded_dict}