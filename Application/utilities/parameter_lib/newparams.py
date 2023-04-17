BASIC_ATOM_TYPES   = {"C_sp":["CZ","CY"],
                      "C_sp2":["C4","C5","C ","CA","CB","CC","CD","CK","CP","CM","CS","CN","CQ","CR","CV","CW","C*","CI"],
                      "C_sp3":["CX","CT"],
                      "N_sp":["NY"],
                      "N_sp2":["N ","NA","NB","NC","N2","N*"],
                      "N_sp3":["N3","NT"],
                      "O_sb":["OW","OH","OS"],
                      "O_db":["OP","O2","O "],
                      "H":["HW","H4","H5","HA","HC","H1","H2","H3","HP","HZ","H ","HO","HS"],
                      "F":["F "],
                      "Cl":["Cl"],
                      "Br":["Br"],
                      "I":["I "],
                      "S":["S ","SH"],
                      "P":["P "],
                      "EP":["EP"]
                     }

BOND_LENGTH_DICT = { 'H-O_sb':0.970, 'C_sp2-H':1.120, 'C_sp3-H':1.120, 
                     'C_sp-H':1.120, 'H-N_sp2':1.010, 'H-N_sp3':1.010, 
                     'H-S':1.341, 'H-H':0.737, 'C_sp2-C_sp2':1.340, 
                     'C_sp2-C_sp3':1.500, 'C_sp2-N_sp2':1.279, 
                     'C_sp2-O_db':1.230, 'C_sp2-O_sb':1.360, 
                     'C_sp3-C_sp3':1.540,'C_sp3-N_sp2':1.465, 
                     'C_sp3-O_sb':1.470, 'C_sp3-N_sp3':1.469, 
                     'C_sp3-S':1.810, 'C_sp-C_sp3':1.460, 
                     'C_sp-N_sp':1.154, 'C_sp-C_sp':1.206, 
                     'O_db-P':1.530, 'O_sb-P':1.530,'N_sp2-P':1.840,
                     'S-S':2.038, 'C_sp3-F':1.380,'C_sp2-F':1.380, 
                     'C_sp3-Cl':1.766,'C_sp2-Cl':1.766, 'C_sp3-I':2.166, 
                     'C_sp2-I':2.166, 'Br-C_sp3':1.944, 'Br-C_sp2':1.944, 
                     'EP-O_db':0.200, 'EP-O_sb':0.200, 'EP-N_sp3':0.200,
                     'EP-N_sp2':0.200, 'EP-S':0.700
                    }

BOND_FORCE_CONST_DICT = {'H-O_sb':553.0, 'C_sp2-H':367, 'C_sp3-H':340, 
                         'C_sp-H':400, 'H-N_sp2':434, 'H-N_sp3':434, 
                         'H-S':274, 'H-H':553, 'C_sp2-C_sp2':410, 
                         'C_sp2-C_sp3':317,'C_sp2-N_sp2':454, 
                         'C_sp2-O_db':613,'C_sp2-O_sb':450, 
                         'C_sp3-C_sp3':310,'C_sp3-N_sp2':337, 
                         'C_sp3-O_sb':320,'C_sp3-N_sp3':367, 
                         'C_sp3-S':232,'C_sp-C_sp3':400,'C_sp-N_sp':600, 
                         'C_sp-C_sp':600,'O_db-P':525,'O_sb-P':230,
                         'N_sp2-P':250,'S-S':166,'C_sp3-F':367,
                         'C_sp2-F':386,'C_sp3-Cl':232,'C_sp2-Cl':193,
                         'C_sp3-I':148,'C_sp2-I':171,'Br-C_sp3':159,
                         'Br-C_sp2':172,'EP-O_db':600,'EP-O_sb':600,
                         'EP-N_sp3':600,'EP-N_sp2':600,'EP-S':600
                    }

ANGLE_MEASURE_DICT = {
"H-O_sb-H":   104.520,"H-H-O_sb":   127.740,"C_sp2-C_sp2-O_db":   124.940,"C_sp2-C_sp2-O_sb":   121.667,
"C_sp2-C_sp2-C_sp2":   118.930,"C_sp2-N_sp2-P":   125.100,"N_sp2-P-O_db":   102.380,
"C_sp2-C_sp2-N_sp2":   119.121,"C_sp3-C_sp2-O_db":   118.700,"C_sp3-C_sp2-N_sp2":   118.300,
"C_sp3-C_sp2-C_sp3":   117.000,"C_sp3-C_sp2-O_sb":   110.800,"N_sp2-C_sp2-N_sp2":   119.200,
"N_sp2-C_sp2-O_db":   121.725,"O_db-C_sp2-O_db":   126.000,"O_db-C_sp2-O_sb":   122.500,
"C_sp2-C_sp2-H":   120.177,"C_sp3-C_sp2-H":   116.125,"H-C_sp2-O_db":   119.500,"H-C_sp2-O_sb":   111.286,
"H-C_sp2-N_sp2":   120.737,"C_sp2-C_sp2-C_sp3":   120.843,"C_sp2-C_sp2-F":   121.000,
"C_sp2-C_sp2-Cl":   118.800,"Br-C_sp2-C_sp2":   118.800,"C_sp2-C_sp2-I":   118.800,"H-C_sp2-H":   117.125,
"H-C_sp3-H":   109.500,"H-C_sp3-N_sp2":   109.500,"H-C_sp3-O_sb":   109.500,"C_sp2-C_sp3-H":   109.500,
"C_sp-C_sp3-H":   110.000,"H-C_sp3-S":   109.500,"H-C_sp3-N_sp3":   109.500,"C_sp2-C_sp3-N_sp2":   110.100,
"C_sp2-C_sp3-N_sp3":   111.200,"C_sp2-C_sp3-C_sp3":   112.517,"C_sp2-C_sp3-O_sb":   109.500,
"C_sp3-C_sp3-C_sp3":   109.500,"C_sp3-C_sp3-H":   109.500,"C_sp3-C_sp3-N_sp2":   110.025,
"C_sp3-C_sp3-O_sb":   109.500,"C_sp3-C_sp3-S":   111.650,"C_sp3-C_sp3-N_sp3":   111.200,
"C_sp-C_sp3-C_sp3":   110.000,"O_sb-C_sp3-O_sb":   101.000,"C_sp-C_sp3-O_sb":   110.000,
"N_sp2-C_sp3-O_sb":   109.500,"F-C_sp3-F":   109.100,"F-C_sp3-H":   109.500,"C_sp3-C_sp3-F":   109.000,
"C_sp3-C_sp3-Cl":   108.500,"Cl-C_sp3-H":   108.500,"Br-C_sp3-C_sp3":   108.000,"Br-C_sp3-H":   106.500,
"C_sp3-C_sp3-I":   106.000,"C_sp3-C_sp-N_sp":   180.000,"C_sp-C_sp-C_sp3":   180.000,
"C_sp-C_sp-H":   180.000,"C_sp2-N_sp2-C_sp3":   123.673,"C_sp2-N_sp2-H":   121.636,
"C_sp3-N_sp2-H":   118.160,"C_sp3-N_sp2-C_sp3":   118.000,"H-N_sp2-H":   120.000,
"C_sp2-N_sp2-C_sp2":   114.595,"C_sp3-N_sp3-H":   109.500,"C_sp3-N_sp3-C_sp3":   109.500,
"H-N_sp3-H":   109.500,"C_sp2-O_sb-H":   111.500,"C_sp3-O_sb-H":   108.500,"H-O_sb-P":   108.500,
"C_sp2-O_sb-C_sp3":   117.000,"C_sp3-O_sb-C_sp3":   109.500,"C_sp3-O_sb-P":   120.500,
"C_sp2-O_sb-P":   120.500,"P-O_sb-P":   120.500,"O_db-P-O_sb":   108.230,"O_db-P-O_db":   119.900,
"O_sb-P-O_sb":   102.600,"C_sp3-S-C_sp3":    98.900,"C_sp3-S-S":   103.700,"C_sp3-S-H":    96.000,
"H-S-H":    92.070,"C_sp2-N_sp2-EP":   123.600,"C_sp3-N_sp3-EP":   109.500,"EP-N_sp3-H":   109.500,
"C_sp2-O_db-EP":   120.000,"EP-O_db-EP":   120.000,"C_sp2-O_sb-EP":   112.125,"C_sp3-O_sb-EP":   109.500,
"EP-O_sb-H":   109.500,"EP-O_sb-EP":   109.500,"C_sp3-S-EP":    90.000,"EP-O_sb-P":   109.500,
"EP-S-EP":   180.000,"EP-S-H":    90.000,"EP-S-S":    96.700,"O_sb":   112.668,"H":   127.740,
"C_sp2":   119.322,"N_sp2":   119.873,"P":   109.009,"C_sp3":   109.920,"C_sp":   180.000,
"N_sp3":   109.500,"S":   111.737,"O_db":   120.000}

ANGLE_FORCE_CONST_DICT = {
"H-O_sb-H":   100.000,"H-H-O_sb":     0.000,"C_sp2-C_sp2-O_db":    80.000,"C_sp2-C_sp2-O_sb":    75.000,
"C_sp2-C_sp2-C_sp2":    63.000,"C_sp2-N_sp2-P":    76.700,"N_sp2-P-O_db":    42.900,"C_sp2-C_sp2-N_sp2":    70.000,
"C_sp3-C_sp2-O_db":    75.000,"C_sp3-C_sp2-N_sp2":    70.000,"C_sp3-C_sp2-C_sp3":    63.000,
"C_sp3-C_sp2-O_sb":    68.000,"N_sp2-C_sp2-N_sp2":    70.000,"N_sp2-C_sp2-O_db":    80.000,
"O_db-C_sp2-O_db":    80.000,"O_db-C_sp2-O_sb":    80.000,"C_sp2-C_sp2-H":    50.000,
"C_sp3-C_sp2-H":    50.000,"H-C_sp2-O_db":    50.000,"H-C_sp2-O_sb":    50.000,"H-C_sp2-N_sp2":    50.000,
"C_sp2-C_sp2-C_sp3":    70.000,"C_sp2-C_sp2-F":    70.000,"C_sp2-C_sp2-Cl":    70.000,"Br-C_sp2-C_sp2":    70.000,
"C_sp2-C_sp2-I":    70.000,"H-C_sp2-H":    35.000,"H-C_sp3-H":    35.000,"H-C_sp3-N_sp2":    50.000,
"H-C_sp3-O_sb":    50.000,"C_sp2-C_sp3-H":    50.000,"C_sp-C_sp3-H":    50.000,"H-C_sp3-S":    50.000,
"H-C_sp3-N_sp3":    50.000,"C_sp2-C_sp3-N_sp2":    63.000,"C_sp2-C_sp3-N_sp3":    80.000,
"C_sp2-C_sp3-C_sp3":    61.083,"C_sp2-C_sp3-O_sb":    52.500,"C_sp3-C_sp3-C_sp3":    40.000,
"C_sp3-C_sp3-H":    50.000,"C_sp3-C_sp3-N_sp2":    72.500,"C_sp3-C_sp3-O_sb":    50.000,
"C_sp3-C_sp3-S":    50.000,"C_sp3-C_sp3-N_sp3":    80.000,"C_sp-C_sp3-C_sp3":    63.000,
"O_sb-C_sp3-O_sb":   160.000,"C_sp-C_sp3-O_sb":    50.000,"N_sp2-C_sp3-O_sb":    50.000,
"F-C_sp3-F":    77.000,"F-C_sp3-H":    50.000,"C_sp3-C_sp3-F":    50.000,"C_sp3-C_sp3-Cl":    50.000,
"Cl-C_sp3-H":    50.000,"Br-C_sp3-C_sp3":    50.000,"Br-C_sp3-H":    50.000,"C_sp3-C_sp3-I":    50.000,
"C_sp3-C_sp-N_sp":    80.000,"C_sp-C_sp-C_sp3":    80.000,"C_sp-C_sp-H":    50.000,
"C_sp2-N_sp2-C_sp3":    64.545,"C_sp2-N_sp2-H":    50.000,"C_sp3-N_sp2-H":    50.000,
"C_sp3-N_sp2-C_sp3":    50.000,"H-N_sp2-H":    35.000,"C_sp2-N_sp2-C_sp2":    70.000,
"C_sp3-N_sp3-H":    50.000,"C_sp3-N_sp3-C_sp3":    50.000,"H-N_sp3-H":    35.000,
"C_sp2-O_sb-H":    51.667,"C_sp3-O_sb-H":    55.000,"H-O_sb-P":    45.000,"C_sp2-O_sb-C_sp3":    60.000,
"C_sp3-O_sb-C_sp3":    60.000,"C_sp3-O_sb-P":   100.000,"C_sp2-O_sb-P":   100.000,"P-O_sb-P":   100.000,
"O_db-P-O_sb":    81.667,"O_db-P-O_db":   140.000,"O_sb-P-O_sb":    45.000,"C_sp3-S-C_sp3":    62.000,
"C_sp3-S-S":    68.000,"C_sp3-S-H":    43.000,"H-S-H":    35.000,"C_sp2-N_sp2-EP":   150.000,
"C_sp3-N_sp3-EP":   150.000,"EP-N_sp3-H":   150.000,"C_sp2-O_db-EP":   150.000,"EP-O_db-EP":   150.000,
"C_sp2-O_sb-EP":   150.000,"C_sp3-O_sb-EP":   150.000,"EP-O_sb-H":   150.000,"EP-O_sb-EP":   150.000,
"C_sp3-S-EP":   150.000,"EP-O_sb-P":   150.000,"EP-S-EP":   150.000,"EP-S-H":   150.000,
"EP-S-S":   150.000,"O_sb":   103.958,"H":     0.000,"C_sp2":    62.307,"N_sp2":    74.741,"P":    82.238,
"C_sp3":    54.763,"C_sp":    70.000,"N_sp3":    80.833,"S":   110.800,"O_db":   150.000}

def GetSimplifiedAtomType(atom):
    for key,value in BASIC_ATOM_TYPES.items():
        if atom in value:
            return str(key)
        
def CalcBondLength(bond_string):
    [a1,a2] = bond_string.split("-")
    a1 = GetSimplifiedAtomType(a1)
    a2 = GetSimplifiedAtomType(a2)
    [a1,a2] = sorted([a1,a2])
    base_bond = f"{a1}-{a2}"
    if base_bond in BOND_LENGTH_DICT.keys():
        return BOND_LENGTH_DICT[base_bond]
    return 1.500

def CalcBondForceConstant(bond_string):
    [a1,a2] = bond_string.split("-")
    a1 = GetSimplifiedAtomType(a1)
    a2 = GetSimplifiedAtomType(a2)
    [a1,a2] = sorted([a1,a2])
    base_bond = f"{a1}-{a2}"
    if base_bond in BOND_FORCE_CONST_DICT.keys():
        return BOND_FORCE_CONST_DICT[base_bond]
    return 300.000

def CalcAngleMeasure(angle):
    [a1,a2,a3] = angle.split("-")
    a1 = GetSimplifiedAtomType(a1)
    a2 = GetSimplifiedAtomType(a2)
    a3 = GetSimplifiedAtomType(a3)
    [a1,a3] = sorted([a1,a3])
    angle_key = f"{a1}-{a2}-{a3}"
    if angle_key in ANGLE_MEASURE_DICT.keys():
        return ANGLE_MEASURE_DICT[angle_key]
    elif a2 in ANGLE_MEASURE_DICT.keys():
        return ANGLE_MEASURE_DICT[a2]
    else:
        return 120.000

def CalcAngleForceConstant(angle):
    [a1,a2,a3] = angle.split("-")
    a1 = GetSimplifiedAtomType(a1)
    a2 = GetSimplifiedAtomType(a2)
    a3 = GetSimplifiedAtomType(a3)
    [a1,a3] = sorted([a1,a3])
    angle_key = f"{a1}-{a2}-{a3}"
    if angle_key in ANGLE_FORCE_CONST_DICT.keys():
        return ANGLE_FORCE_CONST_DICT[angle_key]
    elif a2 in ANGLE_FORCE_CONST_DICT.keys():
        return ANGLE_FORCE_CONST_DICT[a2]
    else:
        return 50.0

def NewBond(bond_string):
    bond_length = CalcBondLength(bond_string)
    force_const = CalcBondForceConstant(bond_string)
    return f"{bond_string:<5}  {force_const:>6.2f}   {bond_length:>5.3f}\n"

def NewAngle(angle_string):
    [a1,a2,a3]=angle_string.split("-")
    angle_measure = CalcAngleMeasure(angle_string)
    force_const = CalcAngleForceConstant(angle_string)
    return f"{angle_string:<8}   {force_const:>5.1f}      {angle_measure:>6.3f}"

def NewDihedral(dihedral_string):
    [a1,a2,a3,a4]=dihedral_string.split("-")

def NewTorsion(torsion_string):
    [a1,a2,a3,a4]=torsion_string.split("-")
