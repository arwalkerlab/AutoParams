BASIC_ATOM_TYPES   = {"C_sp":["CZ","CY","cg","ch","c1"],
                      "C_sp2":["C4","C5","C ","CA","CB","CC","CD","CK","CP","CM","CS","CN","CQ","CR","CV","CW","C*","CI","Ck","Cj","cC","cB","c ","ce","cf","ca"],
                      "C_sp3":["CX","CT","cA","cD","cP","cR","Cg","Cy","2C","3C","Cp","cx","cy"],
                      "N_sp":["NY","n1","no"],
                      "N_sp2":["N ","NA","NB","NC","N2","N*","nN","Ng","nb","nc","nd","n ","ni","nj","ne","nf","n2"],
                      "N_sp3":["N3","NT","nA","n4","nk","nl","np","n3","nq","nh","nm","nn","na"],
                      "O_sb":["OW","OH","OS","oH","oR","oS","oT","Os","Oy","Oh","ow","oh","os","oq"],
                      "O_db":["OP","O2","O ","oC","oO","oP","o "],
                      "H":["HW","H4","H5","HA","HC","H1","H2","H3","HP","HZ","H ","HO","HS","hA","hL","hE","hS","hX","hB","hN","hO","hR","Hc","Hp","Ha","Ho","hc","h1","h2","h3","hx","ha","h4","h5","hn","ho","hs","hw","hp"],
                      "F":["F ","f "],
                      "Cl":["Cl","cl","CL"],
                      "Br":["Br","br","BR"],
                      "I":["I ","i "],
                      "S":["S ","SH","SO","SS","Sm","SF","S*","s ","s2","s4","s6","sh","ss","sp","sq","sx","sy"],
                      "P":["P ","pA","p2","p3","p4","p5","pb","pc","pd","pe","pf","px","py"],
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
                     'EP-N_sp2':0.200, 'EP-S':0.700}

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

DIHEDRAL_SET_DICT = {"C_sp2-C_sp2-C_sp2-O_db":"WW-XX-YY-ZZ   1    2.175       180.0            -2.\nWW-XX-YY-ZZ   1    0.30          0.0             3.",
                   "C_sp3-C_sp2-C_sp2-C_sp3":"WW-XX-YY-ZZ   1    6.65        180.0            -2.\nWW-XX-YY-ZZ   1    1.90        180.0             1.",
                   "C_sp2-C_sp2-C_sp2-X ":"WW-XX-YY-ZZ   1    2.175       180.0            -2.\nWW-XX-YY-ZZ   1    0.30          0.0             3.",
                   "O_db-C_sp2-C_sp2-X ":"WW-XX-YY-ZZ   1    2.175       180.0            -2.\nWW-XX-YY-ZZ   1    0.30          0.0             3.",
                   "C_sp3-C_sp2-C_sp2-X ":"WW-XX-YY-ZZ   1    6.65        180.0            -2.\nWW-XX-YY-ZZ   1    1.90        180.0             1.",
                   "X -C_sp2-C_sp2-X ":"WW-XX-YY-ZZ    4   16.10        180.0             2. ",
                   "H -C_sp2-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.00          0.0            -3.\nWW-XX-YY-ZZ   1    0.250         0.0             1.",
                   "O_sb-C_sp2-C_sp3-C_sp3":"WW-XX-YY-ZZ   1    1.178040    190.97653        -1. \nWW-XX-YY-ZZ   1    0.092102    295.63279        -2. \nWW-XX-YY-ZZ   1    0.962830    348.09535         3. ",
                   "C_sp2-C_sp2-C_sp3-H ":"WW-XX-YY-ZZ   1    0.38        180.0            -3. \nWW-XX-YY-ZZ   1    1.15          0.0             1. ",
                   "N_sp2-C_sp2-C_sp3-N_sp2":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.55        180.0            -3. \nWW-XX-YY-ZZ   1    1.58        180.0            -2.\nWW-XX-YY-ZZ   1    0.45        180.0             1.",
                   "N_sp2-C_sp2-C_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3. \nWW-XX-YY-ZZ   1    0.20          0.0            -2.\nWW-XX-YY-ZZ   1    0.20          0.0             1.",
                   "O_db-C_sp2-C_sp3-H ":"WW-XX-YY-ZZ   1    0.80          0.0            -1. \nWW-XX-YY-ZZ   1    0.00          0.0            -2. \nWW-XX-YY-ZZ   1    0.08        180.0             3. ",
                   "X -C_sp2-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.00          0.0            -3.\nWW-XX-YY-ZZ   1    0.250         0.0             1.",
                   "X -C_sp2-C_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3. \nWW-XX-YY-ZZ   1    0.20          0.0            -2.\nWW-XX-YY-ZZ   1    0.20          0.0             1.",
                   "X -C_sp2-C_sp3-H ":"WW-XX-YY-ZZ   1    0.80          0.0            -1. \nWW-XX-YY-ZZ   1    0.00          0.0            -2. \nWW-XX-YY-ZZ   1    0.08        180.0             3. ",
                   "X -C_sp2-C_sp3-N_sp2":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.55        180.0            -3. \nWW-XX-YY-ZZ   1    1.58        180.0            -2.\nWW-XX-YY-ZZ   1    0.45        180.0             1.",
                   "H -C_sp2-C_sp3-X ":"WW-XX-YY-ZZ   1    0.00          0.0            -3.\nWW-XX-YY-ZZ   1    0.250         0.0             1.",
                   "O_sb-C_sp2-C_sp3-X ":"WW-XX-YY-ZZ   1    1.178040    190.97653        -1. \nWW-XX-YY-ZZ   1    0.092102    295.63279        -2. \nWW-XX-YY-ZZ   1    0.962830    348.09535         3. ",
                   "C_sp2-C_sp2-C_sp3-X ":"WW-XX-YY-ZZ   1    0.38        180.0            -3. \nWW-XX-YY-ZZ   1    1.15          0.0             1. ",
                   "N_sp2-C_sp2-C_sp3-X ":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.55        180.0            -3. \nWW-XX-YY-ZZ   1    1.58        180.0            -2.\nWW-XX-YY-ZZ   1    0.45        180.0             1.",
                   "X -C_sp2-C_sp3-X ":"WW-XX-YY-ZZ   6    0.00          0.0             2. ",
                   "O_db-C_sp2-N_sp2-H ":"WW-XX-YY-ZZ   1    2.50        180.0            -2. \nWW-XX-YY-ZZ   1    2.00          0.0             1. ",
                   "O_db-C_sp2-N_sp2-EP":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "C_sp3-C_sp2-N_sp2-EP":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "X -C_sp2-N_sp2-H ":"WW-XX-YY-ZZ   1    2.50        180.000          -2. \nWW-XX-YY-ZZ   1    2.00          0.0             1. ",
                   "X -C_sp2-N_sp2-EP":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "O_db-C_sp2-N_sp2-X ":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "C_sp3-C_sp2-N_sp2-X ":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "X -C_sp2-N_sp2-X ":"WW-XX-YY-ZZ   2    7.40        180.000           2. ",
                   "X -C_sp2-O_db-X ":"WW-XX-YY-ZZ   4   11.20        180.000           2. ",
                   "O_db-C_sp2-O_sb-H ":"WW-XX-YY-ZZ   1    2.30        180.000          -2.\nWW-XX-YY-ZZ   1    1.90          0.0             1.",
                   "O_db-C_sp2-O_sb-C_sp3":"WW-XX-YY-ZZ   1    2.70        180.000          -2.\nWW-XX-YY-ZZ   1    1.40        180.0             1.",
                   "C_sp2-C_sp2-O_sb-EP":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "H -C_sp2-O_sb-EP":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "O_db-C_sp2-O_sb-X ":"WW-XX-YY-ZZ   1    2.30        180.000          -2.\nWW-XX-YY-ZZ   1    1.90          0.0             1.",
                   "C_sp2-C_sp2-O_sb-X ":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "H -C_sp2-O_sb-X ":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "X -C_sp2-O_sb-H ":"WW-XX-YY-ZZ   1    2.30        180.000          -2.\nWW-XX-YY-ZZ   1    1.90          0.0             1.",
                   "X -C_sp2-O_sb-C_sp3":"WW-XX-YY-ZZ   1    2.70        180.000          -2.\nWW-XX-YY-ZZ   1    1.40        180.0             1.",
                   "X -C_sp2-O_sb-EP":"WW-XX-YY-ZZ   1    0.000       180.000           2.0",
                   "X -C_sp2-O_sb-X ":"WW-XX-YY-ZZ   2    4.00        180.000           2.",
                   "X -C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   9    1.40          0.0             3.",
                   "X -C_sp3-C_sp-X ":"WW-XX-YY-ZZ   3    0.00          0.0             1. ",
                   "O_sb-C_sp3-N_sp2-C_sp2":"WW-XX-YY-ZZ   1    0.98023     110.0984         -1.\nWW-XX-YY-ZZ   1    1.38071      13.7765         -2.\nWW-XX-YY-ZZ   1    0.60481     176.3635         -3.\nWW-XX-YY-ZZ   1    0.30675      17.8068          4.",
                   "C_sp2-C_sp3-N_sp2-C_sp2":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.42          0.0            -3.\nWW-XX-YY-ZZ   1    0.27          0.0            -2.\nWW-XX-YY-ZZ   1    0.00          0.0             1.",
                   "C_sp3-C_sp3-N_sp2-C_sp2":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3.\nWW-XX-YY-ZZ   1    2.00          0.0            -2.\nWW-XX-YY-ZZ   1    2.00          0.0             1.",
                   "X -C_sp3-N_sp2-C_sp2":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.42          0.0            -3.\nWW-XX-YY-ZZ   1    0.27          0.0            -2.\nWW-XX-YY-ZZ   1    0.00          0.0             1.",
                   "O_sb-C_sp3-N_sp2-X ":"WW-XX-YY-ZZ   1    0.98023     110.0984         -1.\nWW-XX-YY-ZZ   1    1.38071      13.7765         -2.\nWW-XX-YY-ZZ   1    0.60481     176.3635         -3.\nWW-XX-YY-ZZ   1    0.30675      17.8068          4.",
                   "C_sp2-C_sp3-N_sp2-X ":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.42          0.0            -3.\nWW-XX-YY-ZZ   1    0.27          0.0            -2.\nWW-XX-YY-ZZ   1    0.00          0.0             1.",
                   "C_sp3-C_sp3-N_sp2-X ":"WW-XX-YY-ZZ   1    0.00          0.0            -4. \nWW-XX-YY-ZZ   1    0.40          0.0            -3.\nWW-XX-YY-ZZ   1    2.00          0.0            -2.\nWW-XX-YY-ZZ   1    2.00          0.0             1.",
                   "X -C_sp3-N_sp2-X ":"WW-XX-YY-ZZ   6    0.00          0.0             2. ",
                   "C_sp3-C_sp3-N_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.30          0.0            -3.\nWW-XX-YY-ZZ   1    0.48        180.0             2.",
                   "H -C_sp3-N_sp3-EP":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "C_sp3-C_sp3-N_sp3-EP":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "X -C_sp3-N_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.30          0.0            -3.\nWW-XX-YY-ZZ   1    0.48        180.0             2.",
                   "X -C_sp3-N_sp3-EP":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "C_sp3-C_sp3-N_sp3-X ":"WW-XX-YY-ZZ   1    0.30          0.0            -3.\nWW-XX-YY-ZZ   1    0.48        180.0             2.",
                   "H -C_sp3-N_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "X -C_sp3-N_sp3-X ":"WW-XX-YY-ZZ    6    1.80          0.0             3. ",
                   "C_sp2-C_sp3-O_sb-C_sp3":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.100       180.0             2.",
                   "C_sp3-C_sp3-O_sb-H ":"WW-XX-YY-ZZ   1    0.16          0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1.",
                   "C_sp3-C_sp3-O_sb-C_sp3":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.1         180.0             2.",
                   "C_sp3-C_sp3-O_sb-C_sp2":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.80        180.0             1.",
                   "C_sp3-C_sp3-O_sb-EP":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "O_sb-C_sp3-O_sb-C_sp3":"WW-XX-YY-ZZ   1    0.10          0.0            -3.\nWW-XX-YY-ZZ   1    0.85        180.0            -2.\nWW-XX-YY-ZZ   1    1.35        180.0             1.",
                   "N_sp2-C_sp3-O_sb-C_sp3":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.65          0.0             2.",
                   "H -C_sp3-O_sb-EP":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "X -C_sp3-O_sb-H ":"WW-XX-YY-ZZ   1    0.16          0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1.",
                   "X -C_sp3-O_sb-C_sp3":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.1         180.0             2.",
                   "X -C_sp3-O_sb-C_sp2":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.80        180.0             1.",
                   "X -C_sp3-O_sb-EP":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "C_sp2-C_sp3-O_sb-X ":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.100       180.0             2.",
                   "C_sp3-C_sp3-O_sb-X ":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.80        180.0             1.",
                   "O_sb-C_sp3-O_sb-X ":"WW-XX-YY-ZZ   1    0.10          0.0            -3.\nWW-XX-YY-ZZ   1    0.85        180.0            -2.\nWW-XX-YY-ZZ   1    1.35        180.0             1.",
                   "N_sp2-C_sp3-O_sb-X ":"WW-XX-YY-ZZ   1    0.383         0.0            -3.\nWW-XX-YY-ZZ   1    0.65          0.0             2.",
                   "H -C_sp3-O_sb-X ":"WW-XX-YY-ZZ   1    0.000         0.000           3.0",
                   "X -C_sp3-O_sb-X ":"WW-XX-YY-ZZ   3    0.50          0.0             3. ",
                   "X -C_sp3-S -X ":"WW-XX-YY-ZZ   3    1.00          0.0             3. ",
                   "C_sp2-O_sb-P -O_db":"WW-XX-YY-ZZ   1    0.0         180.0            -3.\nWW-XX-YY-ZZ   1    0.0         180.0             1.",
                   "C_sp2-O_sb-P -O_sb":"WW-XX-YY-ZZ   1    0.185181     31.79508        -1.\nWW-XX-YY-ZZ   1    1.256531    351.95960        -2.\nWW-XX-YY-ZZ   1    0.354858    357.24748         3.",
                   "C_sp3-O_sb-P -O_sb":"WW-XX-YY-ZZ   1    0.25          0.0            -3.\nWW-XX-YY-ZZ   1    1.20          0.0             2.",
                   "X -O_sb-P -O_db":"WW-XX-YY-ZZ   1    0.0         180.0            -3.\nWW-XX-YY-ZZ   1    0.0         180.0             1.",
                   "X -O_sb-P -O_sb":"WW-XX-YY-ZZ   1    0.185181     31.79508        -1.\nWW-XX-YY-ZZ   1    1.256531    351.95960        -2.\nWW-XX-YY-ZZ   1    0.354858    357.24748         3.",
                   "C_sp2-O_sb-P -X ":"WW-XX-YY-ZZ   1    0.185181     31.79508        -1.\nWW-XX-YY-ZZ   1    1.256531    351.95960        -2.\nWW-XX-YY-ZZ   1    0.354858    357.24748         3.",
                   "C_sp3-O_sb-P -X ":"WW-XX-YY-ZZ   1    0.25          0.0            -3.\nWW-XX-YY-ZZ   1    1.20          0.0             2.",
                   "X -O_sb-P -X ":"WW-XX-YY-ZZ   3    0.75          0.0             3.",
                   "C_sp2-C_sp3-C_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1.",
                   "H -C_sp3-C_sp3-H ":"WW-XX-YY-ZZ   1    0.15          0.0             3. ",
                   "H -C_sp3-C_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.16          0.0             3.",
                   "O_sb-C_sp3-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.144         0.0            -3.\nWW-XX-YY-ZZ   1    1.175         0.0             2.",
                   "F -C_sp3-C_sp3-F ":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    1.20        180.0             1.",
                   "Cl-C_sp3-C_sp3-Cl":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.45        180.0             1.",
                   "Br-C_sp3-C_sp3-Br":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.00        180.0             1.",
                   "H -C_sp3-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1.",
                   "H -C_sp3-C_sp3-F ":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.19          0.0             1.",
                   "H -C_sp3-C_sp3-Cl":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1.",
                   "H -C_sp3-C_sp3-Br":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.55          0.0             1.",
                   "H -C_sp3-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.25          0.0             1.",
                   "N_sp2-C_sp3-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4.",
                   "C_sp2-C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1.",
                   "H -C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.16          0.0             3.",
                   "F -C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    1.20        180.0             1.",
                   "Cl-C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.45        180.0             1.",
                   "Br-C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.00        180.0             1.",
                   "N_sp2-C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4.",
                   "C_sp3-C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1.",
                   "O_sb-C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4.",
                   "X -C_sp3-C_sp3-C_sp2":"WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1.",
                   "X -C_sp3-C_sp3-H ":"WW-XX-YY-ZZ   1    0.16          0.0             3.",
                   "X -C_sp3-C_sp3-F ":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    1.20        180.0             1.",
                   "X -C_sp3-C_sp3-Cl":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.45        180.0             1.",
                   "X -C_sp3-C_sp3-Br":"WW-XX-YY-ZZ   1    0.000         0.0            -3.\nWW-XX-YY-ZZ   1    0.00        180.0             1.",
                   "X -C_sp3-C_sp3-N_sp2":"WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4.",
                   "X -C_sp3-C_sp3-C_sp3":"WW-XX-YY-ZZ   1    0.180         0.0            -3.\nWW-XX-YY-ZZ   1    0.250       180.0            -2.\nWW-XX-YY-ZZ   1    0.200       180.0             1.",
                   "X -C_sp3-C_sp3-O_sb":"WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4.",
                   "X -C_sp3-C_sp3-X ":"WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4.",
                   "C_sp2-N_sp2-P -O_db":"WW-XX-YY-ZZ   1    0.12          0.0             3.",
                   "X -N_sp2-P -O_db":"WW-XX-YY-ZZ   1    0.12          0.0             3.",
                   "C_sp2-N_sp2-P -X ":"WW-XX-YY-ZZ   1    0.12          0.0             3.",
                   "X -N_sp2-P -X ":"WW-XX-YY-ZZ   1    0.12          0.0             3.",
                   "C_sp3-S -S -C_sp3":"WW-XX-YY-ZZ   1    3.50          0.0            -2.\nWW-XX-YY-ZZ   1    0.60          0.0             3.",
                   "C_sp3-S -S -EP":"WW-XX-YY-ZZ   1    0.00          0.0             3.",
                   "C_sp3-S -S -X ":"WW-XX-YY-ZZ   1    3.50          0.0            -2.\nWW-XX-YY-ZZ   1    0.60          0.0             3.",
                   "X -S -S -C_sp3":"WW-XX-YY-ZZ   1    3.50          0.0            -2.\nWW-XX-YY-ZZ   1    0.60          0.0             3.",
                   "X -S -S -X ":"WW-XX-YY-ZZ   1    0.00          0.0             3.",
                   "C_sp3-C_sp-C_sp-H ":"WW-XX-YY-ZZ   1    0.00          0.0             1.",
                   "X -C_sp-C_sp-H ":"WW-XX-YY-ZZ   1    0.00          0.0             1.",
                   "C_sp3-C_sp-C_sp-X ":"WW-XX-YY-ZZ   1    0.00          0.0             1.",
                   "X -C_sp-C_sp-X ":"WW-XX-YY-ZZ   1    0.00          0.0             1.",
                  }


def GetSimplifiedAtomType(atom):
    for key,value in BASIC_ATOM_TYPES.items():
        if atom in value:
            return str(key)
    return "X "
        
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
    
def NewMass(mass_string):
    return f"{mass_string:<2} 999.999\n"

def NewBond(bond_string):
    bond_length = CalcBondLength(bond_string)
    force_const = CalcBondForceConstant(bond_string)
    return f"{bond_string:<5}  {force_const:>6.2f}   {bond_length:>5.3f}\n"

def NewAngle(angle_string):
    [a1,a2,a3]=angle_string.split("-")
    angle_measure = CalcAngleMeasure(angle_string)
    force_const = CalcAngleForceConstant(angle_string)
    return f"{angle_string:<8}   {force_const:>5.1f}      {angle_measure:>6.3f}\n"

def NewDihedral(dihe_string):
    [a1,a2,a3,a4] = dihe_string.split("-")
    a1 = GetSimplifiedAtomType(a1)
    a2 = GetSimplifiedAtomType(a2)
    a3 = GetSimplifiedAtomType(a3)
    a4 = GetSimplifiedAtomType(a4)
    if [a3,a2] == sorted([a2,a3]):
        [a1,a2,a3,a4] = [a4,a3,a2,a1]
        
    ## check for complete key in dihedral dictionary
    if f"{a1:<2}-{a2:<2}-{a3:<2}-{a4:<2}" in DIHEDRAL_SET_DICT.keys():
        dihedral_string = DIHEDRAL_SET_DICT[f"{a1:<2}-{a2:<2}-{a3:<2}-{a4:<2}"]
       
    ## check for first three in dihedral dictionary
    elif f"{a1:<2}-{a2:<2}-{a3:<2}-X " in DIHEDRAL_SET_DICT.keys():
        dihedral_string = DIHEDRAL_SET_DICT[f"{a1:<2}-{a2:<2}-{a3:<2}-X "]
    
    ## check for last three in dihedral dictionary
    elif f"X -{a2:<2}-{a3:<2}-{a4:<2}" in DIHEDRAL_SET_DICT.keys():
        dihedral_string = DIHEDRAL_SET_DICT[f"X -{a2:<2}-{a3:<2}-{a4:<2}"]
    
    ## check for middle two in dihedral dictionary
    elif f"X -{a2:<2}-{a3:<2}-X " in DIHEDRAL_SET_DICT.keys():
        dihedral_string = DIHEDRAL_SET_DICT[f"X -{a2:<2}-{a3:<2}-X "]
    
    ## return a basic dihedral parameter. (sp3-C to sp3-C single bond)
    else:
        dihedral_string = "WW-XX-YY-ZZ   1    0.000         0.000          -1.\nWW-XX-YY-ZZ   1    1.490         0.000          -2.\nWW-XX-YY-ZZ   1    0.156         0.000          -3.\nWW-XX-YY-ZZ   1    0.000         0.000           4."
    [a1,a2,a3,a4] = dihe_string.split("-")
    dihedral_string = dihedral_string.replace("WW",a1)
    dihedral_string = dihedral_string.replace("XX",a2)
    dihedral_string = dihedral_string.replace("YY",a3)
    dihedral_string = dihedral_string.replace("ZZ",a4)
    dihedral_string = dihedral_string.strip() + "\n"
    return dihedral_string

def NewTorsion(torsion_string):
    [a1,a2,a3,a4]=torsion_string.split("-")

def NewNonBonded(nonbonded_string):
    return f"{nonbonded_string:<2}          2.00    0.5000\n"