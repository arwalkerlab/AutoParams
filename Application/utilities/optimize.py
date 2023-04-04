import os
from glob import glob as G
import subprocess as S

def OptimToMainPDB(pdbfile):
    workdir = os.path.dirname(os.path.abspath(pdbfile))
    scrdir = os.path.abspath(G(os.path.join(workdir,"scr.*/"))[0])
    optimpdb = os.path.abspath(G(os.path.join(scrdir,"optim.pdb"))[0])
    with open(optimpdb) as f:
        lines = f.readlines()
        start = 0
        end   = 0
        for i,line in enumerate(lines):
            if "MODEL" in line:
                start = i
            if "ENDMDL" in line:
                end = i
        newpdb = "".join(x for x in lines[start+1:end])    
    with open(pdbfile,"w") as f:
        f.write(newpdb)

def TCOpt(folder,pdbfile,charge):
    startdir = os.getcwd()
    os.chdir(folder)
    pdbfile = G("*.pdb")[0]
    with open("tc.in","w") as f:
        f.write("""#Input files
coordinates     """+str(pdbfile)+"""

# QM
basis           6-31gss
method          b3lyp
convthre        1e-7
threall         1e-14
scf             diis+a
charge          """+str(charge)+"""
spinmult        1
purify          no

#Runtype
run                  minimize
new_minimizer        no
min_coordinates      cartesian

#Computer
gpus            1
gpumem          256
""")
    S.call("terachem -i tc.in 1> opt.out 2> opt.err",shell=True)
    os.chdir(startdir)