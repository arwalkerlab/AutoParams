#!/usr/bin/env python
from ..defaults import *
import os
import subprocess as S


Default_Kwargs = {  "coordinates":"input.xyz",
                    "charge":"0",
                    "spinmult":"1",
                    "basis":"6-31gss",
                    "method":"b3lyp",
                    "convthre":"1e-7",
                    "threall":"1e-14",
                    "precision":"mixed",
                    "maxit":"200",
                    "scf":"diis+a",
                    "gpus":"1",
                    "gpumem":"256",
                    "scrdir":"scr/",
                    "run":"minimize",
                    "new_minimizer":"no",
                    "min_coordinates":"cartesian"}

def WriteTCInput(pdbfile,**kwargs):
    workpath = os.path.dirname(os.path.abspath(pdbfile))
    tcinput = os.path.join(workpath,"opt.in")
    maxkeylen = 20
    settings = Default_Kwargs
    settings["coordinates"] = pdbfile
    for key in settings.keys():
        if len(key) > maxkeylen:
            maxkeylen = len(key)
    for key,val in kwargs.items():
        settings[key] = str(val)
    with open(tcinput,"w") as f:
        for key,val in settings.items():
            f.write(f"{key:<{maxkeylen}} {val}\n")

def OptimizePDB(pdbfile,charge=0,mult=1,method="b3lyp/6-31gss"):
    currdir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(pdbfile)))
    WriteTCInput(pdbfile,charge=charge,spinmult=mult,method=method.split('/')[0],basis=method.split('/')[1])
    with open("TC_Opt.sh","w") as f:
        f.write("""#!/bin/bash
terachem -i opt.in 1> opt.out 2> opt.err
if grep -q "Job finished:" opt.out; then
    OPTIMPDB=$(ls scr/optim.pdb)
    line1=$(grep -n "MODEL" $OPTIMPDB | tail -n 1 | grep -Eo '^[^:]+')
    line2=$(grep -n "ENDMDL" $OPTIMPDB | tail -n 1 | grep -Eo '^[^:]+')
    tail -n +$((line1+1)) $OPTIMPDB | head -n $((line2-line1-1)) > $MAINPDB
    rm -rf ./scr/
fi
""")
    S.call("sh TC_Opt.sh",shell=True)

    os.chdir(currdir)
    return 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--pdbfile",dest="pdbfile",required=True)
    parser.add_argument("-c","--charge",dest="charge",default=0)
    parser.add_argument("-s","--spin",dest="spin",default=1)
    parser.add_argument("-m","--method",dest="method",default="scf/6-31g**")
    args = parser.parse_args()
    OptimizePDB(args.pdbfile,args.charge,args.spin,args.method)