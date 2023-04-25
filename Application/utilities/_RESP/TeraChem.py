import subprocess as S
import os
from glob import glob as G

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
                    "run":"energy",
                    "resp":"yes",
}

def WriteTCRESP(**kwargs):
    respsettings=Default_Kwargs
    for key,value in kwargs.items():
        respsettings[key]=value
    with open("resp.in","w") as f:
        for key in respsettings.keys():
            f.write(f"{key:<20}{respsettings[key]:<}\n")
    return

def readRESPcharges(filename):
    ## TERACHEM-specific read of RESP charges
    charges = []
    i = 0
    with open(filename,'r') as f:
        lines = f.readlines()
        for line in lines:
            if "ESP restraint charges:" in line:
                break
            if "ESP unrestraint charges:" in line:
                break
            i = i + 1
        for line in lines[i+3:]:
            if "-----------------------------------------------------------" in line:
                break
            charge = float(line[38:48])
            charges.append(str(charge))
    return charges

def GetRESPCharges(pdbfile,charge,mult,jobfolder,level_of_theory="b3lyp",basis_set="6-31gss"):
    ## Terachem Version
    currdir = os.getcwd()
    os.chdir(jobfolder)
    WriteTCRESP(coordinates=pdbfile,charge=charge,spinmult=mult,method=level_of_theory,basis=basis_set)
    S.call("terachem -i resp.in 1> resp.out 2> resp.err",shell=True)
    resp_charges = readRESPcharges("resp.out")
    os.chdir(currdir)
    return resp_charges