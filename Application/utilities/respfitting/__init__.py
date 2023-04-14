import subprocess as S
import os
from glob import glob as G

def GetRESPCharges(pdbfile,charge,mult,jobfolder):
    ## This version uses Psi4/PsiRESP with a custom script included alongside
    ## the main Flask Application to properly run PsiRESP externally.
    resp_charges = []
    S.call(f"PsiRESPJob -p {pdbfile} -c {charge} -m {mult} -d {jobfolder}",shell=True)
    currdir = os.getcwd()
    os.chdir(jobfolder+"/single_point/")
    S.call("sh run_single_point.sh",shell=True)
    os.chdir(currdir)
    S.call(f"PsiRESPJob -p {pdbfile} -c {charge} -m {mult} -d {jobfolder}",shell=True)
    if not G("resp.out"):
        return False
    for line in open("resp.out","r").readlines():
        resp_charges.append(line.strip())
    return resp_charges