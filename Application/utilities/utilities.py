from flask import Flask, flash, request, redirect, url_for, render_template, session, send_from_directory,send_file
from werkzeug.utils import secure_filename
from glob import glob as G
import os
import time
import subprocess as S
import uuid
import parmed

### Import all my default settings
from .defaults import *

def random_job_identifier():
    id = uuid.uuid4()
    return id.hex

def RenewSmilesDB():
    S.call(f"touch {SMILES_DB}",shell=True)
    os.chdir(DATABASE_DIR)
    for folder in G("*/"):
        smi_file = G(folder+"*.smi")[0]
        with open(SMILES_DB,"a") as f:
            f.write(open(smi_file).read())
            f.write("\n")
    tmp_lines = []
    for line in open(SMILES_DB).readlines():
        if line.strip() == "":
            continue
        if line not in tmp_lines:
            tmp_lines.append(line)
    with open(SMILES_DB,"w") as f:
        for line in tmp_lines:
            f.write(line)
    os.chdir(MAIN_DIR)

def RefreshDB():
    global DATABASE_MOLS_SEEN
    DATABASE_MOLS_SEEN = []
    RenewSmilesDB()
    os.chdir(DATABASE_DIR)
    with open(TEMPLATES_DIR+"__db_dataset.html","w") as f:
        all_db_folders = G("*/")
        all_db_folders.sort(key=lambda x: os.path.getmtime(x))
        for file_loc in all_db_folders:
            f.write(DataBaseEntry(file_loc))
    CleanUploadsFolder()
    os.chdir(MAIN_DIR)

def DataBaseEntry(location):
    global DATABASE_MOLS_SEEN
    html_code = """<div class="grid-item"><br>"""
    smiles = ""
    if G(location+"*.smi"):
        smiles = open(G(location+'*.smi')[0]).read().strip()
    else:
        return ""
    if smiles in DATABASE_MOLS_SEEN:
        return ""
    DATABASE_MOLS_SEEN.append(smiles)
    if G(location+"*.png"):
        chemdraw = G(location+"*.png")[0]
        html_code+=f"<img src='database/{chemdraw}' width=\"98%\" title=\"SMILES: {smiles}\"><br>\n"
    if G(location+"settings.txt"):
        html_code += "<a href=\"{{ url_for('db_download', filename=\""
        html_code +=G(location+"settings.txt")[0]
        html_code +="\")}}\">AutoParams Settings File</a><br>\n"
    if G(location+"*.frcmod"):
        html_code += "<a href=\"{{ url_for('db_download', filename=\""
        html_code +=G(location+"*.frcmod")[0]
        html_code +="\")}}\">FRCMOD</a><br>\n"
    if G(location+"*.mol2"):
        html_code += "<a href=\"{{ url_for('db_download', filename=\""
        html_code +=G(location+"*.mol2")[0]
        html_code +="\")}}\">MOL2</a><br>\n"
    if G(location+"tleap.in"):
        html_code += "<a href=\"{{ url_for('db_download', filename=\""
        html_code +=G(location+"tleap.in")[0]
        html_code +="\")}}\">tleap Input File</a><br>\n"
    if G(location+"*.prmtop"):
        html_code += "<a href=\"{{ url_for('db_download', filename=\""
        html_code +=G(location+"*.prmtop")[0]
        html_code +="\")}}\">PRMTOP</a><br>\n"
    if G(location+"*.inpcrd"):
        html_code += "<a href=\"{{ url_for('db_download', filename=\""
        html_code +=G(location+"*.inpcrd")[0]
        html_code +="\")}}\">INPCRD</a><br>\n"
    html_code+="</div>\n"
    return html_code

def CheckSMILESinDB(SMILES):
    if not G(SMILES_DB):
        S.call(f"touch {SMILES_DB}",shell=True)
    lines = open(SMILES_DB,"r").read().split("\n")
    if SMILES in lines:
        return True
    return False

def MaybeCopy(source,dest):
    if G(source):
        S.call(f"cp {source} {dest}",shell=True)

def FileExists(filepath):
    return os.path.exists(filepath)

def CleanUploadsFolder():
    def file_age_in_days(filename):
        st = os.stat(filename)
        age_in_seconds = (time.time() - st.st_mtime)
        return int(age_in_seconds/86400)
    os.chdir(MAIN_DIR+'/uploads/')
    for directory in G("*/"):
        if file_age_in_days(directory) > 6:
            S.call(f"rm -r {directory}",shell=True)
        