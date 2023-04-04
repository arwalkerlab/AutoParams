from flask import Flask, flash, request, redirect, url_for, render_template, session, send_from_directory
from werkzeug.utils import secure_filename
from glob import glob as G
import os
import subprocess as S
import uuid
import parmed

### Import all my default settings
from .defaults import *

def random_job_identifier():
    id = uuid.uuid4()
    return id.hex

def DataBaseEntry(location):
    chemdraw = G(location+"*.png")[0]
    frcmod   = G(location+"*.frcmod")[0]
    mol2     = G(location+"*.mol2")[0]
    tleap    = G(location+"tleap.in")[0]
    prmtop   = G(location+"*.prmtop")[0]
    inpcrd   = G(location+"*.inpcrd")[0]
    settings = G(location+"settings.txt")[0]
    smiles   = open(G(location+"*.smi")[0]).read().split()[0]
    html     = """
<fieldset style="border-width:2px; border-style:solid; border-color:#78E2A0; padding: 1em;">
<legend>"""+smiles+"""</legend>
<img src='database/"""+chemdraw+"""' width="500" }}">
<a href="{{ url_for('db_download', filename="""+frcmod+""") }}">frcmod File</a><br>
<a href="{{ url_for('db_download', filename="""+mol2+""" }}">mol2 File</a><br>
</fieldset>
    """
