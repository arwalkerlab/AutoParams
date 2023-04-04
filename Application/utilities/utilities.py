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