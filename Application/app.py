import sys
sys.path.append("/opt/conda/envs/psi4flask/lib/python3.9/site-packages/")

from utilities.utilities import *
from utilities.main_job import *
from utilities.initialize_server import InitializeApp

app = Flask(__name__, instance_path=DATABASE_DIR)

# Clean startup from utilities.initialize_server.py
InitializeApp(app)

### Starting Page
@app.route('/',methods=['GET', 'POST'])
def start_page():
    if request.method=="POST":
        jobid = random_job_identifier()
        session["jobid"] = jobid
        ### Collect Form Data
        pdbfile = request.files['PDBfile']
        optimize_bool = bool(request.form.get('optimize_job') != "on")
        override_db_bool = bool(request.form.get('override_db') == "on")
        charge = request.form.get('charge')
        multiplicity = request.form.get('spin')
        restype = request.form.get('moltype')
        level_of_theory = request.form.get('level_of_theory')
        basis_set = request.form.get('basis_set')
        headconnect = request.form.get('headconnect')
        tailconnect = request.form.get('tailconnect')
        cap_atoms = request.form.get('capatoms')
        if cap_atoms != "NONE":
            cap_atoms = [x.strip() for x in cap_atoms.split(",")]
        else:
            cap_atoms = []
        connections = []
        if headconnect == "NONE":
            connections.append("0")
        else:
            connections.append(headconnect)
        if tailconnect == "NONE":
            connections.append("0")
        else:
            connections.append(tailconnect)
        ### Initialize Job
        CURRENT_JOBS[jobid] = MainJob(jobid,charge,multiplicity,restype,optimize_bool,override_db_bool,level_of_theory,basis_set,connections,cap_atoms)

        ### Upload PDB and begin processing.
        continue_job = CURRENT_JOBS[jobid].UploadPDBFile(pdbfile)
        if not continue_job:
            new_logfile = f"logfiles/{session['jobid']}.html"
            S.call(f"touch templates/{new_logfile}",shell=True)
            new_param_set = f"logfiles/params_{session['jobid']}.html"
            with open("templates/" + new_param_set,"w") as f:
                frcmod = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["FRCMOD"]
                mol2 = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["MOL2"]
                prmtop = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["PRMTOP"]
                inpcrd = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["INPCRD"]
                tleapinput = UPLOAD_FOLDER + session["jobid"] + "/tleap.in"
                if G(frcmod):
                    f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['FRCMOD']) }}">FRCMOD</a><br>""")
                if G(mol2):
                    f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['MOL2']) }}">MOL2</a><br>""")
                if G(prmtop):
                    f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['PRMTOP']) }}">PRMTOP</a><br>""")
                if G(inpcrd):
                    f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['INPCRD']) }}">INPCRD</a><br>""")
                if G(tleapinput):
                    f.write("""<a href="{{ url_for('download', filename='"""+tleapinput+"""' }}">INPCRD File</a><br>""")

            return render_template("finished.html", curr_job = CURRENT_JOBS[session['jobid']] ,logfile=new_logfile, paramset = new_param_set)
        return render_template("jobqueue.html",curr_job=CURRENT_JOBS[jobid])
    ###  make a folder in 'uploads/' using that job-identifier.
    ###  render_template('loading.html')
    return render_template("upload_page.html")

### Job Running Page
@app.route('/process_files/',methods=['GET', 'POST'])
def process_files():
    continue_job = CURRENT_JOBS[session["jobid"]].CheckPDBQuality()
    if not continue_job:
        CURRENT_JOBS[session["jobid"]].LogJobMessage("PDB of insufficient quality.")
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].Optimize()
    if not continue_job:
        CURRENT_JOBS[session["jobid"]].LogJobMessage("Unable to optimize structure.")
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].RESPCharges()
    if not continue_job:
        CURRENT_JOBS[session["jobid"]].LogJobMessage("Unable to generate RESP charges.")
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].Parametrize()
    if not continue_job:
        CURRENT_JOBS[session["jobid"]].LogJobMessage("Unable to generate parameters.")
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].TestParams()
    if not continue_job:
        CURRENT_JOBS[session["jobid"]].LogJobMessage("Unable to successfully test parameters.")
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].AddResultsToDB()
    if not continue_job:
        CURRENT_JOBS[session["jobid"]].LogJobMessage("Unable to add parameters to database.")
        return "Error Encountered."
    CURRENT_JOBS[session["jobid"]].LogJobMessage("Parametrization Complete!")
    return "Files Processed"

### Job Finished Page
@app.route('/finished',methods=['GET', 'POST'])
def show_finished():
    orig_logfile = CURRENT_JOBS[session["jobid"]].file_list["JobLog"]
    new_logfile = f"logfiles/{session['jobid']}.html"
    new_param_set = f"logfiles/params_{session['jobid']}.html"
    with open("templates/" + new_param_set,"w") as f:
        frcmod = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["FRCMOD"]
        mol2 = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["MOL2"]
        prmtop = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["PRMTOP"]
        inpcrd = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["INPCRD"]
        tleapinput = UPLOAD_FOLDER + session["jobid"] + "/tleap.in"
        if G(frcmod):
            f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['FRCMOD']) }}">FRCMOD</a><br>""")
        if G(mol2):
            f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['MOL2']) }}">MOL2</a><br>""")
        if G(prmtop):
            f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['PRMTOP']) }}">PRMTOP</a><br>""")
        if G(inpcrd):
            f.write("""<a href="{{ url_for('download', filename=curr_job.file_list['INPCRD']) }}">INPCRD</a><br>""")
        if G(tleapinput):
            f.write("""<a href="{{ url_for('download', filename='"""+tleapinput+"""' }}">INPCRD File</a><br>""")

    S.call(f"cp {orig_logfile} templates/{new_logfile}",shell=True)
   
    return render_template("finished.html", curr_job = CURRENT_JOBS[session['jobid']] ,logfile=new_logfile,paramset = new_param_set)

### Get Download Link
@app.route('/upload/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename, as_attachment=True)

@app.route('/database/<path:filename>', methods=['GET', 'POST'])
def db_download(filename):
    return send_from_directory(app.config['DATABASE_FOLDER'], filename, as_attachment=True)

### Database Page
@app.route('/database',methods=["GET","POST"])
def database():
    RefreshDB()
    return render_template("database.html")