import sys
sys.path.append("/opt/conda/envs/psi4flask/lib/python3.9/site-packages/")
from utilities.utilities import *
from utilities.main_job import *

app = Flask(__name__)
app.config['SECRET_KEY'] = "01123581321_Dockerized"
app.config['UPLOAD_FOLDER'] = 'uploads/'
main_start_dir = os.getcwd()
MAIN_UPLOAD_FOLDER = os.path.join(main_start_dir,"uploads/")

### Starting Page
@app.route('/',methods=['GET', 'POST'])
def start_page():
    if request.method=="POST":
        jobid = random_job_identifier()
        session["jobid"] = jobid
        CURRENT_JOBS[jobid] = MainJob(jobid)
        if bool(request.form.get('optimize_job')):
            CURRENT_JOBS[jobid]._opt_complete = False
        else:
            CURRENT_JOBS[jobid]._opt_complete = True
        CURRENT_JOBS[jobid].UploadPDBFile(request.files['PDBfile'])
        return render_template("jobqueue.html")
    ###  make a folder in 'uploads/' using that job-identifier.
    ###  render_template('loading.html')
    return render_template("upload_page.html")

### Job Running Page
@app.route('/process_files/',methods=['GET', 'POST'])
def process_files():
    continue_job = CURRENT_JOBS[session["jobid"]].CheckPDBQuality()
    if not continue_job:
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].Optimize()
    if not continue_job:
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].RESPCharges()
    if not continue_job:
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].Parametrize()
    if not continue_job:
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].TestParams()
    if not continue_job:
        return "Error Encountered."
    continue_job = CURRENT_JOBS[session["jobid"]].AddResultsToDB()
    if not continue_job:
        return "Error Encountered."
    return "Files Processed"

### Job Finished Page
@app.route('/finished',methods=['GET', 'POST'])
def show_finished():
    pdb = CURRENT_JOBS[session["jobid"]].file_list["Working PDB"].split("/")[-1]
    charge = CURRENT_JOBS[session["jobid"]]._charge
    moltype = CURRENT_JOBS[session["jobid"]]._restype
    resname = CURRENT_JOBS[session["jobid"]]._resname
    frcmod = CURRENT_JOBS[session["jobid"]].file_list["FRCMOD"]
    mol2 = CURRENT_JOBS[session["jobid"]].file_list["MOL2"]
    imagefile = CURRENT_JOBS[session["jobid"]].file_list["ChemDraw"].split("/")[-1]
    
    return render_template("finished.html",pdb=pdb,
            charge=charge,
            moltype=moltype,
            resid=resname,
            frcmod = frcmod,
            mol2 = mol2,
            imagename=imagefile)

### Get Download Link
@app.route('/upload/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename, as_attachment=True)

### Database Page
@app.route('/database',methods=["GET","POST"])
def database():
    return render_template("database.html")

# @app.route('/',methods=['GET', 'POST'])
# def start_page():
#     os.chdir(main_start_dir)
#     global spiff_job
#     if request.method == 'POST':
#         if 'file' not in request.files:
#             flash('No file part')
#             return redirect(request.url)
#         file = request.files['file']
#         if file.filename == '':
#             flash('No selected file')
#             return redirect(request.url)
#         if file and allowed_file(file.filename):
#             filename = secure_filename(file.filename)
#             jobnum = len(G(app.config['UPLOAD_FOLDER']+"*/"))+1
#             jobdir = os.path.join(app.config['UPLOAD_FOLDER'], f"JOB_{jobnum:>08}/") 
#             S.call(f"mkdir -p {jobdir}",shell=True)
#             pdbsavepath = os.path.join(jobdir, filename)
#             file.save(pdbsavepath)
#             imagename = filename.replace(".pdb",".png")
#             PDBtoChemDraw(pdbsavepath,f"static/{imagename}")
#             charge = request.form['charge']
#             moltype = request.form['moltype']
#             resid = request.form['resid']
#             spiff_job = SPIFF_Job(f"JOB_{jobnum:>08}/{filename}",charge,moltype,resid)
#             return render_template('loading.html',
#             pdb=filename,
#             charge=charge,
#             moltype=moltype,
#             resid=resid,
#             imagename=imagename
#             )
#     return render_template('upload_page.html')

# @app.route('/finished',methods=['GET', 'POST'])
# def show_finished():
#     global spiff_job
#     n_db_dirs=len(G(f"{main_start_dir}/database/*/"))+1
#     dbdir = f"{main_start_dir}/database/job_{n_db_dirs:>09}/"
#     os.makedirs(dbdir,exist_ok=True)
#     S.call(f"cp {spiff_job.frcmod} {spiff_job.mol2} {dbdir}",shell=True)
#     os.chdir(main_start_dir)
#     return render_template("finished.html",pdb=spiff_job.pdb,
#             charge=spiff_job.charge,
#             moltype=spiff_job.moltype,
#             resid=spiff_job.resid,
#             frcmod = spiff_job.frcmod,
#             mol2 = spiff_job.mol2,
#             imagename=spiff_job.pdb.replace(".pdb",".png"))

# @app.route('/process_files/',methods=['GET', 'POST'])
# def process_files():
#     global spiff_job
#     currdir = os.getcwd()
#     os.chdir(app.config['UPLOAD_FOLDER'])
#     spiff_job.run()
#     os.chdir(currdir)
#     return "Files Processed"


