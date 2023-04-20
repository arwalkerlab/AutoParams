import sys
sys.path.append("/opt/conda/envs/psi4flask/lib/python3.9/site-packages/")

from utilities.utilities import *
from utilities.main_job import *

app = Flask(__name__)
app.config['SECRET_KEY'] = "01123581321_Dockerized"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['DATABASE_FOLDER'] = DATABASE_DIR
main_start_dir = os.getcwd()
MAIN_UPLOAD_FOLDER = os.path.join(main_start_dir,"uploads/")
RefreshDB()

### Starting Page
@app.route('/',methods=['GET', 'POST'])
def start_page():
    if request.method=="POST":
        jobid = random_job_identifier()
        session["jobid"] = jobid
        CURRENT_JOBS[jobid] = MainJob(jobid)
        CURRENT_JOBS[jobid]._charge = request.form.get('charge')
        CURRENT_JOBS[jobid]._multiplicity = request.form.get('spin')
        CURRENT_JOBS[jobid]._restype = request.form.get('moltype')
        if request.form.get('optimize_job') == "on":
            CURRENT_JOBS[jobid]._opt_complete = False
        else:
            CURRENT_JOBS[jobid]._opt_complete = True
        print("Current optimization boolean: ",CURRENT_JOBS[jobid]._opt_complete)
        continue_job = CURRENT_JOBS[jobid].UploadPDBFile(request.files['PDBfile'])
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
                print(frcmod)
                print(mol2)
                print(prmtop)
                print(inpcrd)
                print(tleapinput)
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
    orig_logfile = CURRENT_JOBS[session["jobid"]].file_list["JobLog"]
    new_logfile = f"logfiles/{session['jobid']}.html"
    new_param_set = f"logfiles/params_{session['jobid']}.html"
    with open("templates/" + new_param_set,"w") as f:
        frcmod = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["FRCMOD"]
        mol2 = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["MOL2"]
        prmtop = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["PRMTOP"]
        inpcrd = UPLOAD_FOLDER + CURRENT_JOBS[session["jobid"]].file_list["INPCRD"]
        tleapinput = UPLOAD_FOLDER + session["jobid"] + "/tleap.in"
        print(frcmod)
        print(mol2)
        print(prmtop)
        print(inpcrd)
        print(tleapinput)
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