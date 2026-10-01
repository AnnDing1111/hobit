# -*- coding: utf-8 -*-
"""
Flask front end for HobIT: a local web page wrapping the same pipeline the
command line runs. Each submission gets its own task folder under results/,
the analysis runs in a background thread, and the browser polls /status until
a marker file appears.

Here define:
    class HobitArgs
        __init__(task_dir, mode, params, genome_path=None, query_path=None,
                 protein_path=None, hmm_path=None, ref_gff_path=None)
    manage_session()
    inject_translations()
    index()
    toggle_lang()
    run_analysis()
    status(task_id)
    download(task_id, filename)
    _save_upload(task_dir, key, filename)
    _run_task(task_dir, mode, args_obj, genome_path, protein_path, use_esm)
    _bundled_data_file(filename)
    _run_esm_scoring(task_dir, args_obj)
    find_free_port()
    open_browser(port)
    main()

@author: Yurou Ding <yd54@st-andrews.ac.uk>
"""

import os
import socket
import subprocess
import threading
import time
import traceback
import uuid
import webbrowser
from importlib import resources
from pathlib import Path
from flask import (Flask, jsonify, render_template, request, send_file,
                   session, url_for)

from hobit.direct_hmm import run_protein_pipeline
from hobit.hobit_main import process_single_genome

app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET_KEY", "hobit_secret_" + str(uuid.uuid4()))

UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

# inputs and marker files share the task folder with the outputs,
# so these are kept out of the download list
_INPUT_FILENAMES = {
    'completed.txt', 'error.txt', 'ai_warning.txt',
    'genome.fasta', 'query.fasta', 'proteins.fasta',
    'reference.gff3', 'model.hmm'
}

TRANSLATIONS = {
    'en': {
        'title': 'HobIT: Homeobox Identification Tool',
        'tab_genome': 'Genome Mode',
        'tab_protein': 'Protein Mode',
        'tab_help': 'Help & Manual',
        'lbl_genome': 'Genome Sequence (.fasta)',
        'lbl_query': 'Query Reference DB (.fasta) [Optional]',
        'lbl_gff': 'Reference Annotation (.gff3/gtf) [Optional]',
        'lbl_protein': 'Protein Sequence (.fasta)',
        'lbl_hmm': 'HMM Model (.hmm) [Optional]',
        'lbl_params': 'Analysis Parameters',
        'btn_run': 'RUN ANALYSIS',
        'status_running': 'Processing... Please wait.',
        'status_done': 'Analysis Complete!',
        'status_error': 'Error occurred.',
        'help_intro': 'HobIT is a tool for identifying Homeobox genes.',
        'col_scaffold': 'Scaffold Name',
        'col_family': 'Predicted Family',
        'col_score': 'Bit Score'
    },
    'zh': {
        'title': 'HobIT: 同源异型盒基因鉴定工具',
        'tab_genome': '基因组模式',
        'tab_protein': '蛋白质模式',
        'tab_help': '帮助与说明',
        'lbl_genome': '基因组序列 (.fasta)',
        'lbl_query': '查询序列 (.fasta) [可选]',
        'lbl_gff': '参考注释 GFF (.gff3/gtf) [可选]',
        'lbl_protein': '蛋白质序列 (.fasta)',
        'lbl_hmm': 'HMM 模型 (.hmm) [可选]',
        'lbl_params': '分析参数',
        'btn_run': '开始分析',
        'status_running': '正在处理……请耐心等待',
        'status_done': '分析完成！',
        'status_error': '发生错误',
        'help_intro': 'HobIT是一个用于鉴定Homeobox基因的工具。',
        'col_scaffold': '基因组Scaffold名称',
        'col_family': '预测家族',
        'col_score': '比对得分'
    }
}


# THE ARGUMENT OBJECT HANDED TO THE PIPELINE

class HobitArgs:
    # stand-in for the argparse Namespace that hobit_main builds; every attribute
    # either mode might read is set here, since the pipeline reads them directly

    def __init__(self, task_dir, mode, params, genome_path=None, query_path=None,
                 protein_path=None, hmm_path=None, ref_gff_path=None):
        self.out = str(task_dir)
        self.mode = mode
        self.threads = params['threads']
        self.max_evalue = params['max_evalue']
        self.verbose = True
        self.quiet = False

        self.window = params['window']
        self.cluster = params['cluster']
        self.exonerate = params['exonerate'] if mode == 'genome' else False
        self.exonerate_model = 'protein2genome'
        self.filter_repeats = params['filter_repeats']

        self.genome = str(genome_path) if genome_path else None
        self.query = str(query_path) if query_path else ""
        self.protein = str(protein_path) if protein_path else None
        self.hmm = str(hmm_path) if hmm_path else None
        self.compare_annotation = str(ref_gff_path) if ref_gff_path else None
        self.input_files = [self.genome] if mode == 'genome' else [self.protein]

        # the GUI has no separate switches for these, so they follow --exonerate
        self.gff = self.exonerate
        self.out_exons = self.exonerate
        self.out_introns = self.exonerate
        self.fasta_out = ['protein', 'cds', 'gene'] if self.exonerate else ['protein']

        self.pickle = None
        self.xml = None
        self.no_mkblastdb = False
        self.no_xlsx = False
        self.mafft = False


# FLASK ROUTES

@app.before_request
def manage_session():
    if 'lang' not in session:
        session['lang'] = 'en'


@app.context_processor
def inject_translations():
    return {'t': TRANSLATIONS[session['lang']], 'lang': session['lang']}


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/toggle_lang')
def toggle_lang():
    session['lang'] = 'zh' if session['lang'] == 'en' else 'en'
    return jsonify({'status': 'ok'})


@app.route('/run', methods=['POST'])
def run_analysis():
    # save the uploads, start the analysis in a background thread, return its ID
    task_id = str(uuid.uuid4())
    task_dir = Path(RESULTS_FOLDER) / task_id
    task_dir.mkdir(parents=True, exist_ok=True)

    mode = request.form.get('mode', 'genome')
    use_esm = request.form.get('use_esm') == 'true'

    genome_path = _save_upload(task_dir, 'genome_file', 'genome.fasta')
    query_path = _save_upload(task_dir, 'query_file', 'query.fasta')
    ref_gff_path = _save_upload(task_dir, 'ref_gff_file', 'reference.gff3')
    protein_path = _save_upload(task_dir, 'protein_file', 'proteins.fasta')
    hmm_path = _save_upload(task_dir, 'hmm_file', 'model.hmm')

    # form values arrive as strings and are entirely user-controlled
    try:
        params = {
            'window': int(request.form.get('window', 5000)),
            'cluster': int(request.form.get('cluster', 200)),
            'max_evalue': float(request.form.get('max_evalue', 0.1)),
            'threads': int(request.form.get('threads', 4)),
            'exonerate': request.form.get('exonerate') == 'true',
            'filter_repeats': request.form.get('filter_repeats') == 'true'
        }
    except ValueError:
        return jsonify({'task_id': None, 'msg': 'Invalid parameters'}), 400

    # fall back to the bundled data files when the user uploaded none
    if not query_path and mode == 'genome':
        query_path = _bundled_data_file('all_HDs.fa')
    if not hmm_path and mode == 'protein':
        hmm_path = _bundled_data_file('HD.hmm')

    args_obj = HobitArgs(task_dir, mode, params, genome_path, query_path,
                         protein_path, hmm_path, ref_gff_path)

    threading.Thread(
        target=_run_task,
        args=(task_dir, mode, args_obj, genome_path, protein_path, use_esm)
    ).start()

    return jsonify({'task_id': task_id})


@app.route('/status/<task_id>')
def status(task_id):
    # report running / complete / error, listing the downloadable outputs
    task_dir = Path(RESULTS_FOLDER) / task_id

    if (task_dir / "error.txt").exists():
        return jsonify({'status': 'error', 'msg': (task_dir / "error.txt").read_text()})

    if not (task_dir / "completed.txt").exists():
        return jsonify({'status': 'running'})

    # walk the whole task folder: genome mode nests its outputs one level deeper
    files = []
    for root, _, filenames in os.walk(task_dir):
        for fname in filenames:
            if fname in _INPUT_FILENAMES:
                continue

            rel_path = os.path.relpath(os.path.join(root, fname), task_dir)
            rel_path = rel_path.replace('\\', '/')
            download_url = url_for('download', task_id=task_id, filename=rel_path)

            files.append({'name': fname, 'url': download_url})

    response = {'status': 'complete', 'files': files}

    ai_warning_path = task_dir / "ai_warning.txt"
    if ai_warning_path.exists():
        response['ai_warning'] = ai_warning_path.read_text()

    return jsonify(response)


@app.route('/download/<task_id>/<path:filename>')
def download(task_id, filename):
    return send_file(os.path.join(RESULTS_FOLDER, task_id, filename), as_attachment=True)


# HELPERS
def _save_upload(task_dir, key, filename):
    # save one uploaded file into the task folder; None if the field was left empty
    uploaded = request.files.get(key)

    if uploaded and uploaded.filename:
        dest = task_dir / filename
        uploaded.save(str(dest))
        return dest
    return None


def _run_task(task_dir, mode, args_obj, genome_path, protein_path, use_esm):
    # runs in the background thread; everything is caught and written to error.txt,
    # because this thread has nowhere to raise to and the browser only sees marker files
    try:
        if mode == 'genome':
            if not genome_path:
                raise ValueError("Genome file is required for genome mode.")
            process_single_genome(args_obj, str(genome_path))
        else:
            if not protein_path:
                raise ValueError("Protein file is required for protein mode.")
            run_protein_pipeline(args_obj, str(protein_path))

        if use_esm:
            _run_esm_scoring(task_dir, args_obj)

        (task_dir / "completed.txt").touch()

    except Exception:
        err_msg = traceback.format_exc()
        print(err_msg)
        (task_dir / "error.txt").write_text(err_msg)


def _bundled_data_file(filename):
    # locate a file in hobit/data, falling back to a __file__-relative path
    try:
        return resources.files('hobit.data').joinpath(filename)
    except Exception:
        fallback = Path(__file__).parent / 'data' / filename
        return fallback if fallback.exists() else None


def _run_esm_scoring(task_dir, args_obj):
    # optional ESM2 plugin: its absence or failure is a warning in ai_warning.txt,
    # not a failed analysis, since the HobIT results are fine without it
    candidates = list(task_dir.rglob("hobit_protein.fasta"))

    if not candidates:
        msg = "ESM2 AI scoring skipped: hobit_protein.fasta not found in results."
        print(f">>> {msg}")
        (task_dir / "ai_warning.txt").write_text(msg)
        return

    print(f">>> ESM2 AI scoring: found candidate file {candidates[0].name}")

    try:
        from hobit.esm_plugin import run_ai_scoring
        run_ai_scoring(str(candidates[0].parent), args_obj.query)
    except Exception:
        msg = ("AI scoring plugin (ESM2) encountered an error and was skipped.\n\n"
               + traceback.format_exc())
        print(f"WARNING: {msg}")
        (task_dir / "ai_warning.txt").write_text(msg)


# SERVER STARTUP

def find_free_port():
    # prefer 5000-5009 (macOS AirPlay often takes 5000), then let the OS pick one
    for port in range(5000, 5010):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', port)) != 0:
                return port

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        return s.getsockname()[1]


def open_browser(port):
    # wait for Flask to bind, then open the interface
    time.sleep(1.5)
    url = f"http://127.0.0.1:{port}"
    print(f"\n{'-' * 50}")
    print(f"Web Interface Ready: {url}")
    print(f"{'-' * 50}\n")

    # under WSL there is no Linux browser, so the URL is handed to Windows
    try:
        with open('/proc/version') as f:
            if 'microsoft' in f.read().lower():
                print(">>> WSL environment detected - invoking Windows browser.")
                subprocess.run(['explorer.exe', url])
                return
    except OSError:
        pass

    print(">>> Opening default browser...")
    try:
        webbrowser.open(url)
    except webbrowser.Error:
        print(f"Could not open browser automatically. Please navigate to: {url}")


def main():
    # entry point for the hobit-gui command
    port = find_free_port()
    print(f">>> Port automatically configured to: {port}")
    threading.Thread(target=open_browser, args=(port,)).start()
    app.run(host='0.0.0.0', port=port, debug=False, use_reloader=False)


if __name__ == '__main__':
    main()
