#!/bin/env python3

# -*- coding: utf-8 -*-

import os
from datetime import datetime

current_file = __file__

last_modified_timestamp = os.path.getmtime(current_file)

last_modified_datetime = datetime.fromtimestamp(last_modified_timestamp)
# print("Last modified datetime of this script:", last_modified_datetime)

import random
import time
import traceback
import warnings
import sys
import subprocess
from Bio import SeqIO
from datetime import datetime
import argparse
import shutil
import os
import math
from concurrent.futures import ThreadPoolExecutor

def warn_with_traceback(msg, cat, fname, lno, file=None, line=None):
    """
    Custom warning handler that prints a traceback with the warning message.
    """
    logfh = sys.stderr if file is None else file
    traceback.print_stack(file=logfh)
    logfh.write(warnings.formatwarning(msg, cat, fname, lno, line))
    logfh.flush()
    
def get_default_cpu_count():
    """
    Detects logical threads and returns half, minimum 1.
    Prioritizes affinity-aware counts for container/HPC compatibility.
    """
    try:
        # Respects CPU affinity (important for Docker/Slurm)
        logical_threads = len(os.sched_getaffinity(0))
    except AttributeError:
        # Fallback for Windows/macOS
        logical_threads = os.cpu_count() or 1
    
    return max(1, math.floor(logical_threads / 2)), logical_threads

warnings.showwarning = warn_with_traceback

start_time = datetime.now()
cwd_path = os.path.abspath(os.getcwd())
params = dict()
DEFAULT_CPU, TOTAL_THREADS = get_default_cpu_count()

VERSION = "1.0"

HELP_MSG = 'smith_sa v{} - Sequence Mapper for Insertion boundaries in Target Hosts – Sequence Alignment\n'.format(VERSION)
HELP_MSG += '(c) 2026. Pola, G.L. & Arthur Gruber\n'
HELP_MSG += 'For the latest version access: https://github.com/gruberlab/smith_sa\n'
HELP_MSG += 'Usage: smith_sa.py -q <file name> -d <file name> -run local\n'
HELP_MSG += '\tsmith_sa.py -q <file name> -d <file name> -run local -tab <file name>\n'
HELP_MSG += '\tsmith_sa.py -q <file name> -run remote -tab <file name>\n'
HELP_MSG += '\tsmith_sa.py -q <file name> -run remote\n'
HELP_MSG += '\nMandatory parameters:\n'
HELP_MSG += '-q \t<file name>: DNA sequence containing the inserted element\n'
HELP_MSG += '-d \t<file name>: BLAST database (only valid when using -run local)\n'
HELP_MSG += '-run \t<local|remote>: whether BLAST search should run locally or remotely\n'
HELP_MSG += '\nOptional parameters:\n'
HELP_MSG += '-cpu \t<integer>: Number of threads to execute the local BLASTN search (default: {} [half of {} detected threads])\n'.format(DEFAULT_CPU, TOTAL_THREADS)
HELP_MSG += '-color \t<RGB number codes>: The color of the element to be used as a qualifier in the annotation file, separated by commas (default: 255,0,0)\n'
HELP_MSG += '-enddist \t<integer>: Maximum allowed distance (in base pairs, bp) between the 5' or 3' end of an alignment block and the query terminus (default: 50)\n'
HELP_MSG += '-minlen \t<integer>: Minimum accepted length of the element (bp) (default: 4000)\n'
HELP_MSG += '-maxlen \t<integer>: Maximum accepted length of the element (bp) (default: 50000)\n'
HELP_MSG += '-mincov \t<integer>: Minimum % query coverage per subject (default: 30)\n'
HELP_MSG += '-maxcov \t<integer>: Maximum % query coverage per subject (default: 90)\n'
HELP_MSG += '-max_remote_proc \t<integer>: Maximum number of BLAST processes in a remote server (default: 1)\n'
HELP_MSG += '-max_batch_size \t<integer>: Maximum length of concatenated query sequences per batch for remote BLAST searches (default: 1000000)\n'
HELP_MSG += '-org \t<Taxonomy identifier>: TaxID(s) to restrict the database for BLASTN searches (when using -run remote)\n'
HELP_MSG += '-out \t<name>: Output directory name (default: output_dir1)\n'
HELP_MSG += '-tab \t<filename>: Filename of a previous result of a BLASTN search\n'

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-q')
parser.add_argument('-d')
parser.add_argument('-out')
parser.add_argument('-org')
parser.add_argument('-color')
parser.add_argument('-cpu')
parser.add_argument('-tab')
parser.add_argument('-enddist')
parser.add_argument('-minlen')
parser.add_argument('-maxlen')
parser.add_argument('-mincov')
parser.add_argument('-maxcov')
parser.add_argument('-run')
parser.add_argument('-max_remote_proc')
parser.add_argument('-web_inter_batch_delay')
parser.add_argument('-max_batch_size')
parser.add_argument('-version', action='store_true')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()

def is_fasta(fpath):
    """
    Checks if a file is a valid, non-empty FASTA file.
    Returns:
    bool: True if the file contains at least one valid FASTA record.
    """
    try:
        with open(fpath, "r") as handle:
            return any(SeqIO.parse(handle, "fasta"))
    except (IOError, ValueError):
        return False

def rename(i, name, typ):
    """
    Generates a unique file or directory name by appending an index if the
    original name exists.
    """
    path_checker = os.path.isdir if typ == 'dir' else os.path.isfile
    if not path_checker(name):
        return name
    base, ext = os.path.splitext(name)
    # Start numbering from 2 if the original exists.
    i = 2
    while True:
        new_name = "{}_{}{}".format(base, i, ext)
        if not path_checker(new_name):
            return new_name
        i += 1

def blast_parse(blast_fpath):
    """
    Parses a BLAST tabular output file with header comments.
    """
    if not os.path.exists(blast_fpath):
        sys.stderr.write("ERROR: BLAST results file not found: {}\n".format(blast_fpath))
        return [], [], []
    hits = []
    header_qids = []
    cols_str = ''
    try:
        with open(blast_fpath, "r") as f:
            for line in f.readlines():
                if not line.strip():
                    continue
                if line.startswith('#'):
                    if line.startswith('# Fields:'):
                        cols_str = line.split(':', 1)[1].strip()
                    elif line.startswith('# Query:'):
                        header_qids.append(line.split(':', 1)[1].strip())
                else:
                    hits.append(line.strip().split('\t'))
        parsed_cols = [col.strip() for col in cols_str.split(',')] if cols_str else []
        seen_qids = set()
        qids_ordered = [qid for qid in header_qids if qid not in seen_qids and not seen_qids.add(qid)]
        if not qids_ordered and hits:
            if log:
                log.write("\nNOTE: No '# Query:' headers found. Inferring order from hits.\n")
            for hit in hits:
                if hit[0] not in seen_qids:
                    qids_ordered.append(hit[0])
                    seen_qids.add(hit[0])
        return qids_ordered, parsed_cols, hits
    except IOError as e:
        sys.stderr.write("ERROR: Could not read BLAST file {}: {}\n".format(blast_fpath, e))
        return [], [], []

def validate_args(args):
    """
    Validates and processes all command-line arguments.
    """
    is_valid = True
    params = {}

    # --- Path and Mode Validation ---
    if not (args.q and os.path.isfile(args.q) and is_fasta(args.q)):
        sys.stderr.write("ERROR: A valid query FASTA file must be provided with -q.\n")
        is_valid = False
    else:
        params['q'] = os.path.realpath(args.q)

    if args.run and str(args.run).lower() in ['local', 'web']:
        params['run'] = str(args.run).lower()
    elif not args.tab:
        sys.stderr.write("ERROR: Must specify a run mode (-run local|web) or a BLAST table (-tab).\n")
        is_valid = False

    if params.get('run') == 'local' and not args.tab:
        if not (args.d):
            sys.stderr.write("ERROR: A valid database FASTA file is required with -d for local runs.\n")
            is_valid = False
        else:
            params['d'] = os.path.realpath(args.d)

    if args.tab:
        if not os.path.isfile(args.tab):
            sys.stderr.write("ERROR: BLAST table file (-tab) not found.\n")
            is_valid = False
        else:
            params['tab'] = os.path.realpath(args.tab)

    # --- Numeric Parameter Validation ---
    def validate_integer(arg_val, name, default, min_val=None, max_val=None):
        if arg_val is None:
            return default, True
        try:
            val = int(arg_val)
            if (min_val is not None and val < min_val) or \
               (max_val is not None and val > max_val):
                sys.stderr.write("ERROR: -{} must be between {} and {}.\n".format(name, min_val, max_val))
                return None, False
            return val, True
        except (ValueError, TypeError):
            sys.stderr.write("ERROR: -{} must be an integer.\n".format(name))
            return None, False

    params['cpu'], is_valid = validate_integer(args.cpu, 'cpu', DEFAULT_CPU, min_val=1)
    params['enddist'], is_valid = validate_integer(args.enddist, 'enddist', 50, min_val=0) if is_valid else (None, False)
    params['minlen'], is_valid = validate_integer(args.minlen, 'minlen', 4000, min_val=0) if is_valid else (None, False)
    params['maxlen'], is_valid = validate_integer(args.maxlen, 'maxlen', 50000, min_val=0) if is_valid else (None, False)
    params['mincov'], is_valid = validate_integer(args.mincov, 'mincov', 30, min_val=0, max_val=100) if is_valid else (None, False)
    params['maxcov'], is_valid = validate_integer(args.maxcov, 'maxcov', 90, min_val=0, max_val=100) if is_valid else (None, False)
    params['max_remote_proc'], is_valid = validate_integer(args.max_remote_proc, 'max_remote_proc', 1, min_val=1, max_val=10)
    params['max_batch_size'], is_valid = validate_integer(args.max_batch_size, 'max_batch_size', 1000000, min_val=10000)

    if is_valid and params['minlen'] >= params['maxlen']:
        sys.stderr.write("ERROR: -minlen must be less than -maxlen.\n")
        is_valid = False

    if is_valid and params['mincov'] >= params['maxcov']:
        sys.stderr.write("ERROR: -mincov must be less than -maxcov.\n")
        is_valid = False

    # --- Web-Specific and Other Parameters ---
    if params.get('run') == 'web':
        params['org'] = ' OR '.join(['txid{}[ORGN]'.format(tid.strip()) for tid in args.org.split(',')]) if args.org else None

    # Color validation
    if args.color is None:
        params['color'] = [255, 0, 0]
    else:
        try:
            colors = [int(c.strip()) for c in args.color.split(',')]
            if len(colors) == 3 and all(0 <= c <= 255 for c in colors):
                params['color'] = colors
            else:
                raise ValueError
        except (ValueError, TypeError):
            sys.stderr.write("ERROR: -color must be three comma-separated integers between 0 and 255.\n")
            is_valid = False

    # --- Output Directory Setup ---
    if is_valid:
        try:
            out_path_base = args.out if args.out else 'output_dir1'
            # Ensure path is absolute or relative to the original working directory
            out_path_base = os.path.join(cwd_path, out_path_base)
            final_out_path = rename(1, out_path_base, 'dir')
            params["out"] = final_out_path
            os.makedirs(final_out_path)
            print("Creating output directory (-out) {}...".format(final_out_path))
        except Exception as e:
            sys.stderr.write("ERROR: Could not create output directory: {}\n".format(e))
            is_valid = False

    return is_valid, params

def split_fasta(fasta_fpath):
    """
    Splits a FASTA file into memory-efficient batches based on total base pairs.
    """
    log.write("INFO: Starting to split FASTA file '{}' into batches...\n".format(fasta_fpath))
    MAX_BP_PER_BATCH = params.get('max_batch_size', 1000000)
    parts_dir = os.path.join(params["out"], "query_parts")
    log.write("INFO: Creating directory for split query files: '{}'\n".format(parts_dir))
    try:
        os.makedirs(parts_dir)
    except OSError as e:
        if e.errno != 17:  # errno 17 is 'File exists'
            raise

    out_files = []
    batch_num = 1
    current_batch_records = []
    current_batch_bp = 0

    try:
        for record in SeqIO.parse(fasta_fpath, "fasta"):
            record_len = len(record.seq)
            if record_len > MAX_BP_PER_BATCH:
                warn_msg = (
                    "WARN: Sequence '{}' ({} bp) exceeds MAX_BP_PER_BATCH ({}) and will be put in its own batch.\n"
                    .format(record.id, record_len, MAX_BP_PER_BATCH)
                )
                log.write(warn_msg)

            if current_batch_records and \
               (record_len > MAX_BP_PER_BATCH or current_batch_bp + record_len > MAX_BP_PER_BATCH):
                out_fpath = os.path.join(parts_dir, "query_part_{}.fasta".format(batch_num))
                log_msg = (
                    "INFO: Writing batch {} ({} sequences, {} bp) to file '{}'\n"
                    .format(batch_num, len(current_batch_records), current_batch_bp, out_fpath)
                )
                log.write(log_msg)
                SeqIO.write(current_batch_records, out_fpath, "fasta")
                out_files.append(out_fpath)
                batch_num += 1
                current_batch_records = []
                current_batch_bp = 0

            current_batch_records.append(record)
            current_batch_bp += record_len

        if current_batch_records:
            out_fpath = os.path.join(parts_dir, "query_part_{}.fasta".format(batch_num))
            log_msg = (
                "INFO: Writing final batch {} ({} sequences, {} bp) to file '{}'\n"
                .format(batch_num, len(current_batch_records), current_batch_bp, out_fpath)
            )
            log.write(log_msg)
            SeqIO.write(current_batch_records, out_fpath, "fasta")
            out_files.append(out_fpath)

        summary_msg = (
            "INFO: FASTA splitting complete. Created {} batch file(s) in '{}'.\n"
            .format(len(out_files), parts_dir)
        )
        log.write(summary_msg)
        return out_files
    except (IOError, ValueError) as e:
        err_msg = "ERROR: Failed during FASTA splitting for '{}': {}\n".format(fasta_fpath, e)
        sys.stderr.write(err_msg)
        if log:
            log.write(err_msg)
        raise

def build_blastn_cmd(blast_params, blast_args, n_threads):
    """
    Build the list of arguments to call blastn via subprocess.
    """
    cmd = ["blastn"]

    cmd.extend(["-query", blast_args["query"]])
    cmd.extend(["-out", blast_args["out"]])
    cmd.extend(["-outfmt", blast_args["outfmt"]])
    cmd.extend(["-task", blast_args.get("task", "megablast")])
    cmd.extend(["-max_target_seqs", str(blast_args.get("max_target_seqs", 100))])
    cmd.extend(["-evalue", str(blast_args.get("evalue", "1e-5"))])

    run_mode = blast_params.get("run")

    if run_mode == "web":
        cmd.extend(["-db", "nt"])
        cmd.append("-remote")
        if "org" in blast_params and blast_params["org"]:
            cmd.extend(["-entrez_query", blast_params["org"]])
    elif run_mode == "local":
        cmd.extend(["-db", blast_params.get("d")])
        if n_threads and n_threads > 1:
            cmd.extend(["-num_threads", str(n_threads)])
    else:
        raise RuntimeError("Invalid run type: {}".format(run_mode))

    return cmd

def run_blast_batch(batch_fpath, batch_idx, total_batches, blast_params, out_fpath=None):
    """
    Runs a BLAST search for a single batch, with retries and robust logging.
    """
    global blast_log_fh

    log_suffix = "batch {}/{}".format(batch_idx + 1, total_batches)

    if out_fpath:
        batch_out_fpath = out_fpath
    else:
        blast_parts_dir = blast_params.get(
            'current_blastn_parts_dir',
            os.path.join(blast_params["out"], "blastn_parts")
        )
        if not os.path.exists(blast_parts_dir):
            os.makedirs(blast_parts_dir)
        batch_out_fpath = os.path.join(
            blast_parts_dir,
            "blast_batch_{}.tab".format(batch_idx + 1)
        )

    blast_args = {
        'query': batch_fpath,
        "outfmt": "7 qseqid sseqid qcovs qlen slen qstart qend evalue bitscore",
        "out": batch_out_fpath,
        "task": "megablast",
        "max_target_seqs": 100,
        "evalue": "1e-5",
    }

    # threads only for local mode
    if blast_params.get('run') == 'local':
        n_threads = max(1, int(blast_params.get('cpu', 1)))
    else:
        n_threads = 1

    try:
        cmd_list = build_blastn_cmd(blast_params, blast_args, n_threads)
    except RuntimeError as e:
        err_msg = "ERROR: {} for {}.".format(e, log_suffix)
        sys.stderr.write(err_msg + "\n")
        if log:
            log.write("\n" + err_msg + "\n")
            log.flush()
        if blast_log_fh:
            blast_log_fh.write(err_msg + "\n")
            blast_log_fh.flush()
        return 'FAILURE', None, False

    max_attempts = 3
    base_backoff_delay = 15
    attempt = 0
    cpu_warn = False
    FATAL_ERR_PATTERNS = [
        "Failed to connect",
        "Connection refused",
        "Service not found",
        "Name or service not known",
        "failure in name resolution",
        "invalid output format",
        "Unknown output format",
    ]

    while attempt < max_attempts:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        cmd_string = " ".join(cmd_list)

        if blast_log_fh:
            blast_log_fh.write(
                "\n--- {} BLAST Attempt {}/{}, {} at {} ---\nCMD: {}\n".format(
                    blast_params.get('run', '').upper(),
                    attempt + 1,
                    max_attempts,
                    log_suffix,
                    ts,
                    cmd_string,
                )
            )
            blast_log_fh.flush()

        try:
            proc = subprocess.Popen(
                cmd_list,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = proc.communicate()
            retcode = proc.returncode

            if retcode != 0:
                raise subprocess.CalledProcessError(
                    retcode, cmd_list, output=stdout, stderr=stderr
                )

            if stderr and blast_log_fh:
                blast_log_fh.write("STDERR:\n{}\n".format(stderr))
                blast_log_fh.flush()
                if "consumed a large amount of server CPU time" in stderr:
                    cpu_warn = True
                    blast_log_fh.write("NCBI CPU usage warning detected.\n")
                    blast_log_fh.flush()

            if blast_log_fh:
                blast_log_fh.write(
                    "SUCCESS: BLAST for {} complete.\n".format(log_suffix)
                )
                blast_log_fh.flush()

            return 'SUCCESS', batch_out_fpath, cpu_warn

        except subprocess.CalledProcessError as e:
            err_detail = e.stderr or e.output or str(e)

            if (
                "BLAST query/options error" in err_detail
                or "Argument" in err_detail
                or any(p in err_detail for p in ["unknown output format", "invalid output format"])
            ):
                err_msg = (
                    "ERROR: [{}]: Permanent BLAST error detected. "
                    "Not retrying. See blastn.log.\n".format(log_suffix)
                )
                sys.stderr.write(err_msg)
                if log:
                    log.write("\n" + err_msg)
                    log.flush()
                if blast_log_fh:
                    blast_log_fh.write(err_msg)
                    blast_log_fh.flush()
                return 'FAILURE', None, cpu_warn

            if blast_log_fh:
                blast_log_fh.write(
                    "ERROR: BLAST for {} failed (attempt {}/{}), code {}.\n"
                    "Detail: {}\n".format(
                        log_suffix,
                        attempt + 1,
                        max_attempts,
                        e.returncode,
                        err_detail.strip(),
                    )
                )
                blast_log_fh.flush()

            attempt += 1
            if attempt < max_attempts:
                wait_time = (base_backoff_delay * (2 ** (attempt - 1))) + random.uniform(0, 5)
                retry_msg = (
                    "BLAST for [{}] failed. Retrying in {:.1f}s (attempt {}/{})..."
                    .format(log_suffix, wait_time, attempt + 1, max_attempts)
                )
                sys.stderr.write(retry_msg + "\n")
                if log:
                    log.write("\n" + retry_msg + "\n")
                if blast_log_fh:
                    blast_log_fh.write(retry_msg + "\n")
                    blast_log_fh.flush()
                if log:
                    log.flush()
                time.sleep(wait_time)
            else:
                final_fail_msg = (
                    "ERROR: [{}] BLAST failed after {} retries."
                    .format(log_suffix, max_attempts)
                )
                sys.stderr.write(final_fail_msg + "\n")
                if log:
                    log.write("\n" + final_fail_msg + "\n")
                if blast_log_fh:
                    blast_log_fh.write(final_fail_msg + "\n")
                    blast_log_fh.flush()
                return 'FAILURE', None, cpu_warn

    return 'FAILURE', None, cpu_warn

def run_blast(params, query_fpath, out_fname="blastn.tab", t_start=None):
    """
    Orchestrates the BLAST search, including splitting, execution, and merging.
    """
    final_blast_tab = os.path.join(params["out"], out_fname)
    query_parts_dir = os.path.join(params["out"], "query_parts")
    blast_parts_dir = os.path.join(params["out"], "blastn_parts")
    params['current_blastn_parts_dir'] = blast_parts_dir
    t_start = t_start or datetime.now()
    try:
        max_bp = params.get('max_batch_size', 1000000)
        total_bp = sum(len(rec.seq) for rec in SeqIO.parse(query_fpath, "fasta"))
        num_seqs = len(list(SeqIO.parse(query_fpath, "fasta")))
        if total_bp <= max_bp or num_seqs <= 1 or params.get('run') == 'local':
            if log:
                log.write("INFO: Single query or multiple query is within max batch size limits. Skipping splitting and using original query file.\n")
            status, out_path, _ = run_blast_batch(query_fpath, 0, 1, params, out_fpath=final_blast_tab)
            if status != 'SUCCESS':
                return datetime.now() - t_start, None
            return datetime.now() - t_start, out_path
        else:
            batch_fpaths = split_fasta(query_fpath)
            if not batch_fpaths:
                if log:
                    log.write("WARN: No FASTA batches were created from query.\n")
                return datetime.now() - t_start, None
            if not os.path.exists(blast_parts_dir):
                os.makedirs(blast_parts_dir)
            result_files = []
            total_batches = len(batch_fpaths)

            num_workers = params.get('max_remote_proc', 1)
            base_delay = params.get('web_inter_batch_delay', 60)

            effective_delay = base_delay * num_workers

            if params.get('run') == 'web' and num_workers > 1:
                log.write("INFO: Parallel Web BLAST enabled with {} workers.\n".format(num_workers))
                log.write("INFO: Scaled inter-batch delay set to {}s.\n".format(effective_delay))
                
                with ThreadPoolExecutor(max_workers=num_workers) as executor:
                    future_to_batch = {}
                    for i, batch_path in enumerate(batch_fpaths):
                        # Staggered submission
                        if i > 0:
                            time.sleep(effective_delay)
                        
                        future = executor.submit(
                            run_blast_batch, batch_path, i, total_batches, params
                        )
                        future_to_batch[future] = batch_path

                    # Collect results as they finish
                    for future in future_to_batch:
                        status, res_path, cpu_warn = future.result()
                        if status == 'SUCCESS':
                            result_files.append(res_path)
                        # Dynamic backoff logic if CPU warnings are detected
                        if cpu_warn:
                            params['web_inter_batch_delay'] = min(600, params['web_inter_batch_delay'] * 2)
            
            if not result_files:
                sys.stderr.write("ERROR: All BLAST batches failed.\n")
                return datetime.now() - t_start, None
            with open(final_blast_tab, 'wb') as outfile:
                for i, part_file in enumerate(sorted(result_files)):
                    with open(part_file, 'rb') as infile:
                        if i == 0:
                            shutil.copyfileobj(infile, outfile)
                        else:
                            for line in infile:
                                if not line.startswith(b'#'):
                                    outfile.write(line)
            return datetime.now() - t_start, final_blast_tab
    finally:
        if 'batch_fpaths' in locals():
            if os.path.isdir(query_parts_dir):
                shutil.rmtree(query_parts_dir)
            if os.path.isdir(blast_parts_dir):
                shutil.rmtree(blast_parts_dir)

def get_fasta_query_ids(fasta_fpath):
    """
    Extracts all unique sequence IDs from a FASTA file.
    """
    try:
        with open(fasta_fpath, "r") as handle:
            return {record.id for record in SeqIO.parse(handle, "fasta") if record.id}
    except (IOError, ValueError) as e:
        err_msg = "ERROR: Cannot parse query FASTA {}: {}\n".format(fasta_fpath, e)
        sys.stderr.write(err_msg)
        if log:
            log.write(err_msg)
        return set()

def get_missing_queries(blast_table, expected_qids):
    """
    Identifies query IDs from a set that are absent in a BLAST table.
    """
    if not os.path.exists(blast_table):
        if log:
            log.write("\nWARN: BLAST table '{}' not found.\n".format(blast_table))
        return list(expected_qids)
    found_qids = set()
    try:
        with open(blast_table, "r") as f:
            for line in f:
                sline = line.strip()
                if not sline:
                    continue
                if sline.startswith("# Query:"):
                    try:
                        found_qids.add(sline.split(":", 1)[1].strip())
                    except IndexError:
                        continue
                elif not sline.startswith("#"):
                    found_qids.add(sline.split("\t", 1)[0])
    except IOError as e:
        sys.stderr.write("ERROR reading BLAST table {}: {}\n".format(blast_table, e))
        return list(expected_qids)
    missing_ids = expected_qids - found_qids
    return sorted(list(missing_ids)) if missing_ids else False

def write_missing_queries(missing_qids, orig_fasta, out_fasta):
    """
    Filters a FASTA file, writing out only sequences with IDs in a given set.
    """
    if not missing_qids:
        return False
    missing_set = set(missing_qids)
    try:
        with open(orig_fasta, "r") as in_handle, open(out_fasta, "w") as out_handle:
            missing_records = (
                record for record in SeqIO.parse(in_handle, "fasta")
                if record.id in missing_set
            )
            count = SeqIO.write(missing_records, out_handle, "fasta")
        if count > 0:
            if log:
                log.write("Wrote {} missing queries to {}\n".format(count, out_fasta))
            return True
        else:
            if log:
                log.write(
                    "WARN: Could not find any records for {} missing IDs in {}\n"
                    .format(len(missing_qids), orig_fasta)
                )
            return False
    except (IOError, ValueError) as e:
        sys.stderr.write("ERROR writing missing queries to {}: {}\n".format(out_fasta, e))
        if log:
            log.write(traceback.format_exc())
        return False

def merge_intervals(reads):
    """
    Merges overlapping or adjacent integer intervals.
    Args:
    reads (list of lists): A list of [start, end] pairs representing alignment blocks.
    Returns:
    list of lists: A list of merged [start, end] intervals.
    """
    if not reads:
        return []
    valid_reads = []
    for i, read in enumerate(reads):
        try:
            if len(read) == 2:
                valid_reads.append(sorted([int(read[0]), int(read[1])]))
            else:
                if log:
                    log.write("WARN: Skipping invalid alignment format [{}]: {}\n".format(i, read))
        except (ValueError, TypeError):
            if log:
                log.write("WARN: Skipping non-integer coordinate [{}]: {}\n".format(i, read))
    if not valid_reads:
        return []
    valid_reads.sort(key=lambda interval: interval[0])
    merged = []
    current_start, current_end = valid_reads[0]
    for next_start, next_end in valid_reads[1:]:
        if next_start <= current_end + 1:
            current_end = max(current_end, next_end)
        else:
            merged.append([current_start, current_end])
            current_start, current_end = next_start, next_end
    merged.append([current_start, current_end])
    return merged

def one_block(contigs, q_len, end_dist_thresh, max_cov_thresh):
    """
    Analyzes a single alignment block to find a terminal element or classify it.
    Returns a tuple: (status, (start, end, length), details_string).
    - status: 'TERMINAL_ELEMENT', 'FULL_COVERAGE', or 'INTERNAL_BLOCK'.
    """
    if not contigs or not contigs[0]:
        return 'NO_CONTIGS', (0, 0, 0), "N/A"
    block_start, block_end = contigs[0]
    block_len = block_end - block_start + 1
    block_coverage = (block_len / float(q_len)) * 100
    if block_coverage > max_cov_thresh:
        details = "{:.2f}%".format(block_coverage)
        return 'FULL_COVERAGE', (-1, -1, -1), details
    start_dist = block_start - 1
    end_dist = q_len - block_end
    is_close_to_start = (start_dist <= end_dist_thresh)
    is_close_to_end = (end_dist <= end_dist_thresh)
    if is_close_to_start and not is_close_to_end:
        estart, eend, elen = block_end + 1, q_len, q_len - block_end
        details = "{}bp from 5' end".format(start_dist)
        return 'TERMINAL_ELEMENT', (estart, eend, elen), details
    elif not is_close_to_start and is_close_to_end:
        estart, eend, elen = 1, block_start - 1, block_start - 1
        details = "{}bp from 3' end".format(end_dist)
        return 'TERMINAL_ELEMENT', (estart, eend, elen), details
    else:
        details = "{}bp from ends".format(start_dist)
        return 'INTERNAL_BLOCK', (0, 0, 0), details

def unzip_pairs(joined):
    """
    Splits a list of pairs into two separate lists.
    """
    if not joined:
        return [], []
    list1, list2 = zip(*joined)
    return list(list1), list(list2)

def find_largest_gap_between_blocks(starts, ends):
    """
    Finds the largest gap between consecutive sorted intervals.
    """
    if len(starts) != len(ends) or len(starts) < 2:
        return 0, 0, 0
    largest_gap = -1
    gap_start, gap_end = 0, 0
    for i in range(len(ends) - 1):
        try:
            current_end, next_start = int(ends[i]), int(starts[i+1])
            gap = next_start - current_end - 1
            if gap > largest_gap:
                largest_gap, gap_start, gap_end = gap, current_end + 1, next_start - 1
        except (ValueError, TypeError):
            continue
    return (gap_start, gap_end, largest_gap) if largest_gap >= 0 else (0, 0, 0)

def extract_seq(seq_content, qid, estart, eend):
    """
    Extracts a subsequence from a string containing FASTA-formatted data.
    Note: Coordinates are 1-based and inclusive. Use eend=-1 for "to the end".
    """
    header = ">{}".format(qid)
    header_pos = seq_content.find(header)
    if header_pos == -1:
        if log:
            log.write("WARN: Header '{}' not found for extraction.\n".format(header))
        return ""
    seq_block_start = seq_content.find('\n', header_pos) + 1
    if seq_block_start == 0:
        return ""  # Header found but no sequence data follows
    next_header_pos = seq_content.find('>', seq_block_start)
    if next_header_pos == -1:
        seq_block = seq_content[seq_block_start:]
    else:
        seq_block = seq_content[seq_block_start:next_header_pos]
    # Consolidate sequence by removing line breaks
    full_sequence = seq_block.replace('\n', '').replace('\r', '')
    try:
        # Convert 1-based inclusive coordinates to 0-based slice indices
        start_index = int(estart) - 1
        end_index = int(eend) if eend != -1 else None
        if start_index < 0:
            start_index = 0
        return full_sequence[start_index:end_index]
    except (ValueError, TypeError):
        sys.stderr.write("ERROR: Invalid coordinates for extraction: {}, {}\n".format(estart, eend))
        return ""

def write_element_files(qid, qseq_content, elem_start, elem_end, output_dir, color_rgb, log_fh):
    """
    Writes FASTA and GenBank feature files for a discovered element.
    """
    try:
        qid_outdir = os.path.join(output_dir, str(qid))
        if not os.path.exists(qid_outdir):
            os.makedirs(qid_outdir)

        element_sequence = extract_seq(qseq_content, qid, elem_start, elem_end)
        if not element_sequence:
            log_fh.write("[{}] WARN: Extracted sequence is empty. Skipping file output.\n".format(qid))
            return

        # --- Write FASTA file ---
        fasta_fpath = os.path.join(qid_outdir, "{}_element.fasta".format(qid))
        try:
            with open(fasta_fpath, 'w') as fh:
                fh.write(">{}_element_{}-{}\n".format(qid, elem_start, elem_end))
                # Use textwrap for clean, standard line wrapping
                import textwrap
                fh.write("\n".join(textwrap.wrap(element_sequence, 60)) + "\n")
            log_fh.write("[{}] Wrote element FASTA: {}\n".format(qid, os.path.basename(fasta_fpath)))
        except IOError as e:
            err_msg = "ERROR writing element FASTA for {}: {}\n".format(qid, e)
            sys.stderr.write(err_msg)
            log_fh.write("[{}] {}\n".format(qid, err_msg))

        # --- Write GenBank feature table ---
        gb_fpath = os.path.join(qid_outdir, "{}_element.gb".format(qid))
        try:
            full_sequence = extract_seq(qseq_content, qid, 1, -1)  # -1 to get full sequence
            with open(gb_fpath, 'w') as fh:
                fh.write("FEATURES             Location/Qualifiers\n")
                fh.write("     misc_feature    {}..{}\n".format(elem_start, elem_end))
                fh.write("                     /label=\"element\"\n")
                fh.write("                     /color=\"{} {} {}\"\n".format(*color_rgb))
                fh.write("ORIGIN\n")
                for i in range(0, len(full_sequence), 60):
                    line_prefix = "{:>9} ".format(i + 1)
                    line_seq = full_sequence[i:i+60]
                    # Group into blocks of 10
                    chunks = ["".join(line_seq[j:j+10]) for j in range(0, len(line_seq), 10)]
                    fh.write("{}{}\n".format(line_prefix, " ".join(chunks)))
                fh.write("//\n")
            log_fh.write("[{}] Wrote feature table: {}\n".format(qid, os.path.basename(gb_fpath)))
        except IOError as e:
            err_msg = "ERROR writing feature table for {}: {}\n".format(qid, e)
            sys.stderr.write(err_msg)
            log_fh.write("[{}] {}\n".format(qid, err_msg))
    except Exception as e:
        err_msg = "ERROR during output generation for query {}: {}\n".format(qid, e)
        sys.stderr.write(err_msg)
        if log_fh:
            log_fh.write(err_msg + traceback.format_exc())

def rejection_reason(reasons_list, params):
    """
    Determines the single most significant rejection reason for a query from a list of all its hit rejections.
    """
    if not reasons_list:
        return "REJECTED: No valid element found in any hit", "N/A"

    priority_order = {
        "REJECTED: Single block covers too much of query (>{maxcov}%)".format(maxcov=params['maxcov']): 1,
        "REJECTED: Single block too far from query ends (>{enddist}bp)".format(enddist=params['enddist']): 2,
        "REJECTED: No significant gap found between blocks": 3,
        "REJECTED: Element length not in range ({minlen}-{maxlen}bp)".format(
            minlen=params['minlen'], maxlen=params['maxlen']
        ): 4,
    }

    reasons_list.sort(key=lambda x: (priority_order.get(x[0], 99), x[1]), reverse=True)
    best_reason_key, _, best_details_str = reasons_list[0]
    all_details_for_reason = [r[2] for r in reasons_list if r[0] == best_reason_key]

    # Sort details numerically where possible before creating the summary string
    try:
        sorted_details = sorted(
            list(set(all_details_for_reason)),
            key=lambda x: -float(x.split('(')[-1].split('%')[0].split('bp')[0])
        )
    except (ValueError, IndexError):
        sorted_details = sorted(list(set(all_details_for_reason)))

    final_details = ", ".join(sorted_details[:10])
    if len(all_details_for_reason) > 10:
        final_details += "..."

    return best_reason_key, final_details

def main(params, qseq_content, log, blast_log_fh):
    """
    Main analysis pipeline for finding insertion elements.
    """
    # 1. RUN BLAST (if not provided)
    final_blast_tab = os.path.join(params["out"], "blastn.tab")
    all_fasta_qids = params['qid_original_from_fasta']

    if params.get('tab') is None:
        log.write('\nStarting BLASTn search...\n')
        log.flush()
        init_blast_runtime, blast_tab_path = run_blast(params, params['q'], out_fname="blastn.tab")
        print('BLASTn finished in: {}'.format(init_blast_runtime))
        log.write('BLASTn finished in: {}\n'.format(init_blast_runtime))
        log.flush()
        if not blast_tab_path or not os.path.exists(blast_tab_path):
            sys.exit("CRITICAL: BLAST failed to produce an output table. See blastn.log. Exiting.")
    else:
        log.write("\nUsing pre-computed BLAST table: {}\n".format(params['tab']))
        log.flush()
        shutil.copy(params['tab'], final_blast_tab)

    # 2. RUN SUPPLEMENTARY BLAST FOR MISSING QUERIES
    rerun_count = 0
    MAX_COMPLEMENTARY_RUNS = 5
    while rerun_count < MAX_COMPLEMENTARY_RUNS:
        missing_qids = get_missing_queries(final_blast_tab, all_fasta_qids)
        if not missing_qids:
            log.write("\nBLAST results complete for all {} queries.\n".format(len(all_fasta_qids)))
            log.flush()
            break
        rerun_count += 1
        log.write(
            "\nCheck #{}: {} queries missing. Running complementary BLAST.\n"
            .format(rerun_count, len(missing_qids))
        )
        log.write(
            "Missing: {}{}\n".format(', '.join(missing_qids[:10]), '...' if len(missing_qids) > 10 else '')
        )
        log.flush()
        tmp_missing_fasta = os.path.join(params["out"], "missing_queries_temp.fasta")
        if not write_missing_queries(missing_qids, params['q'], tmp_missing_fasta):
            log.write("ERROR: Failed to write missing queries. Halting checks.\n")
            log.flush()
            break
        supp_params = params.copy()
        if supp_params.get('run') != 'local':
            supp_params['run'] = 'web'
        log.write("Supplementary BLAST mode: '{}'.\n".format(supp_params['run']))
        log.flush()
        tmp_supp_blast_tab = os.path.join(params["out"], "blastn_extra_temp.tab")
        supp_blast_runtime, supp_tab_path = run_blast(
            supp_params, tmp_missing_fasta, out_fname=os.path.basename(tmp_supp_blast_tab)
        )
        log.write('Supplementary BLAST #{} finished in: {}\n'.format(rerun_count, supp_blast_runtime))
        log.flush()
        if supp_tab_path and os.path.exists(supp_tab_path):
            log.write("Appending new results to {}.\n".format(os.path.basename(final_blast_tab)))
            log.flush()
            with open(final_blast_tab, 'ab') as master_fh, open(supp_tab_path, 'rb') as extra_fh:
                for line in extra_fh:
                    if not line.startswith(b'#'):
                        master_fh.write(line)
        else:
            log.write("WARN: Supplementary BLAST produced no results.\n")
            log.flush()
        if os.path.exists(tmp_missing_fasta):
            os.remove(tmp_missing_fasta)
        if os.path.exists(tmp_supp_blast_tab):
            os.remove(tmp_supp_blast_tab)
    else:
        final_missing = get_missing_queries(final_blast_tab, all_fasta_qids)
        if final_missing:
            log.write(
                "\nWARN: After {} retries, {} queries still missing.\n"
                .format(MAX_COMPLEMENTARY_RUNS, len(final_missing))
            )
            log.flush()

    if not os.path.exists(final_blast_tab):
        sys.exit("CRITICAL: Final BLAST table '{}' not found. Exiting.".format(final_blast_tab))

    # 3. PARSE AND PROCESS RESULTS
    params['tab'] = final_blast_tab
    _, blast_cols, blast_hits = blast_parse(params['tab'])
    log.write("\nProcessing final BLAST table: {}\n".format(os.path.basename(params['tab'])))
    log.flush()
    print("Processing final BLAST table: {}...".format(os.path.basename(params['tab'])))

    query_count, element_count = 0, 0
    analysis_summary = {
        "VALID: Element found": [],
        "REJECTED: No hits in final table.": [],
        "REJECTED: All hits filtered out by coverage.": [],
        "REJECTED: Single block covers too much of query (>{maxcov}%)".format(maxcov=params['maxcov']): [],
        "REJECTED: Single block too far from query ends (>{enddist}bp)".format(enddist=params['enddist']): [],
        "REJECTED: No significant gap found between blocks": [],
        "REJECTED: Element length not in range ({minlen}-{maxlen}bp)".format(
            minlen=params['minlen'], maxlen=params['maxlen']
        ): [],
        "REJECTED: No valid element found in any hit": [],
    }

    q_id_col = 'query id'
    s_id_col = 'subject id'
    qcov_col = '% query coverage per subject'
    qlen_col = 'query length'
    qstart_col = 'q. start'
    qend_col = 'q. end'

    summary_fpath = os.path.join(params["out"], 'elements.txt')
    coordinates_fpath = os.path.join(params["out"], 'elements.csv')

    with open(summary_fpath, 'w') as summ_fh:
        summ_fh.write('smith_sa v{}\n'.format(VERSION))
        summ_fh.write("\nQuery file: {}".format(params["q"]))
        if args.tab:
            summ_fh.write("\nOriginal BLAST table: {}".format(args.tab))
        summ_fh.write("\nFinal BLAST table: {}".format(params["tab"]))
        db_info = params.get("d", "nt") if params.get('run') == 'local' else 'nt'
        summ_fh.write("\nDatabase ({}): {}".format(params.get('run', 'N/A'), db_info))
        if 'org' in params:
            summ_fh.write("\nTaxids: {}".format(args.org or "N/A"))
        summ_fh.write("\nElement length (bp): {}-{}".format(params["minlen"], params["maxlen"]))
        summ_fh.write("\nQuery coverage (%): {}-{}".format(params["mincov"], params["maxcov"]))
        summ_fh.write("\nMax end distance (bp): {}".format(params["enddist"]))
        if 'cpu' in params and params.get('run') == 'local':
            summ_fh.write("\nCPU threads: {}".format(params["cpu"]))
        summ_fh.write(
            "\n\nQuery ID\tElement Subject ID\tSubject Avg. % Query Cov\tNo. of Contigs in Subject\t"
            "Element identification\tElement 5' coordinate\tElement 3' coordinate\tElement length\tValid\n"
        )

        blast_hits_data = []
        if blast_cols and blast_hits:
            for row in blast_hits:
                hit_dict = {
                    str(col).strip(): str(row[i]).strip() if i < len(row) else ""
                    for i, col in enumerate(blast_cols)
                }
                for key in [qcov_col, qlen_col, qstart_col, qend_col]:
                    try:
                        hit_dict[key] = float(hit_dict.get(key, 0.0))
                    except (ValueError, TypeError):
                        hit_dict[key] = 0.0
                blast_hits_data.append(hit_dict)

        for qid in sorted(list(all_fasta_qids)):
            query_count += 1
            element_found_for_qid = False
            log.write("\n----------\n[{}] Processing Query\n".format(qid))
            hits_for_qid = [h for h in blast_hits_data if h.get(q_id_col) == qid]
            if not hits_for_qid:
                log.write("[{}] No hits found.\n".format(qid))
                summ_fh.write(
                    "{}\tno hits in final table\t\t\tno\t\t\t\tno\n".format(qid)
                )
                analysis_summary["REJECTED: No hits in final table."].append(qid)
                continue

            log.write("[{}] Found {} total hits.\n".format(qid, len(hits_for_qid)))
            q_len = int(hits_for_qid[0].get(qlen_col, 0))
            log.write("[{}] Query length: {} bp\n".format(qid, q_len))
            sids_in_query = {h[s_id_col] for h in hits_for_qid if s_id_col in h}
            candidates = []
            all_hits_outside_coverage = True

            for sid in sids_in_query:
                hits_for_sid = [h for h in hits_for_qid if h.get(s_id_col) == sid]
                covs = [h[qcov_col] for h in hits_for_sid if isinstance(h.get(qcov_col), float)]
                avg_cov = sum(covs) / len(covs) if covs else 0
                if params['mincov'] <= avg_cov:
                    all_hits_outside_coverage = False
                # Only add to candidates if it's potentially valid, not just above mincov
                if params['mincov'] <= avg_cov <= params['maxcov']:
                    reads = [[h[qstart_col], h[qend_col]] for h in hits_for_sid]
                    contigs = merge_intervals(reads)
                    candidates.append({
                        'subject_id': sid,
                        'avg_cov': avg_cov,
                        'num_contigs': len(contigs),
                        'contigs': contigs,
                    })

            if not candidates and all_hits_outside_coverage:
                log.write("[{}] All hits filtered out by coverage.\n".format(qid))
                summ_fh.write(
                    "{}\tno valid hits\t\t\tno\t\t\t\tno\n".format(qid)
                )
                analysis_summary["REJECTED: All hits filtered out by coverage."].append(
                    "{}: N/A".format(qid)
                )
                continue

            candidates.sort(key=lambda x: x['avg_cov'], reverse=True)
            rejection_reasons_for_query = []

            for cand in candidates:
                curr_sid = cand['subject_id']
                curr_avg_cov = cand['avg_cov']
                num_contigs = cand['num_contigs']
                contigs = cand['contigs']
                is_valid = False

                if num_contigs == 1:
                    status, coords, details = one_block(
                        contigs, q_len, params['enddist'], params['maxcov']
                    )
                    elem_start, elem_end, elem_len = coords
                    log.write("[{}] Analyzing Single Block hit ({}): Coords={}-{}, Len={}, Cov={:.2f}%, Status={}\n"
                              .format(qid, curr_sid, elem_start, elem_end, elem_len, curr_avg_cov, status))

                    if status == 'TERMINAL_ELEMENT':
                        if params['minlen'] <= elem_len <= params['maxlen']:
                            is_valid = True
                        else:
                            rejection_reasons_for_query.append((
                                "REJECTED: Element length not in range ({minlen}-{maxlen}bp)".format(
                                    minlen=params['minlen'], maxlen=params['maxlen']
                                ),
                                elem_len,
                                "{} ({}..{}, {}bp, {:.2f}%)".format(
                                    curr_sid, elem_start, elem_end, elem_len, curr_avg_cov
                                ),
                            ))
                    elif status == 'FULL_COVERAGE':
                        rejection_reasons_for_query.append((
                            "REJECTED: Single block covers too much of query (>{maxcov}%)".format(
                                maxcov=params['maxcov']
                            ),
                            curr_avg_cov,
                            "{} (-1bp, {})".format(curr_sid, details),
                        ))
                    elif status == 'INTERNAL_BLOCK':
                        rejection_reasons_for_query.append((
                            "REJECTED: Single block too far from query ends (>{enddist}bp)".format(
                                enddist=params['enddist']
                            ),
                            0,
                            "{} ({}, {:.2f}%)".format(curr_sid, details, curr_avg_cov),
                        ))

                elif num_contigs > 1:
                    starts, ends = unzip_pairs(contigs)
                    elem_start, elem_end, elem_len = find_largest_gap_between_blocks(starts, ends)
                    
                    log.write("[{}] Analyzing Multi-Block hit ({}): GapStart={}, GapEnd={}, GapLen={}, Cov={:.2f}%\n"
                              .format(qid, curr_sid, elem_start, elem_end, elem_len, curr_avg_cov))

                    if elem_len > 0:
                        if params['minlen'] <= elem_len <= params['maxlen']:
                            is_valid = True
                        else:
                            rejection_reasons_for_query.append((
                                "REJECTED: Element length not in range ({minlen}-{maxlen}bp)".format(
                                    minlen=params['minlen'], maxlen=params['maxlen']
                                ),
                                elem_len,
                                "{} ({}..{}, {}bp, {:.2f}%)".format(
                                    curr_sid, elem_start, elem_end, elem_len, curr_avg_cov
                                ),
                            ))
                    else:
                        rejection_reasons_for_query.append((
                            "REJECTED: No significant gap found between blocks",
                            0,
                            "{} ({:.2f}%)".format(curr_sid, curr_avg_cov),
                        ))

                if is_valid:
                    log.write(
                        "[{}] Candidate found in {}: {}..{} ({} bp, {:.2f}%).\n"
                        .format(qid, curr_sid, elem_start, elem_end, elem_len, curr_avg_cov)
                    )
                    summ_fh.write(
                        "{}\t{}\t{:.2f}\t{}\tyes\t{}\t{}\t{}\tyes\n".format(
                            qid, curr_sid, curr_avg_cov, num_contigs,
                            elem_start, elem_end, elem_len
                        )
                    )
                    with open(coordinates_fpath, 'a') as coord_fh:
                        coord_fh.write("{};{};{}\n".format(qid, elem_start, elem_end))
                    write_element_files(
                        qid, qseq_content, elem_start, elem_end,
                        params["out"], params['color'], log
                    )
                    analysis_summary["VALID: Element found"].append(
                        "{}: {} ({}..{}, {}bp, {:.2f}%)".format(
                            qid, curr_sid, elem_start, elem_end, elem_len, curr_avg_cov
                        )
                    )
                    element_found_for_qid = True
                    break

            if element_found_for_qid:
                element_count += 1
            else:
                primary_reason, details = rejection_reason(rejection_reasons_for_query, params)
                log.write(
                    "[{}] No valid element found. Primary rejection reason: {}\n"
                    .format(qid, primary_reason)
                )
                analysis_summary[primary_reason].append(
                    "{} | {}".format(qid, details)
                )
                summ_fh.write(
                    "{}\tno valid element identified\t\t\tno\t\t\t\tno\n".format(qid)
                )

        summ_fh.write('\n\n----- Summary -----\n')
        summ_fh.write('Processed {} queries.\n'.format(query_count))
        summ_fh.write('Found {} valid elements.\n'.format(element_count))
        summ_fh.write('Total execution time: {}\n'.format(datetime.now() - start_time))

    return query_count, element_count, analysis_summary

if __name__ == "__main__":
    log = None
    blast_log_fh = None
    query_count, element_count = 0, 0
    analysis_summary = {}

    try:
        if not len(sys.argv) > 1 or args.help:
            print(HELP_MSG)
            sys.exit(0)
        elif args.version:
            print(VERSION)
            sys.exit(0)

        is_valid, params = validate_args(args)
        if not is_valid:
            sys.stderr.write("ERROR: Invalid arguments. See help (-h) for details.\n")
            sys.exit(1)

        log_fpath = os.path.join(params["out"], "run.log")
        log = open(log_fpath, 'w')

        blast_log_fpath = os.path.join(params["out"], "blastn.log")
        blast_log_fh = open(blast_log_fpath, 'a')
        blast_log_fh.write("\n=== New Run at {} ===\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        blast_log_fh.flush()

        log.write("smith_sa v{}\n".format(VERSION))
        log.write("Last modified datetime of this script: {}\n".format(str(last_modified_datetime)))
        log.write("(c) 2021. Arthur Gruber & Giuliana Pola\n")
        log.write("https://github.com/GiulianaPola/smith_sa\n")
        log.write("Run started at {}\n".format(start_time.strftime("%Y-%m-%d, %H:%M:%S")))
        user = os.getlogin() if hasattr(os, "getlogin") else os.environ.get("USER", "N/A")
        log.write("User: {}\n".format(user))
        log.write("Workspace: {}\n".format(cwd_path))
        log.write("Output dir: {}\n".format(params["out"]))
        log.write("Command line: {}\n".format(" ".join(sys.argv)))
        for pname, pval in sorted(params.items()):
            if pname not in ["qid_original_from_fasta", "current_blastn_parts_dir"]:
                log.write("{}: {}\n".format(pname, pval))
        log.flush()

        with open(params["q"], "r") as qfh:
            qseq_content = qfh.read()
        all_fasta_qids = get_fasta_query_ids(params["q"])
        if not all_fasta_qids:
            raise RuntimeError("No query IDs found in {}. Exiting.".format(params["q"]))
        log.write(
            "Found {} unique queries in {}: {}\n".format(
                len(all_fasta_qids), os.path.basename(params["q"]), ",".join(list(all_fasta_qids))
            )
        )
        log.flush()
        params["qid_original_from_fasta"] = all_fasta_qids

        query_count, element_count, analysis_summary = main(params, qseq_content, log, blast_log_fh)

    except Exception as err:
        traceback_str = traceback.format_exc()
        sys.stderr.write("ERROR: Unhandled exception in main block.\n" + traceback_str)
        if log:
            log.write("ERROR: Unhandled exception.\n" + traceback_str)
            log.flush()
    finally:
        if log and hasattr(log, "write") and not log.closed:
            log.write("\n----- Run Summary -----\n")
            log.write("Processed {} queries.\n".format(query_count))
            log.write("Found {} valid elements.\n".format(element_count))
            log.write("----- Analysis Breakdown by Rejection Reason -----\n")
            reason_order = [
                "VALID: Element found",
                "REJECTED: No hits in final table.",
                "REJECTED: All hits filtered out by coverage.",
                "REJECTED: Single block covers too much of query (>{maxcov}%)".format(
                    maxcov=params.get("maxcov", 90)
                ),
                "REJECTED: Single block too far from query ends (>{enddist}bp)".format(
                    enddist=params.get("enddist", 50)
                ),
                "REJECTED: No significant gap found between blocks",
                "REJECTED: Element length not in range ({minlen}-{maxlen}bp)".format(
                    minlen=params.get("minlen", 4000),
                    maxlen=params.get("maxlen", 50000),
                ),
                "REJECTED: No valid element found in any hit",
            ]
            for reason in reason_order:
                valuelist = analysis_summary.get(reason, [])
                count = len(valuelist)
                if count == 0:
                    continue
                if reason.startswith("VALID"):
                    log.write("- {} {}\n".format(count, reason.replace("VALID: ", "")))
                    log.write("  {}\n".format(", ".join(valuelist)))
                else:
                    displaylist = [v.split("|")[0].strip() for v in valuelist]
                    log.write("- {} {}\n".format(count, reason.replace("REJECTED: ", "")))
                    log.write("  {}\n".format(", ".join(displaylist)))
            log.write("Run finished at: {}\n".format(datetime.now().strftime("%Y-%m-%d, %H:%M:%S")))
            log.write("Total execution time: {}\n".format(datetime.now() - start_time))
            log.write("\nEnd of execution")
            log.flush()
            log.close()

        if blast_log_fh and not blast_log_fh.closed:
            blast_log_fh.write("=== End of Run at {} ===\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            blast_log_fh.flush()
            blast_log_fh.close()

        print("\nEnd of execution")