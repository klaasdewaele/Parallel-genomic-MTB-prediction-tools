#! /data/leuven/348/vsc34807/miniconda3/bin/python
import subprocess
import sys
import os
from os import walk
from datetime import datetime
from time import strftime
import json
import pandas as pd
import pprint
import argparse
from os.path import expanduser
import timeit
from collections import Counter
import fileinput
import shutil
import re
from IPython.display import display, HTML
import numpy as np
from mailmerge import MailMerge
from who_treatment import DstProfile, regimens

drugs = {'isoniazid': 'INH', 'rifampicin': 'RIF', 'streptomycin': 'STR', 'ethambutol': 'EMB',
              'pyrazinamide': 'PYR', 'fluoroquinolones': 'FLQ', 'ofloxacin': 'OFL', 'ciprofloxacin': 'CIP',
              'levofloxacin': 'LEV', 'moxifloxacin': 'MOX', 'gatifloxacin': 'GAT', 'aminoglycosides': 'AGL',
              'amikacin': 'AMK', 'kanamycin': 'KAN', 'para-aminosalicylic_acid': 'PAS', 'capreomycin': 'CAP',
              'ethionamide': 'ETH', 'cyclosporine': 'CYC', 'linezolid': 'LIN', 'prothionamide': 'PRO',
              'rifabutin': 'RIB'}

drug_dict = {'isoniazid': 'INH', 'rifampicin': 'RIF', 'streptomycin': 'STR', 'ethambutol': 'EMB',
             'pyrazinamide': 'PYR', 'fluoroquinolones': 'FLQ', 'ofloxacin': 'OFL', 'ciprofloxacin': 'CIP',
             'levofloxacin': 'LEV', 'moxifloxacin': 'MOX', 'gatifloxacin': 'GAT', 'aminoglycosides': 'AGL',
             'amikacin': 'AMK', 'kanamycin': 'KAN', 'para-aminosalicylic_acid': 'PAS', 'capreomycin': 'CAP',
             'ethionamide': 'ETH', 'cyclosporine': 'CYC', 'linezolid': 'LIN', 'prothionamide': 'PRO',
             'rifabutin': 'RIB', 'inh': 'INH', 'rif': 'RIF', 'pza': 'PYR', 'emb': 'EMB', 'str': 'STR', 'eth': 'ETH',
             'kan': 'KAN', 'cap': 'CAP', 'amk': 'AMK', 'cip': 'CIP', 'levo': 'LEV', 'oflx': 'OFL', 'pas': 'PAS',
             'rib': 'RIB', 'mox': 'MOX', 'gat': 'GAT', 'cyc': 'CYC', 'lin': 'LIN', 'LZD': 'LIN', 'pro': 'PRO',
             'FQ': 'FLQ', 'RMP': 'RIF', 'SM': 'STR', 'CPR': 'CAP', 'PZA': 'PYR', 'INH': 'INH', 'RIF': 'RIF',
             'EMB': 'EMB', 'ETH': 'ETH', 'OFL': 'OFL', 'CIP': 'CIP', 'LEV': 'LEV', 'MOX': 'MOX',
             'GAT': 'GAT', 'AMK': 'AMK', 'AMI': 'AMK', 'KAN': 'KAN', 'PAS': 'PAS', 'CAP': 'CAP',
             'BDQ': 'BED', 'AGL': 'AGL'}

# Possible associations between genes and drugs: Table 6 from https://apps.who.int/iris/rest/bitstreams/1352664/retrieve
gene_dict = {'rpoB': 'RIF', 'rpoC': 'RIF', 'Rv2752c': 'RIF, INH', 'rpoA': 'RIF', 'inhA': 'INH, ETH', 'katG': 'INH',
             'ndh': 'INH',
             'fabG1': 'INH, ETH', 'Rv1258c': 'INH, PYR, STR', 'mshA': 'INH', 'ahpC': 'INH', 'ndh': 'INH', 'embB': 'EMB',
             'embA': 'EMB', 'embC': 'EMB', 'ubiA': 'EMB', 'embR': 'EMB', 'pncA': 'PYR', 'PPE35': 'PYR',
             'Rv3236c': 'PYR', 'clpC1': 'PYR', 'panD': 'PYR', 'gyrA': 'FLQ', 'gyrB': 'FLQ', 'mmpL5': 'BED, CLO',
             'Rv1979': 'BED, CLO', 'rplC': 'LIN', 'rrs': 'AMK, STR, LIN, KAN, CAP', 'rrl': 'LIN', 'ddn': 'DEL',
             'eis': 'AMK, KAN',
             'aftB': 'AMK, CAP', 'ccsA': 'AMK', 'whiB6': 'AMK, STR', 'fprA': 'AMK, CAP', 'whiB7': 'AMK, STR, KAN',
             'rpsL': 'STR',
             'gid': 'STR', 'ethA': 'ETH', 'pepQ': 'BED', 'Rv0678': 'BED, CLO', 'mmpS5': 'BED, CLO', 'atpE': 'BED',
             'fgd1': 'DEL', 'fbiA': 'DEL', 'fbiB': 'DEL', 'fbiC': 'DEL', 'Rv2983': 'DEL', 'ethR': 'ETH',
             'Rv3083': 'ETH', 'ndh': 'ETH', 'tlyA': 'CAP', 'ccsA': 'CAP', 'BDQ': 'BED'}

# Similar to gene_dict with drugs as keys - used to determine which genes are shown for each drug
drug_gene_dict = {'INH': ['inhA', 'katG', 'ahpC', 'mshA', 'ndh', 'Rv1258c', 'Rv2752c'],
                  'RIF': ['rpoB', 'rpoA', 'rpoC', 'Rv2752c'],
                  'EMB': ['embA', 'embB', 'embC', 'embR', 'ubiA'],
                  'PYR': ['pncA', 'clpC1', 'panD', 'Rv1258c', 'PPE35', 'Rv3236c'],
                  'FLQ': ['gyrA', 'gyrB'],
                  'CIP': ['gyrA', 'gyrB'],
                  'MOX': ['gyrA', 'gyrB'],
                  'LEV': ['gyrA', 'gyrB'],
                  'GAT': ['gyrA', 'gyrB'],
                  'OFL': ['gyrA', 'gyrB'],
                  'BED': ['pepQ', 'Rv0678', 'mmpL5', 'mmpS5', 'atpE', 'Rv1979c'],
                  'DEL': ['fgd1', 'ddn', 'fbiA', 'fbiB', 'fbiC', 'Rv2983'],
                  'AMK': ['rrs', 'eis', 'whiB7', 'whiB6', 'ccsA', 'fprA', 'aftB'],
                  'STR': ['rrs', 'rpsL', 'gid', 'whiB7', 'Rv1258c', 'whiB6'],
                  'ETH': ['inhA', 'ethA', 'ethR', 'mshA', 'Rv3083', 'ndh'],
                  'KAN': ['rrs', 'eis', 'whiB7'],
                  'CAP': ['rrs', 'tlyA', 'whiB6', 'ccsA', 'fprA', 'aftB'],
                  'PAS': ['thyA', 'dfrA', 'folC', 'folP1', 'folP2', 'thyX'],
                  # https://www.nature.com/articles/s41598-019-48940-5
                  'AGL': ['rrs', 'eis', 'whiB7', 'whiB6', 'ccsA', 'fprA', 'aftB', 'Rv1258c', 'gid', 'rpsL'],
                  # composite from AMK and STR
                  'AMK, STR, LIN, KAN, CAP': ['rrs'],
                  'LIN': ['rplC', 'rrs', 'rrl'],
                  '?': list(gene_dict.keys())
                  }

# use mykrobe_reg_dict to convert drug abbreviations to capitalized full drug names used in Hunt et al. 2019 who_regimen.py
mykrobe_reg_dict = {
    'INH': 'Isoniazid',
    'RIF': 'Rifampicin',
    'PYR': 'Pyrazinamide',
    'EMB': 'Ethambutol',
    'KAN': 'Kanamycin',
    'AMK': 'Amikacin',
    'CAP': 'Capreomycin',
    'STR': 'Streptomycin',
    'OFL': 'Ofloxacin',
    'CIP': 'Ciprofloxacin',
    'MOX': 'Moxifloxacin',
    'BED': 'Bedaquiline',
    'LIN': 'Linezolid',
    'CLO': 'Clofazimide',
    'CYC': 'Cycloserine',
    'TER': 'Terizidone'
}

def find_files(path):
    all_files = os.listdir(path) # List all files in path
    fastq_files = [] # A list that will contain fastq-files 
    id_no_direction_list = []
    id_no_direction_dict = {} 
    # Select fastq-files
    for file in all_files:
        ext = ('fastq.gz', 'fq.gz')
        if file.endswith(ext):
            # Append to fastq_files
            fastq_files.append(file)
            # Get filename without extension
            id_no_extension = re.split(r'(.fastq.gz|.fq.gz|.fastq)', file)[0]
            # Remove direction identifier ('1' or '2' preceded by non-numeric character; i.e. 'R1','_1', '-1')
            id_no_direction = re.split('\D(1|2)$', id_no_extension)[0]
            # Append to list: used to identify single-read files and paired-end files
            id_no_direction_list.append(id_no_direction)
            # Create dictionary that links filenames without direction to original file
            # If file id not yet added to dict: add
            if id_no_direction not in id_no_direction_dict:
                id_no_direction_dict[id_no_direction] = []
                id_no_direction_dict[id_no_direction].append(file)
            else:
                id_no_direction_dict[id_no_direction].append(file)

    print(f'\nFastq-files found: {fastq_files}. \nIdentifying single-end and paired-end read files ...')
    # Identify which files are paired-end vs single-end
    paired_end = []
    single_end = []
    # 'counts' is a dictionary that contains for each file in id_no_direction_list the number of occurrences
    counts = Counter(id_no_direction_list)
    for id, count in counts.items():
        if count == 2:
            paired_end.append(id)
        elif count == 1:
            single_end.append(id_no_direction_dict[id])
        else: print(f'Error for file {id}; could not determine whether file has single-reads or paired-end reads')
    
    print(f'Single-end read files: {single_end}')
    print(f'\nLooking for paired-end read files and formatting ...\n')
    # Rename paired-end files to comply to requirements by MTBseq
    # Problem: paired_end only contains filenames without extension and without direction identifier
    # Look up all paired_end files in id_no_direction_dict; rename; and return list of lists for paired-end files
    formatted = []
    for file in paired_end:
        # Fetch original file format including extension 
        original_files = id_no_direction_dict[file] # 'original_files' is a list of forward and reverse reads: ['sample_R1', 'sample_R2']
        # duo is a list that pairs forward and reverse reads
        duo = []
        # Loop through forward and reverse reads files:
        for reads in original_files:
            # Get file_identifier (original file name without extension)
            id_no_extension = re.split(r'(.fastq.gz|.fq.gz)', reads)[0]
            # Split file with '_' as delimiter
            id_split = id_no_extension.split('_') 
            # Format according to MTBseq-required format FILENAME_LIBRARYID_R1 or _R2.fastq.gz, see https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md
            R1_or_R2 = ('R1', 'R2')
            if ((len(id_split) > 2) and (id_split[-1] in R1_or_R2)):
                duo.append(reads)
                print(f'\tFile {id_no_extension} had correct format')
            else:
                lib_id = ''.join(id_split[1:-1]) # This captures parts of id_split between first element and last 
                new_name = f"{id_split[0]}_{lib_id if lib_id != '' else 'library-id'}_R{id_split[-1][-1]}.fastq.gz"
                os.rename(os.path.join(path, reads), os.path.join(path, new_name))
                print(f'\tFile {id_no_extension} has been renamed to {new_name}')
                duo.append(new_name)
        formatted.append(duo)
    
    print(f'\nPaired-end read files:')
    pprint.pprint(formatted)
    
    if len(fastq_files) != 0: # In case fastq-files have been found
        return formatted, single_end, (f'\nInput setup succeeded')
    else:
        return [], [], (f'\nNo valid files found in {path}. Exiting ...')

def start():
    # DISPLAY PIPELINE INFORMATION
    print('\n<<< LAG open-access MTB diagnosis pipeline >>>')
    CPUs = json.loads(subprocess.run(['grep', '-c', 'processor', '/proc/cpuinfo'], capture_output=True).stdout.strip().decode('utf-8'))
    RAM = subprocess.run(['free', '-h'], capture_output=True).stdout.strip().decode('utf-8').split()[7]
    print(f'\nRunning on {os.uname().sysname} {os.uname().version} @{os.uname().nodename}, no. of CPUs: {CPUs}, RAM:{RAM}')
    
    print('\nImplemented tools:')
    mykrobe_v = subprocess.run(['conda', 'run', '-n', 'mykrobe', 'mykrobe', '--v'], capture_output=True).stdout.strip().decode("utf-8").split()[1]
    tbprofiler_v = subprocess.run(['conda', 'run', '-n', 'tbprofiler', 'tb-profiler', 'version'], capture_output=True).stdout.strip().decode("utf-8").split()[2]
    gentb_v = 'v. 19  Aug 2021'
    mtbseq_v = subprocess.run(['conda', 'run', '-n', 'mtbseq', 'MTBseq', '--version'], capture_output=True).stdout.strip().decode()
    print(f'\tMykrobe version {mykrobe_v}')
    print(f'\tTBProfiler version {tbprofiler_v}')
    print(f'\tGenTB version {gentb_v}')
    print(f'\tMTBseq version {mtbseq_v}')
    return mykrobe_v, tbprofiler_v, gentb_v, mtbseq_v    

def input_setup(input_path):
    # GET VALID INPUT FILES
    if os.path.exists(input_path):
        print(f'\nInput directory: {input_path}') 
        input_dir = input_path
    else:
        print(f'\nNo valid input directory provided')
        input_dir = subprocess.run('pwd', capture_output=True).stdout.strip().decode('utf-8')
        print(f'Examining present directory: {input_dir}')
    formatted, single_end, comment = find_files(input_dir)
    print(comment)
    if 'Exiting' in comment:
        quit()
    
    return input_dir, formatted, single_end

def output_setup(output_path):
    # SET UP OUTPUT DIRECTORY
    session_time = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    if os.path.exists(output_path):
        session_folder = os.path.join(output_path, session_time)
    else:
        print('\nNo valid output directory provided')
        home = '/scratch/leuven/348/vsc34807/output'
        session_folder = os.path.join(home, session_time)
    
    # check whether session_folder exists: avoid duplicate session_folders when executing this script parallel at same time
    if os.path.exists(session_folder):
        seconds_plus_one = str(int(session_folder.split('-')[-1]) + 1)
        session_folder = '-'.join(session_folder.split('-')[:-1] + [seconds_plus_one])
    os.mkdir(session_folder)
    print(f'Created output folder at {session_folder}')

    return session_folder, session_time

def copy_rename_result(path, ext, output_dir, tool):
    # find files at source destination 'path'
    result_list = os.listdir(path)
    # Identify result file: this is a .json file (TBP, Mykrobe) or predict.json file (GenTB) or .tab (MTBseq); these extensions are in ext argument
    try: 
        result = [file for file in result_list if re.search(ext, file)][0] # Pick the first file from the list
        # For MTBseq we must retain info on mincovf, mincovr, minfreq, minphred20, that are in filename: 
        # extract this info from result (that has format: [SampleID]_[LibID]_[*].gatk_position_variants_[mincovf]_[mincovr]_[minfreq]_[minphred20]_[all_vars][snp_vars][lowfreq_vars].tab)
        if tool == 'mtbseq_susceptibility':
            info = result.split('_')
            info = '_'.join(info[(info.index('variants') + 1):(info.index('variants') + 5)] + [info[-1].split('.')[0]]) + '_' # get four parts after 'variants' part + info on all_vars, snp_vars, lowfreq_vars
        else: info = ''
        shutil.copy(os.path.join(path, result), os.path.join(output_dir, f'{tool}_{info}result'))
    except Exception as e:
        print(f'\tNo result file found for {tool}: {e}')

def mykrobe(input_dir, file, args, output_dir):
    file_name = str(file[0].split('.')[0])
    raw_output_folder = os.path.join(output_dir, 'mykrobe_raw')
    time = datetime.now().strftime('%d-%m-%Y, %H:%M:%S')
    os.mkdir(raw_output_folder)
    print(f"\tMykrobe analysis started at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}")
    log = subprocess.run(['conda', 'run', '-n', 'mykrobe', 'mykrobe', 'predict', '--sample', file_name, '--species', 'tb', '--output', f'{file_name}.json', '--format', 'json', '--seq', os.path.join(input_dir, file[0]), os.path.join(input_dir, file[1])], cwd=os.path.join(output_dir, 'mykrobe_raw'), capture_output=True)
    with open(os.path.join(output_dir, 'mykrobe_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))
    copy_rename_result(raw_output_folder, '.json$', output_dir, 'mykrobe') #copy_rename_result uses regex to find .json at end of string (regex module needed for finding correct mtbseq filenames) 
    print(f"\tMykrobe analysis finished at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}, elapsed time: {datetime.now() - datetime.strptime(time, '%d-%m-%Y, %H:%M:%S')}")
    if log.returncode == 0:
        stderr = 'STDERR empty'
    else:
        stderr = log.stderr.decode('utf-8')
    return ' '.join(log.args), str(log.returncode), stderr

def tbprofiler(input_dir, file, args, output_dir):
    file_name = str(file[0].split('.')[0])
    raw_output_folder = os.path.join(output_dir, 'tbprofiler_raw')
    os.mkdir(raw_output_folder)
    time = datetime.now().strftime('%d-%m-%Y, %H:%M:%S')
    print(f"\tTBProfiler analysis started at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}")
    log = subprocess.run(['conda', 'run', '-n', 'tbprofiler', 'tb-profiler', 'profile', '-1', str(os.path.join(input_dir, file[0])), '--dir', str(os.path.join(output_dir, 'tbprofiler_raw'))], capture_output=True)
    with open(os.path.join(output_dir, 'tbprofiler_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))
    copy_rename_result(os.path.join(raw_output_folder, 'results'), '.json$', output_dir, 'tbprofiler')
    print(f"\tTBProfiler analysis finished at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}, elapsed time: {datetime.now() - datetime.strptime(time, '%d-%m-%Y, %H:%M:%S')}")
    if log.returncode == 0:
        stderr = 'STDERR empty'
    else:
        stderr = log.stderr.decode('utf-8')
    return ' '.join(log.args), str(log.returncode), stderr

def mtbseq(input_dir, file, args, output_dir):
    os.mkdir(os.path.join(output_dir, 'mtbseq_raw'))
    time = datetime.now().strftime('%d-%m-%Y, %H:%M:%S')
    print(f"\tMTBseq analysis started at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}")
    # Copy files to mtbseq_raw: cannot execute MTBseq directly on input_dir since MTBseq will process ALL files at once
    for file_to_copy in file:
        input_path = os.path.join(input_dir, file_to_copy)
        subprocess.run(['cp', input_path, os.path.join(output_dir, 'mtbseq_raw', file_to_copy)])
        print(f'\t\tCopied {file_to_copy} to session folder')
    input_dir = os.path.join(output_dir, 'mtbseq_raw')
    log = subprocess.run(['conda', 'run', '-n', 'mtbseq', 'MTBseq', '--step', 'TBfull', '--threads', '8'], cwd=input_dir, capture_output=True)
    #log = subprocess.run(['conda', 'run', '-n', 'mtbseq', 'MTBseq', '--step', 'TBfull', '--threads', '8', '--mincovf', '1', '--mincovr', '1', '--minphred20', '1', '--lowfreq_vars', '--minfreq', '5'], cwd=input_dir, capture_output=True)
    with open(os.path.join(output_dir, 'mtbseq_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))
    # Remove input fastq files:
    for file_to_remove in file:
        subprocess.run(['rm', file_to_remove], cwd=os.path.join(output_dir, 'mtbseq_raw'))
        print(f'\t\tRemoved {file_to_remove} from session folder')
    copy_rename_result(os.path.join(input_dir, 'Called'), 'variants.+tab$', output_dir, 'mtbseq_susceptibility')
    copy_rename_result(os.path.join(input_dir, 'Classification'), '.tab$', output_dir, 'mtbseq_classification')
    copy_rename_result(os.path.join(input_dir, 'Statistics'), '.tab$', output_dir, 'mtbseq_statistics')
    copy_rename_result(os.path.join(input_dir, 'Position_Tables'), '.tab$', output_dir, 'mtbseq_pos_tables')
    
    print(f"\tMTBseq analysis finished at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}, elapsed time: {datetime.now() - datetime.strptime(time, '%d-%m-%Y, %H:%M:%S')}")
    if log.returncode == 0:
        stderr = 'STDERR empty'
    else:
        stderr = log.stderr.decode('utf-8')
    return ' '.join(log.args), str(log.returncode), stderr

def create_snakefile(input_dir, file, output_dir):
    source_config = '/data/leuven/348/vsc34807/software/gentb-snakemake/data/config.yaml'
    source_snakefile = '/data/leuven/348/vsc34807/software/gentb-snakemake/Snakefile'
    source_scripts = '/data/leuven/348/vsc34807/software/gentb-snakemake/scripts/'
    # These are original files downloaded from GenTB github repository, except for:
    # 1) Correction of indentation error in line 382 of varMatchUnk.py
    # 2) Change of hardcoded path in >>> load(“data/predict_rdate/pza_finalpredict_v2_0.RData”) to >>> load(paste(data_dir, “/pza_finalpredict_v2_0.RData”, sep=””)

    # Create copy of config.yaml file that will be modified
    copy_config = os.path.join(output_dir, 'gentb_raw', 'config.yaml')
    copy_snakefile = os.path.join(output_dir, 'gentb_raw', 'Snakefile')
    results_dir = os.path.join(output_dir, 'gentb_raw', 'results')
    os.mkdir(results_dir)
    shutil.copy(source_config, copy_config)
    shutil.copy(source_snakefile, copy_snakefile)
    #print(f'results dir: {results_dir}')
    for line in fileinput.input(copy_config, inplace=1):
        line = line.replace('data/', '/data/leuven/348/vsc34807/software/gentb-snakemake/data/')
        sys.stdout.write(line)
    for line in fileinput.input(copy_snakefile, inplace=1):
        line = line.replace('glob_wildcards("data/fastq/{sample}_1.fastq.gz")', f'["{file}"]')
        line = line.replace('data/config.yaml', copy_config)
        line = line.replace('data/fastq/{sample}_1.fastq.gz', f'{input_dir}/{{sample}}_R1.fastq.gz')
        line = line.replace('data/fastq/{sample}_2.fastq.gz', f'{input_dir}/{{sample}}_R2.fastq.gz')
        line = line.replace('scripts/', source_scripts)
        line = line.replace('results', results_dir)
        line = line.replace('-t 1', '-t 8')
        line = line.replace('--threads 1', '--threads 8')
        sys.stdout.write(line)
    print('\t\tCreated Snakefile')

def gentb(input_dir, file, args, output_dir):
    file_name = '_'.join(file[0].split('.')[0].split('_')[:-1])
    time = datetime.now().strftime('%d-%m-%Y, %H:%M:%S')
    raw_output_folder = os.path.join(output_dir, 'gentb_raw')
    os.mkdir(raw_output_folder)
    print(f"\tGenTB analysis started at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}")
    create_snakefile(input_dir, file_name, output_dir)
    log = subprocess.run(['conda', 'run', '-n', 'gentb', 'snakemake'], cwd=os.path.join(output_dir, 'gentb_raw'), capture_output=True)
    with open(os.path.join(output_dir, 'gentb_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))
    copy_rename_result(os.path.join(raw_output_folder, 'results'), 'predict.json$', output_dir, 'gentb')
    print(f"\tGenTB analysis finished at {datetime.now().strftime('%d-%m-%Y, %H:%M:%S')}, elapsed time: {datetime.now() - datetime.strptime(time, '%d-%m-%Y, %H:%M:%S')}")
    if log.returncode == 0:
        stderr = 'STDERR empty'
    else:
        stderr = log.stderr.decode('utf-8')
    return ' '.join(log.args), str(log.returncode), stderr


def mtbseq_results(sample_name, path):
    # df_res will contain information on variants identified per drug
    df_res = pd.DataFrame()
    # get filename for susceptibility, phylogenetics and qc results, using 'mtbseq_susceptibility', 'mtbseq_qc', ...
    # For mtbseq this is variable, since name depends on thresholds for coverage, depth, phredscore etc: use regex to find file that contains 'mtbseq_susceptibility'
    # find files at source destination 'path'
    result_list = os.listdir(path)
    result = [file for file in result_list if re.search('mtbseq_susceptibility', file)][0] # Pick the first file from the list
    path_res = os.path.join(path, result)
    path_phylo = os.path.join(path, 'mtbseq_classification_result')
    path_qc = os.path.join(path, 'mtbseq_statistics_result')

    # Check whether result file is available: if not exit function and return empty dataframes (df_res)
    try: df_all_variants = pd.read_csv(path_res, delimiter='\t')
    except IOError: return df_res, df_res, df_res

    ### Find resistance-associated variants
    ### Select rows in df_all_variants that are not empty in 'ResistanceSNP'
    resistance_associated_variants = df_all_variants.loc[df_all_variants['ResistanceSNP'] != ' ']    
    # First list drugs for which a resistance-conferring mutation was identified; deduplicate
    drugs_with_variants = set(list(resistance_associated_variants.ResistanceSNP))
    # Yields: {'fluoroquinolones (FQ)', 'isoniazid (INH)', 'ethionamide (ETH)', 'rifampicin (RMP)', 'ethambutol (EMB)', 'streptomycin (SM)', ' ', 'amikacin (AMK) kanamycin (KAN) capreomycin (CPR)'}
    drugs_with_variants_ungrouped = []
    # ungroup to get aminoglycosides and INH/ETH separately, 
    for group in drugs_with_variants:  # re.findall produces: ['FQ'] ['INH'] ['ETH'] ['RMP'] ['EMB'] ['SM'] [] ['AMK', 'KAN', 'CPR']
        drugs_with_variants_ungrouped.extend(re.findall('\((.*?)\)', group))
    # drugs_with_variants_ungrouped yields: ['FQ', 'INH', 'ETH', 'RMP', 'EMB', 'SM', 'AMK', 'KAN', 'CPR']
    for drug in drugs_with_variants_ungrouped:
        # select all rows in Gene/GeneName, Subst, Freq, per drug for which variants were found
        variants_per_drug = resistance_associated_variants.loc[resistance_associated_variants['ResistanceSNP'].str.contains(drug)]
        row = {} # row is a dictionary contains all information on detected variants per drug: create column per gene
        row['sample'] = sample_name
        row['tool'] = 'mtbseq'
        row['drug'] = drug_dict[drug]
        row['cat'] = 'R'
        # Iterate over rows of variants_per_drug:
        for index, variant in variants_per_drug.iterrows():
            if variant['GeneName'] in row: # if a variant was already added for this particular gene: append
                row[variant['GeneName']].append(f"{variant['Subst']} ({variant['Type']})")
                row[variant['GeneName']].append(round((variant['Freq']) / 100, 3))
                F_R_coverage = [variant['CovFor'], variant['CovRev']]
                row[variant['GeneName']].append(F_R_coverage)
            else: # initialize the entry in row per gene as a list, to be able to append multiple variants per gene
                row[variant['GeneName']] = []
                row[variant['GeneName']].append(f"{variant['Subst']} ({variant['Type']})")
                row[variant['GeneName']].append(round((variant['Freq']) / 100, 3))
                F_R_coverage = [variant['CovFor'], variant['CovRev']]
                row[variant['GeneName']].append(F_R_coverage)
            # Convert row dictionary to dataframe and concatenate with previous rows
            new_row = pd.DataFrame.from_dict(row, orient="index")
            new_row = new_row.transpose()
        df_res = pd.concat([df_res, new_row])

    ### Select rows in df_all_variants that are not empty in 'InterestingRegion', but that are empty in 'ResistanceSNP'
    variants_unknown_significance = df_all_variants.loc[(df_all_variants['InterestingRegion'] != ' ') & (df_all_variants['ResistanceSNP'] == ' ')]
    # First list drugs for which a resistance-conferring mutation was identified; deduplicate
    drugs_with_variants = set(list(variants_unknown_significance.InterestingRegion))
    # Yields: {'fluoroquinolones (FQ)', 'isoniazid (INH)', 'ethionamide (ETH)', 'rifampicin (RMP)', 'ethambutol (EMB)', 'streptomycin (SM)', ' ', 'amikacin (AMK) kanamycin (KAN) capreomycin (CPR)'}
    drugs_with_variants_ungrouped = []
    for group in drugs_with_variants:  # re.findall produces: ['FQ'] ['INH'] ['ETH'] ['RMP'] ['EMB'] ['SM'] [] ['AMK', 'KAN', 'CPR']
        drugs_with_variants_ungrouped.extend(re.findall('\((.*?)\)', group)) # drugs_with_variants_ungrouped yields: ['FQ', 'INH', 'ETH', 'RMP', 'EMB', 'SM', 'AMK', 'KAN', 'CPR']
    for drug in drugs_with_variants_ungrouped:
        # select all rows in Gene/GeneName, Subst, Freq, per drug for which variants were found but exclude those that are known to confer resistance
        variants_per_drug = variants_unknown_significance.loc[
            variants_unknown_significance['InterestingRegion'].str.contains(drug)]
        row = {}
        row['sample'] = sample_name
        row['tool'] = 'mtbseq'
        row['drug'] = drug_dict[drug]
        row['cat'] = '?'
        # Iterate over rows of variants_per_drug:
        for index, variant in variants_per_drug.iterrows():
            if variant['GeneName'] in row:
                row[variant['GeneName']].append(f"{variant['Subst']} ({variant['Type']})")
                row[variant['GeneName']].append(round((variant['Freq']) / 100, 3))
                F_R_coverage = [variant['CovFor'], variant['CovRev']]
                row[variant['GeneName']].append(F_R_coverage)
            else:
                row[variant['GeneName']] = []
                row[variant['GeneName']].append(f"{variant['Subst']} ({variant['Type']})")
                row[variant['GeneName']].append(round((variant['Freq']) / 100, 3))
                F_R_coverage = [variant['CovFor'], variant['CovRev']]
                row[variant['GeneName']].append(F_R_coverage)
            new_row = pd.DataFrame.from_dict(row, orient="index")
            new_row = new_row.transpose()
        df_res = pd.concat([df_res, new_row])

    ### Get phylogenetic information
    df_phylo = pd.read_csv(path_phylo, delimiter='\t')
    df_phylo = df_phylo.transpose()
    df_phylo.index.name = 'parameter'
    df_phylo['tool'] = 'MTBseq'
    df_phylo.set_index(['tool'], inplace=True, append=True)
    df_phylo.index = df_phylo.index.reorder_levels(['tool', 'parameter'])

    ### Get qc information
    df_qc = pd.read_csv(path_qc, delimiter='\t')
    df_qc = df_qc.transpose()
    df_qc.index.name = 'parameter'
    df_qc['tool'] = 'MTBseq'
    df_qc.set_index(['tool'], inplace=True, append=True)
    df_qc.index = df_qc.index.reorder_levels(['tool', 'parameter'])
    return df_res, df_phylo, df_qc

def mykrobe_results(sample_name, path):
    # df will contain all info on susceptibility, phylogenetics and qc
    path = os.path.join(path, 'mykrobe_result')
    
    df = pd.DataFrame()
    try:
        with open(path, 'r') as json_load:
            result = json.load(json_load)
    except IOError:
        return df, df

    # Highest level of Mykrobe output .json is sample name without .fastq.gz extension
    for drug, prediction in result[''.join(sample_name.split('.')[:-2])]['susceptibility'].items():
        row = {}
        if prediction['predict'] == 'R':
            for gene, info in prediction['called_by'].items():
                row['sample'] = sample_name
                row['tool'] = 'mykrobe'
                row['drug'] = drug_dict[drug.lower()]
                row['cat'] = prediction['predict']

                gene_ = gene.split('_')[0]  # gene is e.g. gyrA_D94A-GAC7581GCC
                variant = gene.split('_', 1)[1]
                coverage = info['info']['coverage']['alternate']['percent_coverage']
                depth = info['info']['coverage']['alternate']['median_depth']
                if gene_ in row:
                    row[gene_].append(variant)
                    row[gene_].append(int(coverage) / 100)
                    row[gene_].append(depth)
                else:
                    row[gene_] = []
                    row[gene_].append(variant)
                    row[gene_].append(int(coverage) / 100)
                    row[gene_].append(depth)
            new_row = pd.DataFrame.from_dict(row, orient="index")
            new_row = new_row.transpose()
            df = pd.concat([df, new_row])

    # Phylogenetic info
    phylogenetic_info = result[''.join(sample_name.split('.')[:-2])]['phylogenetics']
    row = {}
    row['phylo_group'] = phylogenetic_info['phylo_group']
    row['subcomplex'] = phylogenetic_info['sub_complex']
    row['species'] = phylogenetic_info['species']
    row['lineage'] = phylogenetic_info['lineage']['lineage']
    df_phylo = pd.DataFrame.from_dict(row.items())
    df_phylo.rename(columns={0: 'parameter'}, inplace=True)
    df_phylo['tool'] = 'Mykrobe'
    df_phylo.set_index(['tool', 'parameter'], inplace=True)
    df_phylo.rename(columns={1: 0}, inplace=True)

    return df, df_phylo


def gentb_results(sample_name, path):
    df = pd.DataFrame()
    path = os.path.join(path, 'gentb_result')
    
    try:
        with open(path, 'r') as json_load:
            result = json.load(json_load)
    except IOError:
        return df
    # result[0] is a list of resistance-associated variants, e.g. [[variant_1], [variant_2], [variant_3], … ]
    # result[1] is a dictionary with as key the filepath, and then a list of five lists of X variants detected [[variant1, variant2, variant3, …], [variant1b, …]]
    # 'mutations' contains this list of five lists
    
    try: # To catch KeyErrors when output json is empty or inconsistent results outputted (such as resistance predicted but no variants reported)

        mutations = result[1][list(result[1].keys())[0]]
        for index, res_prob in enumerate(result[0]):
            row = {}
            if res_prob[2] == "1": 
                row['cat'] = 'R'
                for i in range(5):
                    if mutations[i][index]:
                        genes = mutations[i][index].rsplit('_', 1)[1].rsplit('.')  # split '_' and '.': fabG1 and inhA reported together in GenTB but not for other tools
                        variant = mutations[i][index].rsplit('_', 1)[0]
                        for gene in genes:
                            row[gene] = str([variant, np.NaN, np.NaN])

                row['drug'] = drug_dict[res_prob[1]]
                row['tool'] = 'gentb'
                row['sample'] = sample_name
                row['prob'] = res_prob[5]
                new_row = pd.DataFrame.from_dict(row, orient="index")
                new_row = new_row.transpose()
                df = pd.concat([df, new_row])
        print(f'gentb df na concat: \n {df}')
        # Loop through variants with unknown phenotype list (in results[2]) and add all
        row = {}
        mutations = result[2][list(result[2].keys())[0]]
        for i in range(5):
            for mutation in mutations[i]:
                if mutation:
                    genes = mutation.rsplit('_', 1)[1].rsplit('.')
                    variant = mutation.rsplit('_', 1)[0]
                    for gene in genes:
                        row[gene.rsplit('-')[-1]] = str([variant, np.NaN, np.NaN])
                    row['drug'] = gene_dict[gene.rsplit('-')[-1]]
                    row['tool'] = 'gentb'
                    row['sample'] = sample_name
                    row['cat'] = '?'
                    row['prob'] = '?'
                    new_row = pd.DataFrame.from_dict(row, orient="index")
                    new_row = new_row.transpose()
                    df = pd.concat([df, new_row])
        print(f'\ngentb df ? na concat: \n {df}')
    except Exception as e: print(f'\v\t GenTB data extraction failed: {e}')
    return df

def tbprofiler_results(sample_name, path):
    df = pd.DataFrame()
    path = os.path.join(path, 'tbprofiler_result')

    try:
        with open(path, 'r') as json_load:
            result = json.load(json_load)
    except IOError:
        return df, df, df

    # loop through 'dr_variants' entry and then through 'drugs' to identify the drugs with detected mutations
    drugs_with_variants = []
    for dr_variant in result['dr_variants']:
        for confers_resistance_to in dr_variant['drugs']:
            drugs_with_variants.append(confers_resistance_to['drug'])
    # eliminate duplicates
    drugs_with_variants = set(drugs_with_variants)

    # collect data on variants per drug
    for drug in drugs_with_variants:
        row = {}
        for dr_variant in result['dr_variants']:
            for confers_resistance_to in dr_variant['drugs']:
                if drug == confers_resistance_to['drug']:
                    row['sample'] = sample_name
                    row['tool'] = 'tbprofiler'
                    row['drug'] = drug_dict[drug]
                    row['cat'] = 'R'  # or: confers_resistance_to['confers']
                    if dr_variant['gene'] in row:
                        row[dr_variant['gene']].append(
                            f"{dr_variant['change']} ({dr_variant['nucleotide_change']}, {dr_variant['type']})")
                        row[dr_variant['gene']].append(round(dr_variant['freq'], 3))
                        row[dr_variant['gene']].append(result['qc']['median_coverage'])
                    else:
                        row[dr_variant['gene']] = []
                        row[dr_variant['gene']].append(
                            f"{dr_variant['change']} ({dr_variant['nucleotide_change']}, {dr_variant['type']})")
                        row[dr_variant['gene']].append(round(dr_variant['freq'], 3))
                        row[dr_variant['gene']].append(result['qc']['median_coverage'])
        new_row = pd.DataFrame.from_dict(row, orient="index")
        new_row = new_row.transpose()
        df = pd.concat([df, new_row])

    # collect data on variants with unknown phenotype
    for other_variant in result['other_variants']:
        row = {}
        row['sample'] = sample_name
        row['tool'] = 'tbprofiler'
        if other_variant['gene'] in gene_dict:
            row['drug'] = gene_dict[other_variant['gene']]
        else:
            row['drug'] = '?'
        row['cat'] = '?'
        if other_variant['gene'] in row:
            row[other_variant['gene']].append(
                f"{other_variant['change']} ({other_variant['nucleotide_change']}, {other_variant['type']})")
            row[other_variant['gene']].append(round(other_variant['freq'], 3))
            row[other_variant['gene']].append(result['qc']['median_coverage'])
        else:
            row[other_variant['gene']] = []
            row[other_variant['gene']].append(
                f"{other_variant['change']} ({other_variant['nucleotide_change']}, {other_variant['type']})")
            row[other_variant['gene']].append(round(other_variant['freq'], 3))
            row[other_variant['gene']].append(result['qc']['median_coverage'])
        new_row = pd.DataFrame.from_dict(row, orient="index")
        new_row = new_row.transpose()
        df = pd.concat([df, new_row])

    # Loop through entries under Lineage section: find entries with most information
    row = {'lineage': '', 'family': '', 'spoligotype': '', 'rd': ''}
    for entry in result['lineage']:
        if len(entry['lin']) > len(row['lineage']): row['lineage'] = entry['lin']
        if len(entry['family']) > len(row['family']): row['family'] = entry['family']
        if len(entry['spoligotype']) > len(row['spoligotype']): row['spoligotype'] = entry['spoligotype']
        if len(entry['rd']) > len(row['rd']): row['rd'] = entry['rd']

    df_phylo = pd.DataFrame.from_dict(row.items())
    df_phylo['tool'] = 'TBProfiler'
    df_phylo.rename(columns={0: 'parameter'}, inplace=True)
    df_phylo.set_index(['tool', 'parameter'], inplace=True)
    df_phylo.rename(columns={1: 0}, inplace=True)

    # Extract QC info from result['qc']
    row = {}
    row['reads_mapped'] = result['qc']['pct_reads_mapped']
    row['median_coverage'] = result['qc']['median_coverage']
    row['d_r_type'] = result['drtype']
    df_qc = pd.DataFrame.from_dict(row.items())
    df_qc.set_index(0, inplace=True)
    df_qc['tool'] = 'TBProfiler'
    df_qc.index.name = 'parameter'
    df_qc.rename(columns={1: 0}, inplace=True)
    df_qc.set_index(['tool'], inplace=True, append=True)
    df_qc.index = df_qc.index.reorder_levels(['tool', 'parameter'])
    return df, df_phylo, df_qc


drug_list = ['INH', 'RIF', 'EMB', 'PYR', 'STR', 'CAP', 'KAN', 'AMK', 'FLQ', 'OFL', 'LEV', 'MOX', 'CIP', 'GAT', 'ETH', 'PRO', 'CYC', 'TER', 'PAS', 'BED', 'DEL', 'LIN', 'CLO', 'AMC']
tool_list = ['tbprofiler', 'mykrobe', 'gentb', 'mtbseq']
pred_drugs_per_tool = {'mykrobe': ['INH', 'RIF', 'EMB', 'PYR', 'STR', 'CAP', 'KAN', 'AMK', 'OFL', 'MOX', 'CIP'], # Derived from Table 1 Hunt et al. 2019
                      'tbprofiler': ['INH', 'OFL', 'PYR', 'STR', 'CLO', 'LIN', 'AMK', 'FLQ', 'PAS', 'EMB', 'AGL', 'CYC', 'BED', 'CAP', 'ETH', 'RIF', 'DEL', 'KAN', 'LEV', 'CIP', 'MOX'], # Derived from tbdb.csv on TBProfiler github,
                      'gentb': ['INH', 'RIF', 'PYR', 'STR', 'EMB', 'ETH', 'KAN', 'CAP', 'AMK', 'LEV', 'OFL', 'PAS'], # From Gröschel et al 2021 - one drug missing - predictions to 13 drugs with RF?
                      'mtbseq': ['AMK', 'KAN', 'CAP', 'PAS', 'EMB', 'PYR', 'STR', 'LIN', 'KAN', 'ETH', 'FLQ', 'CAP', 'INH', 'RIF']} # Antibiotic column from Github MTB_Resistance_Mediating.txt

aa_dict = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'B': 'Asx', 'C': 'Cys', 'Q': 'Gln',
          'E': 'Glu', 'Z': 'Glx', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
          'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr',
          'V': 'Val', 'X': 'X'}

def format_variants(drug_tool_with_variants):
    # pop off columns with sample id, tool, drug and categorization: put into info_columns
    info_columns = drug_tool_with_variants[['sample', 'tool', 'drug', 'cat', 'prob']]
    drug_tool_with_variants = drug_tool_with_variants.drop(['sample', 'tool', 'drug', 'cat', 'prob'])
    # loop through columns with genes and variants: if np.NaN: change into 'ND', else: loop through list items (individual variants)
    variants_per_gene = []
    for gene, variant in drug_tool_with_variants.iteritems():
        if isinstance(variant, float): # to detect values of gene_column np.NaN: these are float; info is stored as lists (all tools except GenTB) or string-converted lists (GenTB)
            pass
        else:
            if info_columns['tool'] == 'tbprofiler': 
                for i in range(int(len(variant) / 3)): # parse variant per triplets of list-items (triplet 0, 1, 2, …)
                    triplet = variant[i*3:i*3+3] # e.g. ['c.-15C>T (c.-15C>T, upstream_gene_variant)', 1.0, 25]
                    variant_simple = f'{gene}, {triplet[0].split(".")[1].split(" ")[0]} (freq {triplet[1] * 100}%, med_cov genome {triplet[2]})'

            if info_columns['tool'] == 'mykrobe':
                for i in range(int(len(variant) / 3)):
                    triplet = variant[i*3:i*3+3]
                    variant_mykrobe = triplet[0].split("-")[:-1] # ['I194T', 'ATC1674781ACC'] or ['C', '15X', 'C1673425T'] => [I194T] or ['C', '15X']
                    if len(variant_mykrobe) == 1: # In case of [I194T]: turn letters into amino acids
                        aa_start = aa_dict[variant_mykrobe[0][0]]
                        aa_stop = aa_dict[variant_mykrobe[0][-1]]
                        var_pos = variant_mykrobe[0][1:-1]
                        res = f'{aa_start}{var_pos}{aa_stop}'
                    else: # In case of ['C', '15X']: reconstruct to '-C-15X'
                        res = f"-{'-'.join(variant_mykrobe)}"
                    variant_simple = f'{gene}, {res} (freq {triplet[1] * 100}%, med_cov var {triplet[2]})'

            if info_columns['tool'] == 'mtbseq':
                for i in range(int(len(variant) / 3)):
                    triplet = variant[i*3:i*3+3]
                    variant_substring = triplet[0].split(" ")[0] # ['Ile194Thr (atc/aCc) (SNP)'] => Ile194Thr
                    if variant_substring == '': # sometimes mtbseq marks variants without naming a gene or specific variant: these are gene '' or '-'
                        gene = 'undef'
                        variant_substring = 'undef'
                    variant_simple = f'{gene}, {variant_substring} (freq {triplet[1] * 100}%, cov F/R {triplet[2]})'

            if info_columns['tool'] == 'gentb':
                # results from gentb are string representations of lists; this is because dataframe pd.from_dict can't handle lists
                # convert string representation to list:
                variant_list = [] # convert string representation of list to list, remove brackets and ''
                [variant_list.append(variant.replace("'", "")) for variant in variant.strip('][').split(', ')]
                for i in range(int(len(variant_list) / 3)):
                    triplet = variant_list[i*3:i*3+3]
                    variant_gentb = triplet[0].split('_')[3:] # ['SNP_CN_761998_T2192C_L731P', nan, nan] => T2192C
                    #aa_start = aa_dict[variant_gentb[0]]
                    #aa_stop = aa_dict[variant_gentb[-1]]
                    #var_pos = variant_gentb[1:-1]
                    #res = f'{aa_start}{var_pos}{aa_stop}'
                    res = '_'.join(variant_gentb)
                    variant_simple = f'{gene}, {res} (NA)'
            variants_per_gene.append(variant_simple)
    variants_per_gene_string = '- ' + '\n- '.join(variants_per_gene)
    return variants_per_gene_string

def assign_regimen(summary_table):
    regimen_dict_all_tools = {} # will contain per tool the assigned treatment regimen
    regimen_legend_list = [] # is a list containing different predicted regimens across four tools; this will later be deduplicated
    regimen_details = {} # a dictionary that contains regimen def and drugs per tool
    legend_list_dicts = [] # this is a list of regimen_details for each regimen to create the legend
   
    # GenTB doesn't report categorization for fluoroquinolones as a group, neither does it predict moxifloxacin categorization
    # MTBseq reports quinolones as a group,  but does not predict moxifloxacin categorizatoin
    # the who_regimen.py script however uses moxifloxacin categorization for defining XDR-TB: therefore add to summary_table for GenTB and MTBseq an entry for mox if FLQ or CIP reported R
    if summary_table['gentb_CIP'] == 'R': summary_table['gentb_MOX'] = 'R'
    if summary_table['mtbseq_FLQ'] == 'R': summary_table['mtbseq_MOX'] = 'R'
 
    # use mykrobe_reg_dict to convert drug abbreviations to capitalized full drug names used in Mykrobe who_regimen.py
    mykrobe_reg_dict = {
    'INH': 'Isoniazid',
    'RIF': 'Rifampicin',
    'PYR': 'Pyrazinamide',
    'EMB': 'Ethambutol',
    'KAN': 'Kanamycin',
    'AMK': 'Amikacin',
    'CAP': 'Capreomycin',
    'STR': 'Streptomycin',
    'OFL': 'Ofloxacin',
    'CIP': 'Ciprofloxacin',
    'MOX': 'Moxifloxacin',
    'BED': 'Bedaquiline',
    'LIN': 'Linezolid',
    'CLO': 'Clofazimide',
    'CYC': 'Cycloserine',
    'TER': 'Terizidone'
    }

    for tool in ('mykrobe', 'tbprofiler', 'mtbseq', 'gentb'):
        reg_dict_per_tool = {} # contains predictions per tool and drug formatted as: {'Isoniazid': 'R'/'S'/None, …}
        for tool_drug, cat in summary_table.items():
            if (tool in tool_drug) & (tool_drug.split('_')[1] in mykrobe_reg_dict.keys()): 
                formatted_drug = mykrobe_reg_dict[tool_drug.split('_')[1]]
                formatted_cat = 'R' if cat == 'R' else 'S' if cat == 'ND' else 'S' if cat == '?' else None # in case of variant of uncertain significance: still predict S for regimen prediction
                reg_dict_per_tool[formatted_drug] = formatted_cat
        # input formatted predictions per tool in DstProfile class
        pred_reg = DstProfile(reg_dict_per_tool)
        # Get 'definition' of predicted regimen, which is 'DS-TB', 'XDR-TB', etc
        regimen_dict_all_tools[f'{tool}_regimen'] = [*pred_reg.regimen][1]
        # add all different predicted regiments to a list to create a legend
        regimen_legend_list.append(str([*pred_reg.regimen][0])) 
    # Deduplicate this list
    regimen_legend = set(regimen_legend_list)
    # for all unique regimens, append them to a list of dictionaries for parameters regimen definition, mandatory and optional drugs
    for regimen_no in regimen_legend:
        # look up regimen_no in who_treatment.py dictionary 'regimens' and get attributes; put into regimen_details
        regimen_details['regimen_def'] = str(regimens[int(regimen_no)][1])
        regimen_details['regimen_mandatory'] = str(regimens[int(regimen_no)][2])
        regimen_details['regimen_optional'] = str(regimens[int(regimen_no)][3])
        legend_list_dicts.append(regimen_details.copy())
    return regimen_dict_all_tools, legend_list_dicts                
    

def get_susceptibility_dicts(df_susceptibility):
    summary_table = {}
    prob_dict = {}
    overview_table = []
    for drug in drug_list:
        new_row_overview_table_R = {}
        new_row_overview_table_unk = {}
        for tool in tool_list:
            # filter rows based on presence of present drug in df
            df_drug = df_susceptibility.loc[(df_susceptibility['drug'].str.contains(drug)) & (df_susceptibility['tool'].str.contains(tool))]
            if df_drug.empty: # if entry doesn't exist for this tool-drug combination: set to 'not detected' (ND) or dash (-) if no prediction
                if drug in pred_drugs_per_tool[tool]:
                    summary_table[f'{tool}_{drug}'] = 'ND'
                else:
                    summary_table[f'{tool}_{drug}'] = '-'
            else: # here are both 'R' and '?' predictions
                # loop through df_drug: put 'R' predictions in summary_table
                for index, row in df_drug.iterrows():
                    if 'R' in row['cat']:
                        summary_table[f'{tool}_{drug}'] = 'R'
                        new_row_overview_table_R[f'overview_drug'] = f'{drug} (R)'
                        # get list of formatted variants
                        formatted_row = format_variants(row)
                        if tool == 'gentb':
                            prob_dict[f'prob_{drug}'] = int(float(row['prob'])*100)
                        if tool == 'tbprofiler': # view a simplified format of TBProfiler-detected variants in summary_table
                            formatted_row_simple_raw = formatted_row.split('- ')
                            formatted_row_simple = '\n'.join([var.split(' (')[0] for var in formatted_row_simple_raw if var != ''])
                            summary_table[f'variants_{drug}'] = formatted_row_simple
                        new_row_overview_table_R[f'overview_{tool}'] = formatted_row
                    elif '?' in row['cat']: 
                        if f'{tool}_{drug}' not in summary_table.keys(): # if no 'R' prediction but still variant of unknown significance detected: report '?'
                            summary_table[f'{tool}_{drug}'] = '?'
                        new_row_overview_table_unk[f'overview_drug'] = f'{drug} (?)'
                        # get list of formatted variants
                        formatted_row = format_variants(row)
                        new_row_overview_table_unk[f'overview_{tool}'] = formatted_row

        if new_row_overview_table_R:
            overview_table.append(new_row_overview_table_R)
        if new_row_overview_table_unk:
            overview_table.append(new_row_overview_table_unk)
    return summary_table, overview_table, prob_dict

def get_phylo_dict(df_phylo):
    phylo_dict = {}
    phylo_dict['species'] = str(*df_phylo.loc['Mykrobe', 'species'].to_dict()['result'].keys())
    phylo_dict['species_percent_coverage'] = str(df_phylo.loc['Mykrobe', 'species'].to_dict()['result'][phylo_dict['species']]['percent_coverage'])
    phylo_dict['species_median_depth'] = str(df_phylo.loc['Mykrobe', 'species'].to_dict()['result'][phylo_dict['species']]['median_depth'])

    phylo_dict['complex'] = str(*df_phylo.loc['Mykrobe', 'phylo_group'].to_dict()['result'].keys())
    phylo_dict['complex_percent_coverage'] = str(df_phylo.loc['Mykrobe', 'phylo_group'].to_dict()['result'][phylo_dict['complex']]['percent_coverage'])
    phylo_dict['complex_median_depth'] = str(df_phylo.loc['Mykrobe', 'phylo_group'].to_dict()['result'][phylo_dict['complex']]['median_depth'])

    phylo_dict['subcomplex'] = str(*df_phylo.loc['Mykrobe', 'subcomplex'].to_dict()['result'].keys())
    phylo_dict['subcomplex_percent_coverage'] = str(df_phylo.loc['Mykrobe', 'subcomplex'].to_dict()['result'][phylo_dict['subcomplex']]['percent_coverage'])
    phylo_dict['subcomplex_median_depth'] = str(df_phylo.loc['Mykrobe', 'subcomplex'].to_dict()['result'][phylo_dict['subcomplex']]['median_depth'])

    phylo_dict['lineage'] = str(*df_phylo.loc['TBProfiler', 'lineage'].values).split('lineage')[1]
    phylo_dict['family'] = str(*df_phylo.loc['TBProfiler', 'family'].values)
    phylo_dict['spoligotype'] = str(*df_phylo.loc['TBProfiler', 'spoligotype'].values)
    phylo_dict['rd'] = str(*df_phylo.loc['TBProfiler', 'rd'].values)

    return phylo_dict

def get_qc_dict(df_qc):
    qc_dict = {}
    qc_dict['total_reads'] = str(*df_qc.loc['MTBseq','Total Reads'].values)[1:]
    qc_dict['mapped_reads'] = str(*df_qc.loc['MTBseq','Mapped Reads'].values)[1:]
    qc_dict['per_mapped_reads'] = str(*df_qc.loc['MTBseq','% Mapped Reads'].values)[1:]
    qc_dict['genome_size'] = str(*df_qc.loc['MTBseq','Genome Size'].values)[1:]
    qc_dict['genome_gc'] = str(*df_qc.loc['MTBseq','Genome GC'].values)[1:]
    qc_dict['any_total_bases'] = str(*df_qc.loc['MTBseq','(Any) Total Bases'].values)[1:]
    qc_dict['per_any_total_bases'] = str(*df_qc.loc['MTBseq','% (Any) Total Bases'].values)[1:]
    qc_dict['any_gc_content'] = str(*df_qc.loc['MTBseq','(Any) GC-Content'].values)[1:]
    qc_dict['any_coverage_mean'] = str(*df_qc.loc['MTBseq','(Any) Coverage mean'].values)[1:]
    qc_dict['any_coverage_median'] = str(*df_qc.loc['MTBseq','(Any) Coverage median'].values)[1:]
    qc_dict['unambiguous_total_bases'] = str(*df_qc.loc['MTBseq','(Unambiguous) Total Bases'].values)[1:]
    qc_dict['per_unambiguous_total_bases'] = str(*df_qc.loc['MTBseq','% (Unambiguous) Total Bases'].values)[1:]
    qc_dict['unambiguous_gc_content'] = str(*df_qc.loc['MTBseq','(Unambiguous) GC-Content'].values)[1:]
    qc_dict['unambiguous_coverage_mean'] = str(*df_qc.loc['MTBseq','(Unambiguous) Coverage mean'].values)[1:]
    qc_dict['unambiguous_coverage_median'] = str(*df_qc.loc['MTBseq','(Unambiguous) Coverage median'].values)[1:]
    qc_dict['snps'] = str(*df_qc.loc['MTBseq','SNPs'].values)[1:]
    qc_dict['deletions'] = str(*df_qc.loc['MTBseq','Deletions'].values)[1:]
    qc_dict['insertions'] = str(*df_qc.loc['MTBseq','Insertions'].values)[1:]
    qc_dict['uncovered'] = str(*df_qc.loc['MTBseq','Uncovered'].values)[1:]
    qc_dict['substitutions'] = str(*df_qc.loc['MTBseq','Substitutions (Including Stop Codons)'].values)[1:]
    qc_dict['d_r_type'] = str(*df_qc.loc['TBProfiler', 'd_r_type'].values) 
    return qc_dict

def get_SNP_distances(output_dir, sample):
    phylo_dir = '/data/leuven/348/vsc34807/phylo'
    # get sample name: sample variable is file[0], which is e.g. 'ERR1664619_library-id_R1.fastq.gz': retain only identifier and library id
    sample_name = '_'.join(sample.split('_')[:2])
    # Copy mtbseq_pos_tables_result to [SampleID]_[LibID]_[*].gatk_position_table.tab, as required by MTBseq TBjoin module
    shutil.copy(os.path.join(output_dir, sample, 'mtbseq_pos_tables_result'), os.path.join(phylo_dir, 'Position_Tables', f'{sample_name}.gatk_position_table.tab'))
    result_list = os.listdir(os.path.join(output_dir, sample))
    result = [file for file in result_list if re.search('mtbseq_susceptibility', file)][0] # Pick the first file from the list
    info_string = '_'.join(result.split('_')[2:7]) # info_string contains info on mincovf, mincovr, minfreq, minphred contained in filename
    # Copy mtbseq_susceptibility_result to Called/[SampleID]_[LibID]_[*].gatk_position_variants_[mincovf]_[mincovr]_[minfreq]_[minphred20]_[all_vars][snp_vars][lowfreq_vars].tab
    shutil.copy(os.path.join(output_dir, sample, result), os.path.join(phylo_dir, 'Called', f'{sample_name}.gatk_position_variants_{info_string}.tab'))
    # create tab-separated file with sample id in first column and library id in second column, needed for TBjoin and TBamend: append to previous samples.txt, after deduplication
    append_sample = pd.DataFrame(sample.split('_')[0:2]).transpose() # this is new sample that we append
    samples_txt = pd.read_csv(os.path.join(phylo_dir, 'samples.txt'), sep='\t', header=None)
    new_samples_txt = pd.merge(append_sample, samples_txt, how='outer').drop_duplicates()
    new_samples_txt.to_csv(os.path.join(phylo_dir, 'samples.txt'), header=False, index=False, sep='\t')

    session_name = output_dir.split('/')[-1]
    # Call MTBseq TBjoin module: places files in Joint folder
    log = subprocess.run(['conda', 'run', '-n', 'mtbseq', 'MTBseq', '--step', 'TBjoin', '--samples', 'samples.txt', '--project', session_name], cwd=phylo_dir, capture_output=True)
    with open(os.path.join(output_dir, sample, 'TBjoin_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))

    # Call MTBseq Amend module: creates phylo.fasta file in Amend folder - --continue argument executes TBgroups subsequently
    log = subprocess.run(['conda', 'run', '-n', 'mtbseq', 'MTBseq', '--step', 'TBamend', '--samples', 'samples.txt', '--project', session_name, '--continue'], cwd=phylo_dir, capture_output=True)
    with open(os.path.join(output_dir, sample, 'TBamend_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))

    # Execute fasttree on session_name[*]phylo.fasta file; e.g. NONE_joint_cf4_cr4_fr75_ph4_samples4_amended_u95_phylo.fasta
    # Find this file in Amend folder
    fasta_files = os.listdir(os.path.join(phylo_dir, 'Amend'))
    fasta = [os.path.join(phylo_dir, 'Amend', file) for file in fasta_files if re.search(rf"{session_name}.+phylo.fasta", file)]
    latest_fasta = max(fasta, key=os.path.getctime) # get last created fasta file
    # Run fasttree on nucleotide data (-nt) and with GTR+CAT model (-gtr)
    fasta_path = os.path.join(phylo_dir, 'Amend', latest_fasta)
    with open(os.path.join(output_dir, sample, 'tree'), 'w') as outfile:
        subprocess.run(['fasttree', '-gtr', '-nt', fasta_path], stdout=outfile)

    # Render phylogenetic tree using separate script for relies on xvfb wrapper that cannot be called together with Conda ete3 environment call
    tree_path = os.path.join(output_dir, sample, 'tree')
    log = subprocess.run(['conda', 'run', '-n', 'ete3', 'python', '/data/leuven/348/vsc34807/software/render_tree.py', tree_path], capture_output=True)
    with open(os.path.join(output_dir, sample, 'ete_log.txt'), 'w') as logfile:
        logfile.write(log.stdout.decode('utf-8'))
        logfile.write(log.stderr.decode('utf-8'))

    # Read latest Groups file to get SNP distance matrix
    matrix_files = os.listdir(os.path.join(phylo_dir, 'Groups'))
    matrix = [os.path.join(phylo_dir, 'Groups', file) for file in matrix_files if re.search(rf"{session_name}.+.matrix", file)]
    latest_matrix = max(matrix, key=os.path.getctime) # get last created fasta file
    
    matrix_df = pd.DataFrame()
    column_list = [0]
    with open(latest_matrix, 'r') as file:
        for line in file:
            column_list.append(line.split('\t')[0])
            matrix_df = pd.concat([matrix_df, pd.DataFrame(line.split('\t')[:-1]).transpose()])
        # Add one additional line to df for last entry so no length mismatch with column names
        matrix_df = pd.concat([matrix_df, pd.DataFrame([np.NaN]*(len(matrix_df.columns)+1)).transpose()])
    matrix_df.columns = column_list[:-1]

    # remove index column; replace by '0' column containing sample numbers
    matrix_df.set_index(0, inplace=True)
    # threshold for SNP distance
    threshold = 12

    df_threshold = pd.DataFrame()
    index_sample = sample_name
    df_numeric = pd.to_numeric(matrix_df[index_sample])
    df_threshold = matrix_df.loc[df_numeric <= threshold, index_sample]

    SNP_distances = df_threshold.to_string()[1:]
    print(f'snp distances: {SNP_distances}')
    no_close_samples = len(df_threshold)
    print(f'no of close samples: {no_close_samples}')

    return SNP_distances, no_close_samples

    
def docx_mailmerge(summary_dict, overview_dict, phylo_dict, qc_dict, regimen_dict, legend_dict):

    template = '/data/leuven/348/vsc34807/software/report_template.docx'
    document = MailMerge(template)
    document.merge(**summary_dict)
    document.merge_rows('overview_drug', overview_dict)
    document.merge(**phylo_dict)
    document.merge(**qc_dict)
    document.merge(**regimen_dict)
    document.merge_rows('regimen_def', legend_dict)
    document.write(os.path.join(output_dir, f'{file[0].split(".")[0]}_report.docx'))

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, required=False)
parser.add_argument('-o', type=str, required=False)
parser.add_argument('--mykrobe', type=str, required=False)
parser.add_argument('--tbprofiler', type=str, required=False)
parser.add_argument('--mtbseq', type=str, required=False)
parser.add_argument('--gentb', type=str, required=False)
parser.add_argument('--light', type=str, required=False)
parser.add_argument('--phylo', type=str, required=False)

args=parser.parse_args()  

mykrobe_v, tbprofiler_v, mtbseq_v, gentb_v = start()
input_dir, formatted_paired_end, single_end = input_setup(str(args.i))
output_dir, analysis_time = output_setup(str(args.o))

# First proces paired-end files in list 'formatted_paired_end'
for index, file in enumerate(formatted_paired_end):

    file_output_dir = os.path.join(output_dir, file[0])
    os.mkdir(file_output_dir)
    print(f'\n File {index + 1}/{len(formatted_paired_end)}: {file}')
    if str(args.mykrobe) == 'false':
        print('\tSkipped Mykrobe analysis')
        mykrobe_return_code = '1'
        mykrobe_stderr = 'Mykrobe not executed'
        mykrobe_args = 'Mykrobe not executed'
    else:
        mykrobe_args, mykrobe_return_code, mykrobe_stderr = mykrobe(input_dir, file, str(args.mykrobe), file_output_dir)
        if mykrobe_return_code == '0':
            try: res_mykrobe, phylo_mykrobe = mykrobe_results(file[0], file_output_dir)
            except Exception as e: print(f'Mykrobe return code 0 but data extraction failed: {e}')

    if str(args.tbprofiler) == 'false':
        print('\tSkipped TBProfiler analysis')
        tbprofiler_return_code = '1'
        tbprofiler_stderr = 'TBProfiler not executed'
        tbprofiler_args = 'TBProfiler not executed'
    else:
        tbprofiler_args, tbprofiler_return_code, tbprofiler_stderr = tbprofiler(input_dir, file, str(args.tbprofiler), file_output_dir)
        if tbprofiler_return_code == '0':
            try: res_tbprofiler, phylo_tbprofiler, qc_tbprofiler = tbprofiler_results(file[0], file_output_dir)
            except Exception as e: print(f'TBProfiler return code 0 but data extraction failed: {e}')

    if str(args.gentb) == 'false':
        print('\tSkipped GenTB analysis')
        gentb_return_code = '1'
        gentb_stderr = 'GenTB not executed'
        gentb_args = 'GenTB not executed'
    else:
        gentb_args, gentb_return_code, gentb_stderr = gentb(input_dir, file, str(args.gentb), file_output_dir)
        if gentb_return_code == '0':
            try: res_gentb = gentb_results(file[0], file_output_dir)
            except Exception as e: print(f'GenTB return code 0 but data extraction failed: {e}')

    if str(args.mtbseq) == 'false':
        print('\tSkipped MTBseq analysis')
        mtbseq_return_code = '1'
        mtbseq_stderr = 'MTBseq not executed'
        mtbseq_args = 'MTBseq not executed'
    else:
        mtbseq_args, mtbseq_return_code, mtbseq_stderr = mtbseq(input_dir, file, str(args.mtbseq), file_output_dir)
        if mtbseq_return_code == '0':
            try: res_mtbseq, phylo_mtbseq, qc_mtbseq = mtbseq_results(file[0], file_output_dir)
            except Exception as e: print(f'GenTB return code 0 but data extraction failed: {e}')

    # Temporary try construction to catch KeyErrors and errors in docx_mailmerge when tools are disabled
    try:
        # Concatenate susceptibility dataframes
        df_susceptibility = pd.concat([res_tbprofiler, res_mykrobe, res_gentb, res_mtbseq])
        # Create summary_dict as a dictionary containing categorization per tool per drug and the variants detected by TBProfiler
        # create overview_dict as a list of dictionaries for populating overview table
        summary_dict, overview_dict, prob_dict = get_susceptibility_dicts(df_susceptibility)
        # use summary_dict to get predicted treatment regimens per tool (regimen_dict) and a legend of details per predicted regimen (legend_dict)
        regimen_dict, legend_dict = assign_regimen(summary_dict)
        df_susceptibility.to_excel(os.path.join(file_output_dir, 'susceptibility_results.xlsx'))
        df_susceptibility.to_pickle(os.path.join(file_output_dir, 'dataframe_susceptibility_results.pkl'))
        df_susceptibility.to_csv(os.path.join(file_output_dir, f'susceptibility_results.csv'))
        # Put values of prob_dict in into readable format:
        prob_list = []
        for drug, prob in prob_dict.items():
            prob_list.append(f'{drug.split("_")[1]}: {prob}%')
        prob_string = ', '.join(prob_list)
        summary_dict['prob_string'] = prob_string

        # Create phylogenetics dataframe and get dictionary for mailmerge
        df_phylo = pd.concat([phylo_tbprofiler, phylo_mykrobe, phylo_mtbseq])
        df_phylo.rename(columns={0: 'result'}, inplace=True)
        phylo_dict = get_phylo_dict(df_phylo)
        df_phylo.to_excel(os.path.join(file_output_dir, 'phylogenetic_classification.xlsx'))
        df_phylo.to_pickle(os.path.join(file_output_dir, 'dataframe_phylo.pkl'))

        # Create qc dataframe and get dictionary for mailmerge
        df_qc = pd.concat([qc_tbprofiler, qc_mtbseq])
        df_qc.rename(columns={0: 'result'}, inplace=True)
        qc_dict = get_qc_dict(df_qc)
        # Add some global variables to qc_dict
        qc_dict['sample_id'] = file[0].split('.')[0]
        qc_dict['analysis_time'] = analysis_time
        qc_dict['creation_time'] = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
        qc_dict['mykrobe_v'] = mykrobe_v
        qc_dict['mykrobe_args'] = mykrobe_args
        qc_dict['mykrobe_return_code'] = mykrobe_return_code
        qc_dict['mykrobe_stderr'] = mykrobe_stderr
        qc_dict['tbprofiler_v'] = tbprofiler_v
        qc_dict['tbprofiler_args'] = tbprofiler_args
        qc_dict['tbprofiler_return_code'] = tbprofiler_return_code
        qc_dict['tbprofiler_stderr'] = tbprofiler_stderr
        qc_dict['mtbseq_v'] = mtbseq_v
        qc_dict['mtbseq_args'] = mtbseq_args
        qc_dict['mtbseq_return_code'] = mtbseq_return_code
        qc_dict['mtbseq_stderr'] = mtbseq_stderr
        qc_dict['gentb_v'] = gentb_v
        qc_dict['gentb_args'] = gentb_args
        qc_dict['gentb_return_code'] = gentb_return_code
        qc_dict['gentb_stderr'] = gentb_stderr

        df_qc.to_excel(os.path.join(file_output_dir, 'qc.xlsx'))
        df_qc.to_pickle(os.path.join(file_output_dir, 'dataframe_qc.pkl'))

        if args.phylo == 'true':
            SNP_distances, no_close_samples = get_SNP_distances(output_dir, file[0])
            phylo_dict['SNP_distances'] = str(SNP_distances)
            phylo_dict['no_close_samples'] = str(no_close_samples)
            if no_close_samples > 0:
                cluster_warning = f'Clustering detected with {no_close_samples} reference strains! See SNP distances matrix.'
            else:
                cluster_warning = 'This sample was not recognized as closely related to any of the reference strains.'
            phylo_dict['cluster_warning'] = cluster_warning
        else: print('Comparative phylogenetic analysis skipped')

        docx_mailmerge(summary_dict, overview_dict, phylo_dict, qc_dict, regimen_dict, legend_dict)

    # remove 'raw' folders when light argument is true
    if str(args.light) == 'true':
        for tool in ('tbprofiler', 'mykrobe', 'gentb', 'mtbseq'):
            shutil.rmtree(os.path.join(file_output_dir, f'{tool}_raw'))
    
    print('finished!')

