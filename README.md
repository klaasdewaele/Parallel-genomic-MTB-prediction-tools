# Parallel open-access genomic MTB prediction tools

We provide a Python script for automatized execution of four popular genomic prediction tools for *Mycobacterium tuberculosis* diagnosis:
- Mykrobe: https://github.com/Mykrobe-tools/mykrobe
- TBProfiler: https://github.com/jodyphelan/TBProfiler
- MTBseq: https://github.com/ngs-fzb/MTBseq_source
- GenTB: https://github.com/farhat-lab/gentb-snakemake

Tools are executed serially; output of the tools is collated in three .xlsx files: susceptibility, phylogenetic and QC information.

A WHO treatment regimen is assigned to the predicted susceptibility pattern of each tool using `who_treatment.py` from https://github.com/iqbal-lab-org/tb-amr-benchmarking

Additionally, per sample a clinical report is generated starting from a modifiable template.

## Installation instructions

**Hardware requirements**: MTBseq requires at least 20-25 GB of RAM.
**OS requirements**: tested on CentOS Linux 7 and Debian 4.19.

Install Miniconda: https://docs.conda.io/en/latest/miniconda.html

Clone the `software` directory to your working directory:
- contains `report_template.docx`: a template with MergeFields for merging results
- contains `gentb-snakemake`: cloned from https://github.com/farhat-lab/gentb-snakemake with two adaptations:
   1) Correction of an indentation error in `varMatchUnk.py` line 382
   2) Adaptation of a hardcoded path in `pza_finalpredict_v2_0.RData` to allow it to run in current workflow

Create the `phylo`folder in your working directory: `mkdir path_to_your_working_dir/phylo`

Install respective tools and ETE3 toolkit using provided .yml files (see yml_files`):

For example for Mykrobe:
`conda env create -f yml_files/mykrobe_env.yml`

`docx-mailmerge` and `fasttree=2.1.10` are installed in the base environment.

MTBseq uses GATK3.8, which requires registering this software separately; see: https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md

Download the main.py script to your Conda base environment and make it executable.
Finally, running the main.py script will require adapting the respective paths to your environments, home folder, phylo folder and template.


