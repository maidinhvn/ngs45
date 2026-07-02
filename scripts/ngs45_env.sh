#!/bin/sh
# Convenience PATH for running ngs45 against tools already installed on this host
# WITHOUT building the dedicated `ngs45` conda env. For real use, prefer:
#   conda env create -f environment.yml && conda activate ngs45
# This composite PATH is only for local development/testing on this machine.
E=/path/to/.conda/envs
export PATH="$E/easy45/bin:$E/getorganelle_env/bin:$E/bcftools_env/bin:/usr/bin:$PATH"
export PYTHONPATH="/path/to/ngs245/src:$PYTHONPATH"
# use the biopython-equipped interpreter from easy45 for ngs45 itself
alias ngs45="$E/easy45/bin/python -m ngs45.cli"
