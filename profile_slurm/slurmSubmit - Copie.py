#!/usr/bin/env python3

"""
Submit Snakemake jobs to SLURM.

Usage:
snakemake -j 99 --use-conda --cluster-config cluster.json --immediate-submit --notemp --cluster 'python3 slurmSubmit.py {dependencies}'
"""

import os
import sys
import subprocess
from snakemake.utils import read_job_properties

# Read the job script and properties
jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

# Construct the SLURM command
cmdline = ["sbatch"]
for param, val in job_properties.get('cluster', {}).items():
    cmdline.extend([f"--{param}", str(val)])

# Add dependencies if they exist
dependencies = sys.argv[1:-1]
if dependencies:
    cmdline.extend(["--dependency=afterok:{}".format(":".join(dependencies))])

# Add the job script
cmdline.append(jobscript)

# Submit the job and capture the job ID
try:
    result = subprocess.run(
        cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True
    )
    # Extract job ID from SLURM output
    job_id = result.stdout.strip().split()[-1]
    sys.stdout.write(job_id + "\n")
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Error submitting job: {e.stderr}\n")
    sys.exit(1)
