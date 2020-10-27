"""
Script for submitting stochastic release model jobs to Moab
"""
import os

qsub_template = """#!/bin/bash
#PBS -q {queue_name}
#PBS -N {job_id}
#PBS -r n
#PBS -l ncpus=64,mem=2g,walltime=00:30:00
#PBS -o {log_path}/{job_id}.out
#PBS -j oe
source {conda_path}/bin/activate {conda_env}
python {aisynphys_path}/tools/stochastic_release_model.py --cache-path={cache_path} --no-cache --no-gui --workers=64 {expt_id} {pre_cell_id} {post_cell_id}
"""

base_path = os.getcwd()
conda_path = base_path + '/miniconda3'
cache_path = base_path + '/cache'
log_path = base_path + '/log'
aisynphys_path = base_path + '/aisynphys'
for d in [cache_path, log_path, aisynphys_path, conda_path]:
    assert os.path.isdir(d), f'Missing path: {d}'

job_id = 'synphys_model_test'
queue_name = 'aibs_dbg'
#queue_name = 'celltypes'

expt_id = '1530559621.966'
pre_cell_id = '7'
post_cell_id = '6'
job_id = f'{expt_id}_{pre_cell_id}_{post_cell_id}'

qsub = qsub_template.format(
    queue_name=queue_name,
    job_id=job_id,
    base_path=base_path,
    log_path=log_path,
    cache_path=cache_path,
    conda_path=conda_path,
    aisynphys_path=aisynphys_path,
    conda_env='aisynphys',
    expt_id=expt_id,
    pre_cell_id=pre_cell_id,
    post_cell_id=post_cell_id,
)

qsub_file = f'{log_path}/{job_id}.qsub'
open(qsub_file, 'w').write(qsub)

