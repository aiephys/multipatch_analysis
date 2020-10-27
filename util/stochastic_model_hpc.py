"""
Script for submitting stochastic release model jobs to Moab
"""
import os

qsub_template = """
#!/bin/bash
#PBS -q {queue_name}
#PBS -N {job_id}
#PBS -r n
#PBS -l ncpus=64,mem=1g,walltime=00:30:00
#PBS -o {log_dir}/{job_id}.out
#PBS -j oe
source {conda_dir}/bin/activate 
python {aisynphys_dir}/tools/stochastic_release_model.py --cache-dir={cache_dir} --no-cache --no-gui --workers=64 {expt_id} {pre_cell_id} {post_cell_id}
"""

base_dir = os.getcwd()
conda_dir = base_dir + '/miniconda3'
cache_dir = base_dir + '/cache'
log_dir = base_dir + '/log'
aisynphys_dir = base_dir + '/aisynphys'
for d in [cache_dir, log_dir, aisynphys_dir, conda_dir]:
    assert os.path.isdir(d)

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
    base_dir=base_dir,
    log_dir=log_dir,
    cache_dir=cache_dir,
    conda_dir=conda_dir,
    aisynphys_dir=aisynphys_dir,
    expt_id=expt_id,
    pre_cell_id=pre_cell_id,
    post_cell_id=post_cell_id,
)

qsub_file = '{log_dir}/{job_id}.qsub'
open(qsub_file, 'w').write(qsub)

