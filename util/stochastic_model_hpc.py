"""
Script for submitting stochastic release model jobs to Moab
"""
import os, subprocess, re
from aisynphys.database import default_db as db
import aisynphys.config


def optint(x):
    return None if x in (None, '--') else int(x)

def qstat():
    stat = []
    lines = subprocess.check_output(['qstat', '-n']).split(b'\n')[5:]
    while len(lines) > 1:
        l1 = lines.pop(0).decode()
        l2 = lines.pop(0).strip().decode()
        cols = re.split('\s+', l1)

        if l2 != '--':
            m = re.match(r'(n\d+)/(\d+)(-(\d+))?', l2)
            assert m is not None, l2.decode()
            node, start_cpu, _, end_cpu = m.groups()
            start_cpu = int(start_cpu)
            end_cpu = optint(end_cpu)
        else:
            node = None
            start_cpu = None
            end_cpu = None

        stat.append({
            'job_id': cols[0].split('.')[0], 
            'user': cols[1],
            'queue': cols[2],
            'job_name': cols[3],
            'nodes': optint(cols[5]),
            'tasks': optint(cols[6]),
            'mem_req': cols[7],
            'time_req': cols[8],
            'status': cols[9],
            'elapsed_time': cols[10],
            'node': node,
            'start_cpu': start_cpu,
            'end_cpu': end_cpu,
        })

    return stat


qsub_template = """#!/bin/bash
#PBS -q {queue_name}
#PBS -N {job_id}
## Job is rerunable
#PBS -r n
#PBS -l ncpus={ncpus:d},mem={mem},walltime={walltime}
#PBS -o {log_path}/{job_id}.out
#PBS -j oe
source {conda_path}/bin/activate {conda_env}
python {aisynphys_path}/../tools/stochastic_release_model.py --cache-path={cache_path} --no-gui --workers=$PBS_NP {expt_id} {pre_cell_id} {post_cell_id} > {log_path}/{job_id}.log 2>&1
"""

# we are going to distribute a variety of different limit options in order to assist the HCP scheduler
limit_pool = [
    {'ncpus': 88, 'mem': '32g', 'walltime': '2:00:00'},
    {'ncpus': 88, 'mem': '32g', 'walltime': '2:00:00'},
    {'ncpus': 32, 'mem': '32g', 'walltime': '8:00:00'},
    {'ncpus': 32, 'mem': '32g', 'walltime': '8:00:00'},
    {'ncpus': 24, 'mem': '32g', 'walltime': '12:00:00'},
]


base_path = os.getcwd()
conda_path = base_path + '/miniconda3'
log_path = base_path + '/log'

aisynphys_path = os.path.split(aisynphys.__file__)[0]
cache_path = os.path.join(aisynphys.config.cache_path, 'stochastic_model_results')

for d in [cache_path, log_path, aisynphys_path, conda_path]:
    assert os.path.isdir(d), f'Missing path: {d}'

job_id = 'synphys_model_test'
#queue_name = 'aibs_dbg'
queue_name = 'celltypes'

pairs = db.pair_query(synapse=True).all()

queued_job_ids = {job['job_id'] for job in qstat()}

for i,pair in enumerate(pairs):
    expt_id = pair.experiment.ext_id
    pre_cell_id = pair.pre_cell.ext_id
    post_cell_id = pair.post_cell.ext_id

    job_id = f'{expt_id}_{pre_cell_id}_{post_cell_id}'

    if job_id in queued_job_ids or os.path.exists(f'{cache_path}/{job_id}.pkl'):
        print(f'{job_id} => SKIP')
        continue

    format_opts = dict(
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
    limit_opts = limit_pool[i % len(limit_pool)]
    format_opts.update(limit_opts)
    qsub = qsub_template.format(**format_opts)

    qsub_file = f'{log_path}/{job_id}.qsub'
    open(qsub_file, 'w').write(qsub)

    sub_id = subprocess.check_output(f'qsub {qsub_file}', shell=True).decode().strip()
    open(qsub_file, 'a').write(f"\n# {sub_id}\n")

    print(f"{job_id} => {sub_id} {limit_opts}")
