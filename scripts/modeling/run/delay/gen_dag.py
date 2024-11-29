import os
import sys
import numpy as np
import pandas as pd

PROJECT_DIR = os.environ['PROJECT_DIR']
INPUT_DIR = os.environ['INPUT_DIR']


def gen_dag():
    optima_csv = os.path.join(PROJECT_DIR, 'results', 'pnc_optima.csv')
    if not os.path.exists(optima_csv):
        print("Run Figure 2 first")
        return
    optima = pd.read_csv(optima_csv, index_col=0)
    # limit the optima to subsample
    subsample = np.loadtxt(os.path.join(INPUT_DIR, 'pnc_subsample_200.txt'), dtype=str)
    optima = optima.loc[optima['sub'].isin(subsample)]
    # select best runs
    best_runs = optima.groupby('sub')['gof'].idxmax().values
    data = optima.loc[best_runs].set_index('sub')
    dag_filename = 'delay.dag'
    f = open(os.path.join(PROJECT_DIR, 'scripts', 'modeling', 'run', 
                          'cuBNM', 'delay', dag_filename), 'w')
    submit_path = os.path.join(PROJECT_DIR, 'scripts', 'modeling', 'run', 
                          'cuBNM', 'delay', 'run_delay.submit')
    job_count = 0
    for sub, row in data.iterrows():
        f.write(f'JOB job_{job_count} {submit_path}\n')
        f.write(f'VARS job_{job_count} sub="{sub}" SeedMW="{row["SeedMW"]}"\n\n')
        job_count += 1
    f.write('\nCATEGORY ALL_NODES delay')
    f.write('\nMAXJOBS delay 5\n\n') # occupy 5 GPU nodes at a time
    f.close()
    


if __name__=='__main__':
    gen_dag()