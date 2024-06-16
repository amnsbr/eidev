import os
import sys
import numpy as np
import pandas as pd

PROJECT_DIR = os.environ.get('PROJECT_DIR')

def gen_dag(dataset, optima_csv):
    optima = pd.read_csv(optima_csv, index_col=0)
    best_runs = optima.groupby('sub')['gof'].idxmax().values
    data = optima.loc[best_runs].set_index('sub')
    dag_filename = f'perturbation_{dataset}.dag'
    f = open(os.path.join(PROJECT_DIR, 'scripts', 'modeling', 'run', 
                          'cuBNM', 'perturbation', dag_filename), 'w')
    submit_path = os.path.join(PROJECT_DIR, 'scripts', 'modeling', 'run', 
                          'cuBNM', 'perturbation', 'run_perturbation.submit')
    # randomly select 40 subjects
    data = data.sample(n=40, random_state=0)
    job_count = 0
    for sub, row in data.iterrows():
        f.write(f'JOB job_{job_count} {submit_path}\n')
        f.write(f'VARS job_{job_count} dataset="{dataset}" sub="{sub}" SeedMW="{row["SeedMW"]}"\n\n')
        job_count += 1
    f.write(f'\nCATEGORY ALL_NODES {dataset}')
    f.write(f'\nMAXJOBS {dataset} 5\n\n') # occupy 5 GPU nodes at a time
    f.close()
    


if __name__=='__main__':
    dataset = sys.argv[1]
    optima_csv = sys.argv[2]
    gen_dag(dataset, optima_csv)