# Adolescent development of E-I using biophysical network modeling
This repository includes the code associated with the paper "Adolescent maturation of cortical excitation-inhibition balance based on individualized biophysical network modeling", Saberi et al.

## Structure
- `scripts`: Scripts used for image preprocessing and running biophysical network modeling simulation-optimization jobs. Read the details in [`./scripts/`](scripts/).
- `tools`: Software developed in this project for biophysical network modeling simulation-optimization on GPU and CPU
    - [`bnm_cuda`](https://github.com/amnsbr/bnm_cuda): Main C++/CUDA used for model simulation-optimization
    - [`cuBNM`](https://github.com/amnsbr/cuBNM): A Python toolbox containing a similar C++/CUDA core that was used in some of the analyses
- `results`: Jupyter notebooks and helper scripts used to run statistical tests and generate the paper main and supplementary figures. Note that these notebooks cannot be run without the input data which is restricted-access and cannot be shared publicly. However, some of the data such as the maps of statistical effects are shared within this folder.

## Dependencies
- Data: The scripts require input data of [Philadelphia Neurodevelopmental Cohort](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000607.v3.p2) and the [IMAGEN dataset](https://imagen-project.org/the-imagen-dataset/) to exist in pre-specified directories defined under `$PROJECT_DIR`, `$PNC_PROJECT_DIR`, `$IMAGEN_PROJECT_DIR`, `$INPUT_DIR` and `$OUTPUT_DIR` which are expected to be defined as environment variables. However, data used in this project is restricted-access and cannot be shared publicly. Accordingly we have not shared the scripts used for downloading and handling the input data.
- Python 3 (tested with version 3.10) and packages listed in [`./scripts/setup/requirements.txt`](scripts/setup/requirements.txt) are needed to run the Python scripts.
- Singularity (tested with version 3.11.4) is needed to build and run the following containerized software:
    - [<img src="https://img.shields.io/badge/docker-freesurfer/freesurfer:7.1.1-blue.svg?logo=docker">](https://hub.docker.com/layers/freesurfer/freesurfer/7.1.1/images/sha256-922fc8242c6dac65529e11e93c251b946ed6772f734a8ed589695f05df8a37e1?context=explore)
    - [<img src="https://img.shields.io/badge/docker-nipreps/fmriprep:22.0.0-blue.svg?logo=docker">](https://hub.docker.com/layers/nipreps/fmriprep/22.0.0/images/sha256-30bbf9ebc870de988b9ed7a9e32414dfdd98d57794f3b30126e1156079983402?context=explore)
    - [<img src="https://img.shields.io/badge/docker-micalab/micapipe:v0.1.2-blue.svg?logo=docker">](https://hub.docker.com/layers/micalab/micapipe/v0.1.2/images/sha256-b2f94d47c9810105cf020fc0deda6c1b51eeddd1444593963612ba5c7d9a3cfd?context=explore): Note that we used an earlier version (0.1.1) which is currently not available on Docker Hub.
- `scripts/modeling/run` scripts require Nvidia GPUs and a few additional dependencies:
    - `scripts/modeling/run/bnm_cuda`: See [here](https://github.com/amnsbr/bnm_cuda?tab=readme-ov-file#build-dependencies)
    - `scripts/modeling/run/cuBNM`: Additionally requires `cuBNM` to be installed in the Python virtual environment (not installed via `./scripts/setup/requirements.txt`). For `cuBNM` requirements see [here](https://github.com/amnsbr/cuBNM/tree/main?tab=readme-ov-file#installation)
- The preprocessing and modeling scripts were run as jobs on our HTC and HPC clusters which use HTCondor and Slurm for job scheduling, respectively.

## Support
Feel free to contact me (amnsbr\[at\]gmail.com, a.saberi\[at\]fz-juelich.de, saberi\[at\]cbs.mpg.de) if you have any questions.