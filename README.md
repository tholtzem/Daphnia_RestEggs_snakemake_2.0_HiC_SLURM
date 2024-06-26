# Daphnia_RestEggs_snakemake_2.0_HiC_SLURM

This is the continued version of the snakemake workflow [Daphnia_RestEggs_snakemake_pbs_2.0_HiC](https://github.com/tholtzem/Daphnia_RestEggs_snakemake_pbs_2.0_HiC). The previous version was designed for the execution of the mach2 HPC cluster (PBS) which was discontinued. The rest of the analysis was performed on the leo5 HPC cluster with a SLURM batch job sumission system which is shown in the worklfow below. 


======================================================

## Some softwares required manual installation using a singularity image

* angsd (angsd+ngsrelate.sif)
* ngsRelate (angsd+ngsrelate.sif)
* ngsLD (available in ngs+tools.sif)
* pcangsd (available in ngs+tools.sif)

To copy image to your repository

```
rsync -avP /scratch/c770xxxx/bio/angsd+ngsrelate.sif $SCRATCH/
rsync -avP /scratch/c770xxxx/bio/ngs+tools.sif $SCRATCH/

```

To execute image in snakemake rule

```
module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
#angsd v0.939 (used here)
singularity exec --home $PWD:$HOME $SCRATCH/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd
#angsd v0.940
singularity exec --home $PWD:$HOME $SCRATCH/bio/angsd+ngsrelate.sif /opt/angsd/angsd
singularity exec --home $PWD:$HOME $SCRATCH/bio/ngs+tools.sif pcangsd
singularity exec --home $PWD:$HOME $SCRATCH/bio/ngs+tools.sif /opt/ngsTools/ngsLD/ngsLD
```

## Anaconda and Mamba


Mamba (https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++.

```
# Load Anaconda on cluster (here leo5):
module load Anaconda3/2023.10/miniconda-base-2023.10
 

# To use conda commands in your current shell session, first do:
eval "$(/$UIBK_CONDA_DIR/bin/conda shell.bash hook)"

# Users of the c770 group can simply activate the environment on leo5 with:
conda activate daphnia

# !!! In case you want to change the environment or create your own !!! 
## Create environment from yaml file (in envs/):
mamba env create -f envs/s21.yaml -p $SCRATCH/envs/daphnia
conda config --append envs_dirs $SCRATCH/envs

# Software can be installed/updated by modifying the envs/s21.yaml file.
# To use the modified conda environment, update with:
mamba env update --name daphnia --file envs/s21.yaml

```

## Users of the c770 group can get a simple version of the snakemake pipeline for the slurm based cluster (here leo5)

```
mkdir <some_directory>
cd <some_directory>
rsync -avP /scratch/c770xxxx/snakemake ./

```

## Some information on the snakemake pipeline for the slurm based cluster (here leo5)

1. Testing and execution of snakemake in working directory (where your snakefile is located)

2. Slurm specific files and logs are in slurm/ directory
* cluster submission script for main job (slurm/clusterSnakemake.sh) **Do not forget to change the job-name and mail-user!**
* configuration profile for slurm (slurm/config.yaml)
* sbatch log files (output & error files) are sent to slurm/log/

3. The snakefile:
* contains the names of the output files (should correspond to the output in your rule)
* always includes the file rules/common.smk
* includes additional files with rules you want to execute (f.ex. rules/test.smk)

4. The rules/common.smk file:
* calls the input config file (config/config.yaml) and 
* contains python code for loading sample information (and helper functions)

5. The snakemake configuration file (config/config.yaml) contains the paths to adapters, Kraken2 database and references 

6. Sample information and metadata are in list/samples.tsv

7. Data was downloaded from https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz **NOT required for c770 group!**

## Local execution of snakemake (for testing only, else not recommended!)

```
# dry-run
snakemake -n

```

## Run snakemake via cluster execution on leo5

For more information on leo5 and slurm (https://www.uibk.ac.at/zid/systeme/hpc-systeme/common/tutorials/slurm-tutorial.html)

```
sbatch slurm/clusterSnakemake.sh

```

## Check/cancel jobs

```
# check all jobs
sq

# check your jobs only
squ

# cancel job
scancel <JOBID>

# check usage of running job
## ssh to node (NODELIST)
ssh <NODELIST>
## on the compute node
top

# check usage after the job has completed
export SACCT_FORMAT="JobID%20,JobName,User,Partition,NodeList,Elapsed,State,ExitCode,MaxRSS,AllocTRES%32"
sacct -j <JOBID>

```
