# NetPert v1
Ranks genes and their proteins in network pathways between a driver gene and downstream response genes using perturbation theory and response functions.

## Clone from GitHub

Clone the branch `manuscript` from the GitHub repository.

```shell
git clone -b manuscript https://github.com/joelbaderlab/netpert_v1.git
```

Change directory

```shell
cd netpert_v1
```

Create your own branch

```
git checkout -b manuscript_username
```

## File and Directory Organization

`bin`: scripts.

`databases`: fixed data files downloaded from external sources and created by NetPert.

`projects`: project directories containing data and results specific to a project. An example directory is provided:
`projects/manuscript_results` - an example directory for outputting manuscript results.

## Setup: download large files (~1.3 GB size) and create conda environment

Run setup bash script
```shell
bash bin/setup.sh
```

## Activate conda environment

```shell
conda activate netpert_env
```

## Create a new project directory

Make a new project directory for your work

```shell
mkdir projects/dir_name 
```

## Run NetPert analysis

```shell
python ./bin/netpert.py analysis [-h] [-s {human,mouse}] project_dir driver
```
positional arguments:
project_dir	project directory name.
driver	network driver gene symbol.

optional arguments:
-h, --help	show this help message and exit.
-s {human,mouse}    species. Default human.

Required files in project directory:
1. responses file

Create a response file:

Make a file called `responses` in the project directory. The file must have extension .xlsx, .csv, or .tsv. Each row should contain the gene symbol of one response gene. There is to be no header. See `projects/manuscript_results/responses.csv` for an example.  

Output:
NetPert outputs a file called `geneRankings.tsv` to the project directory. The file contains the ranking of each gene in the network. The columns are as follows:
`type` D: driver, R: response, DIR: intermediate having direct connections to the driver and at least one response gene, DI: intermediate connected directly to the driver but not to a response gene, IR: intermediate connected directly to a response gene but not the driver, and I: intermediate directly connected to neither the driver nor a response gene.
`DIIR` DIIR: a subset of DI and IR genes that are on a path of length 3 between driver and a response.
`NetPert_rank` NetPert ranking.
`NetPert_weight` NetPert weight.
`repurposinghub_cpds` Drug Repurposing Hub compounds that target the protein of the specified gene. Compounds with more than 5 targets were omitted. 

Example: output manuscript results

```shell
python ./bin/netpert.py analysis -s mouse ./projects/manuscript_results Twist1
```

## Create a driver script (recommended)

It is recommended to create a driver script (not to be confused with driver gene) in the project directory that carries out every command. See `projects/manuscript_results/runall.sh` for an example of running the NetPert analysis.

Run driver script example:
```shell
bash projects/manuscript_results/runall.sh
```
