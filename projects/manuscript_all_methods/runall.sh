#!/usr/bin/env bash

#activate conda environment
source activate netpert_env

#project directory
project_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#run NetPert, BC, and TieDIE
python ./bin/netpert.py all_methods -s mouse $project_dir Twist1

