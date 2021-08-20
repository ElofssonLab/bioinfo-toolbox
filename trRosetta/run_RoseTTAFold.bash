#!/bin/bash -x

source /opt/miniconda3/etc/profile.d/conda.sh
#source ~/anaconda3/etc/profile.d/conda.sh
conda activate RoseTTAFold
python3 ~/git/RoseTTAFold/network/predict_complex.py -m /scratch3/arnee/RoseTTAFold/weights -Ls $1 $2 -i $3 -o $4
