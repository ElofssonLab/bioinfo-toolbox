#!/bin/bash

NAME=$1
EVAL=$2
N=3

cd loop_jackhmmer
tar -vzxf output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.tar.gz
cd ..
#mkdir ./alignments_dis/${NAME}
#mkdir ./input/${NAME}
#cp ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.a2m ./alignments_dis/${NAME}/${NAME}_e${EVAL}_n${N}_m50_f0.a2m
#cp ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.a3m ./alignments_dis/${NAME}/${NAME}_e${EVAL}_n${N}_m50_f0.a3m
#cp ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}_MI_DIs.txt ./alignments_dis/${NAME}/${NAME}_e${EVAL}_n${N}_m50_f0_MI_DIs.txt
#cp ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.fa ./input/${NAME}/${NAME}.fa
#cp ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.fa ./input/${NAME}/${NAME}_memsatsvm3state.txt

mv ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.a3m ./alignments/jackhammer/

rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.a2m
#rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.a3m
rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.hamming_input
rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.hhr
rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.matlablog
rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}_MI_DIs.txt
rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.stdout
rm ./loop_jackhmmer/output/${NAME}_n${N}_m50_f0_t0.3_g_r_e${EVAL}.weight_table

