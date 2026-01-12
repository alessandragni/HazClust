#!/bin/bash

seeds=($(seq 0 99))
gammapars=(1e-8 1e-6 1e-4 1e-3 1e-2 1e-1 0.2 0.4)
frailties=("gamma")
ks=(20 50) # 5 20 50
cs=(3)
censoring=("administrative" "normal") # "administrative" "exponential" "normal" "uniform"

mkdir -p logs

for seed in "${seeds[@]}"; do
  for gammapar in "${gammapars[@]}"; do
    for frailty in "${frailties[@]}"; do
      for censor in "${censoring[@]}"; do
        for k in "${ks[@]}"; do
          for c in "${cs[@]}"; do
            echo "Submitting job with seed=$seed, gammapar=$gammapar, frailty=$frailty, censor=$censor, k=$k, c=$c"
            output_file="logs/output_seed${seed}_gammapar${gammapar}_frailty${frailty}_censor${censor}_k${k}_c${c}.txt"
            error_file="logs/error_seed${seed}_gammapar${gammapar}_frailty${frailty}_censor${censor}_k${k}_c${c}.txt"
            qsub -v seed=$seed,gammapar=$gammapar,frailty=$frailty,k=$k,c=$c,censor=$censor \
                 -o "$output_file" \
                 -e "$error_file" \
                 sim1COV.job &
          done
        done
      done
    done
  done
done



