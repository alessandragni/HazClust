#!/bin/bash

seeds=($(seq 0 99))
gammapars=(0 1e-6 1e-4 1e-3 1e-2 1e-1 0.2 0.4)
ks=(20 50)
cs=(3)
censoring=("administrative" "normal")
scenarios=("baseline")
type=("tab2")
frailty_intensities=(0.5)

mkdir -p logs/tab2

for seed in "${seeds[@]}"; do
  for gammapar in "${gammapars[@]}"; do
    for k in "${ks[@]}"; do
      for c in "${cs[@]}"; do
        for censor in "${censoring[@]}"; do
          for scenario in "${scenarios[@]}"; do
            for frailty_intensity in "${frailty_intensities[@]}"; do
            
              echo "Submitting: seed=$seed, gamma=$gammapar, k=$k, c=$c, censor=$censor, scenario=$scenario, type=$type, frailty_intensity=$frailty_intensity"
            
              output_file="logs/tab2/out_seed${seed}_g${gammapar}_k${k}_c${c}_${censor}_${scenario}_${frailty_intensity}.txt"
              error_file="logs/tab2/err_seed${seed}_g${gammapar}_k${k}_c${c}_${censor}_${scenario}_${frailty_intensity}.txt"
            
            qsub -v seed=$seed,gammapar=$gammapar,k=$k,c=$c,censor=$censor,scenario=$scenario,type=$type,frailty_intensity=$frailty_intensity \
                 -o "$output_file" \
                 -e "$error_file" \
                 enhanced_simulation.job &
            
            done
          done
        done
      done
    done
  done
done

echo "All jobs submitted!"

