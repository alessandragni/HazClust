#!/bin/bash

seeds=($(seq 0 99))
gammapars=(1e-4 1e-1)
ks=(20) #50
cs=(3)
censoring=("administrative") #start from this only
type=("tabS61")
frailty_intensities=(0.5 1 1.5)

# Scenarios
scenarios=(
  "baseline" "misspec_frailty" "misspec_baseline" "misspec_both"
  "weak_separation" "medium_separation" "imbalanced_moderate"
  "imbalanced_severe" "worst_case"
)



mkdir -p logs/tabS61

for seed in "${seeds[@]}"; do
  for gammapar in "${gammapars[@]}"; do
    for k in "${ks[@]}"; do
      for c in "${cs[@]}"; do
        for censor in "${censoring[@]}"; do
          for scenario in "${scenarios[@]}"; do
            for frailty_intensity in "${frailty_intensities[@]}"; do
            
              echo "Submitting: seed=$seed, gamma=$gammapar, k=$k, c=$c, censor=$censor, scenario=$scenario, type=$typee, frailty_intensity=$frailty_intensity"
            
              output_file="logs/tabS61/out_seed${seed}_g${gammapar}_k${k}_c${c}_${censor}_${scenario}_${frailty_intensity}.txt"
              error_file="logs/tabS61/err_seed${seed}_g${gammapar}_k${k}_c${c}_${censor}_${scenario}_${frailty_intensity}.txt"
            
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