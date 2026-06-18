#!/bin/bash

seeds=($(seq 0 99))
cs=(3)
censoring=("administrative")
type=("tabS62")
frailty_intensities=(0.5 1 1.5)

# Scenarios
scenarios=(
  "baseline" "misspec_frailty" "misspec_baseline" "misspec_both"
  "weak_separation" "medium_separation" "imbalanced_moderate"
  "imbalanced_severe" "worst_case"
)

mkdir -p logs/tabS62

for seed in "${seeds[@]}"; do
  for c in "${cs[@]}"; do
    for censor in "${censoring[@]}"; do
      for scenario in "${scenarios[@]}"; do
        for frailty_intensity in "${frailty_intensities[@]}"; do

          echo "Submitting: seed=$seed, c=$c, censor=$censor, scenario=$scenario, type=$type, frailty_intensity=$frailty_intensity"

          output_file="logs/tabS62/out_seed${seed}_c${c}_${censor}_${scenario}_${frailty_intensity}.txt"
          error_file="logs/tabS62/err_seed${seed}_c${c}_${censor}_${scenario}_${frailty_intensity}.txt"

          qsub -v seed=$seed,c=$c,censor=$censor,scenario=$scenario,type=$type,frailty_intensity=$frailty_intensity \
               -o "$output_file" \
               -e "$error_file" \
               alternative_simulation.job &

        done
      done
    done
  done
done

echo "All jobs submitted!"