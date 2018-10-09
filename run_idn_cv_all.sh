#!/bin/bash

declare -a cv_type=("random" "spatial")
declare -a model_type=("covs" "ml" "all")

for i in "${cv_type[@]}" 
do
   for j in "${data_type[@]}" 
   do
      echo "Run $i $j"
      ./run_idn_cv_single.sh $i $j &
   done
done

