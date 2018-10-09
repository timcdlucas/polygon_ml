#!/bin/bash

###
# If you need to rerun a specific fold that failed
###

if [ $# -ne 3 ]; then
  echo $0: usage: $0 cv_type data_type fold
  exit 1
fi

if [[ "$1" != @(random|spatial) ]]; then
  echo "cv_type must be random or spatial"
  exit 1
fi

if [[ "$2" != @(covs|ml|all) ]]; then
  echo "data_type must be covs, ml or all"
  exit 1
fi

Rscript run_single_cv_fold.R "$3" "$1" "$2" & 


echo 'finished'
