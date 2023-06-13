#!/bin/bash

start=$SECONDS

START=51
END=100

GW_ID='S190425z'
HEALPIX_FILE='GW190425_PublicationSamples_flattened.fits.gz'
NUM_CPU=4
MODEL_TYPE=linear
BATCH_DIR=offaxis

#manual_dirs=(10 17 20 27 30 37 40 47 50 57 60 67 70 73)

msg="Processing dirs [${START},${END}]"
echo "${msg}"

for i in $( seq $START $END ); do
#for i in ${manual_dirs[@]}; do

  msg="Starting job #${i}"
  echo "${msg}"

  docker compose run --rm gw_script python ./web/src/analysis/model_detection_efficiency.py \
  --gw_id $GW_ID \
  --healpix_file $HEALPIX_FILE \
  --sub_dir $i \
  --num_cpu $NUM_CPU \
  --model_type $MODEL_TYPE \
  --batch_dir $BATCH_DIR

  msg="Job #${i} complete"
  echo "${msg}"

done

duration=$(( SECONDS - start ))
msg="**** Batch complete. Directories [${START},${END}] processed in: ${duration} [seconds]. ****"
echo "${msg}"
echo "" #newline

