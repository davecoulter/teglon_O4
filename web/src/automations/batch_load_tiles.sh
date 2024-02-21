#!/bin/bash

#start=$SECONDS
#
#msg="Processing tiles..."
#echo "${msg}"
#echo "" #newline

files=$(find ./web/events/S190425z/observed_tiles -type f -iname "*.ecsv")
IFS='\n' read -r -a array <<< "$files"

for element in "${array[@]}"
do
    echo "this is $element"
done
#echo "$files"
#for i in "${files[@]}"
#do
#   echo "docker compose run --rm gw_script python3 ./web/src/ingestion/load_observed_tiles.py --gw_id S190425z --healpix_file GW190425_PublicationSamples_flattened.fits.gz --tile_file $i"
#done



#echo "" #newline
#duration=$(( SECONDS - start ))
#msg="**** Batch complete. Processed in: ${duration} [seconds]. ****"
#echo "${msg}"
#echo "" #newline