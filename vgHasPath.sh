#!/bin/bash
# test if vg file contains a path with given name

set -e

VG_PATH=$1
NAME=$2

NUM_PATHS=`vg view $VG_PATH -j | jq .path | jq length`

for (( i=0; i<$NUM_PATHS; i++ ))
do
	 PATH_NAME=`vg view $VG_PATH -j | jq .path[${i}] | jq .name`
	 if [ "$PATH_NAME" == "\"$NAME\"" ]; then
		  echo 1
		  exit 0
	 fi
done

echo 0

