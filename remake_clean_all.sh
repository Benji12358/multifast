#!/bin/bash
make clean_all
make
./job_laptop_singleproc.sh $1 $2 $3 $4
