#!/bin/bash
#rm -r TMP/TEST_OPEN_STREAMWISE_1/
rm -r TMP/TEST_Laminar/
make clean_all
make
./job_TEST_laptop.sh
