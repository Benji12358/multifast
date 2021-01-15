#!/bin/bash

IT1=101000
EVERY=1000
IT2=200000

for numero in `seq $IT1 $EVERY $IT2`
    do rm ./field"$numero"/Elec* ./field"$numero"/void* ./field"$numero"/sca* ./field"$numero"/Sca*
done
