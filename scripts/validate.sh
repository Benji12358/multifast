#!/bin/bash

check_results(){

	if [ -d "$REF_PATH" ]; then
		AA_FILE=$DNS_PATH/$DNS_NAME/Log/sensor.csv
		REF_AAFILE=$REF_PATH/$DNS_NAME/Log/sensor.csv

		if [[ ! $IT_TOTEST ]]; then IT_TOTEST=50; fi
		IT_TOTEST=$((IT_TOTEST+1))

		$CHECKSCRIPT $AA_FILE $REF_AAFILE 1e-14 $IT_TOTEST
	else
		error "ERROR: REF simulation doesnt found !!!"
		exit 1;

fi
}


DNS_PATH=$1
DNS_NAME=$2
REF_PATH=$3
IT_TOTEST=$4



SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CHECKSCRIPT=$SCRIPT_DIR/check.sh
TEXTHANDLER_SCRIPT=$SCRIPT_DIR/texthandler.sh


source $TEXTHANDLER_SCRIPT

check_results $DNS_PATH $DNS_NAME $REF_PATH $IT_TOTEST

