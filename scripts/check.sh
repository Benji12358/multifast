
#!/bin/bash

AA_FILE=$1
REF_AAFILE=$2
tolerance=$3
IT_TOTEST=$4


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TEXTHANDLER_SCRIPT=$SCRIPT_DIR/texthandler.sh
CMPSCRIPT=$SCRIPT_DIR/cmp.py


source $TEXTHANDLER_SCRIPT

line=$IT_TOTEST,"$IT_TOTEST"p
ref_lastline=`sed -n $line $REF_AAFILE`
lastline=`sed -n $line $AA_FILE`

IFS=',' read -a VALS <<< "$lastline"
IFS=',' read -a REFVALS <<< "$ref_lastline"


SUCCES=1

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

diff_str=""
str=""

for i in ${!REFVALS[*]}; do 
	A=${REFVALS[$i]}
	B=${VALS[$i]}
	python $CMPSCRIPT $A $B $tolerance
	res=$?


	if (( res==1)); then
		SUCCES=0
		diff_str=$diff_str"${red}${VALS[$i]}${reset}"
	else
		diff_str=$diff_str"${green}${VALS[$i]}${reset}"
	fi

	
	str=$str"${reset}${REFVALS[$i]}"

done

if (( SUCCES==1)); then
	success "---------SUCCES----------"
else
	error "---------FAILURE----------"
fi

echo -e $str
echo -e $diff_str
echo


if (( SUCCES==1)); then
	exit 0
else
	exit 1
fi
