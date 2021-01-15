#!/bin/bash
hcenter() {

text="$2"

cols="$3"

IFS=$'\n'$'\r'
for line in $(echo -e $text); do

	line_length=`echo $line| wc -c`
	half_of_line_length=`expr $line_length / 2`
	center=`expr \( $cols / 2 \) - $half_of_line_length`

	spaces=""

	for i in `seq 1 $center`; do
		spaces="$spaces "
	done

	echo $1"$spaces$line$spaces"$1

	done

}


title(){
	title="$1"
	width=36
	echo
	hcenter "================================" "=======================" $width;
	hcenter "================================" "$title" $width
	hcenter "================================" "=======================" $width;
	echo
}

section(){
	section="$1"
	width=36
	hcenter "--------------------------------" "$section" $width
}

warning(){
	warnmess="$1"
	orange='\033[0;33m'
	NC='\033[0m'

	echo -e ${orange}"$warnmess" ${NC}
	
}

error(){
	errmess="$1"
	red='\033[0;31m'
	NC='\033[0m'

	echo -e ${red}"$errmess" ${NC}
	
}

success(){
	successmess="$1"
	green='\033[0;32m'
	NC='\033[0m'

	echo -e ${green}"$successmess" ${NC}
	
}

if [[ "$(basename -- "$0")" == "texthandler.sh" ]]; then
	title "$1"
fi
