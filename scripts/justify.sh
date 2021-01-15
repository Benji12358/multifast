#!/bin/bash

pad=$(printf '%0.1s' "$4"{1..120})
[ ! $4 ] && pad=$(printf '%0.1s' "_"{1..120})



LEFT="$1"
RIGHT="$2"
padlength=$3

printf '%s' "$LEFT"
printf '%*.*s' 0 $((padlength - ${#LEFT} )) "$pad"


printf '%s\n' "$RIGHT"

