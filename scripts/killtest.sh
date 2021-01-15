WORKDIR=$1
ITSTOP=$2
TEST_EXEC=$3

CURIT_FILE=$WORKDIR/current_it

while [ ! -f $CURIT_FILE ]
do
	sleep 1
done

READIT=$(tail -n 1 $CURIT_FILE)

while [ $READIT -lt $ITSTOP ]
do
	sleep 1
	READIT=$(tail -n 1 $CURIT_FILE)
done

killall $TEST_EXEC
