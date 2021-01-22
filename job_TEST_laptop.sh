#!/bin/bash

launch_simulation(){
	cd $TMPDIR
	
	echo $TEXT_LEVEL1 Running simulation...
	mpirun -np $NPROCS -x LD_LIBRARY_PATH $TMPDIR/DNS_EXEC $DNS_NAME $DNS_TIME $PROW $PCOL 0 > $DNS_OUTPUT 2>temp_arrays.log

	cd -
}

list_dir_to_create(){

	cd $TMPDIR
	
	#./PREPROCESS $DNS_NAME $dir_to_create_list >lol
	seq $IT1 $EVERY $IT2 >	$dir_to_create_list
	echo $TEXT_LEVEL2 The following results directory will be created:
	more $dir_to_create_list

	cd -
}

create_results_dir()
{
	cd $TMPDIR/$DNS_NAME/Results/3D	

	for numero in `seq $IT1 $EVERY $IT2`
    		do mkdir ./field"$numero"
	done

	cd -
}

build_running_env(){

	echo $TEXT_LEVEL1 Building running environnment...

	rm -rf $TMPDIR/*
	mkdir $TMPDIR
	cp -r $SIMULATION_DIR/.recovery/arborescence $TMPDIR/$DNS_NAME

	cp -r $SIMULATION_DIR/Input $TMPDIR/$DNS_NAME
	cp -r $SIMULATION_DIR/Log $TMPDIR/$DNS_NAME
	cp -r $SIMULATION_DIR/.recovery $TMPDIR/$DNS_NAME

	cp $DNS_CODE/DNS_EXEC $TMPDIR
	#cp $DNS_CODE/PREPROCESS $TMPDIR

}

clean_results()
{

	for numero in `more $dir_to_create_list`
    		do field=$TMPDIR/$DNS_NAME/Results/3D/field$numero
		nb_file=`ls -l $field|wc -l`
		echo $TMPDIR/$DNS_NAME/Results/3D/field$numero contient $nb_file fichiers

		if [ $nb_file -ne 7 ]
		then
		   echo $field "est vide";
		   rm -rf $field
		fi
	done

}

save_simulation(){

	echo $TEXT_LEVEL1 Saving results...

	cp $TMPDIR/core* $TMPDIR/$DNS_NAME

	clean_results


	cp $TMPDIR/$DNS_OUTPUT $SIMULATION_DIR/
	cp -r $TMPDIR/$DNS_NAME/Results/3D/* $SIMULATION_DIR/Results/3D
	cp -r $TMPDIR/$DNS_NAME/.recovery $SIMULATION_DIR/
	cp -r $TMPDIR/$DNS_NAME/Log $SIMULATION_DIR/
	
}

TEXT_LEVEL1="---------"
TEXT_LEVEL2="-------------"


#DNS_NAME=TEST2
#DNS_NAME=TEST_Laminar_OPEN
#DNS_NAME=TEST_Laminar_Olivier
#DNS_NAME=TEST_Laminar
DNS_NAME=TEST_Laminar_PostTreatment
DNS_OUTPUT=OUT_DNS
DNS_TIME=345000
#JOB_NAME=job
#DNS_PROCS=4
#PROC_BY_NODE=16
PROW=1
PCOL=1
NPROCS=1


DNS_CODE=/home/benj/WORKSPACE/Codes/DNS/MULTIFAST_MHD_v2
TMPDIR=/home/benj/WORKSPACE/Codes/DNS/MULTIFAST_MHD_v2/TMP/$DNS_NAME
SIMULATION_DIR=/home/benj/WORKSPACE/Codes/DNS/MULTIFAST_MHD_v2/Simulations/$DNS_NAME
IT1=0
EVERY=100
IT2=5000


dir_to_create_list=$TMPDIR/mkdir_dir_list

build_running_env
#list_dir_to_create
create_results_dir
launch_simulation
#save_simulation
