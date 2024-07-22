#!/bin/bash

launch_simulation(){
	cd $TMPDIR_host
	
	echo $TEXT_LEVEL1 Running simulation...

	mpirun -np `cat $OAR_FILE_NODES|wc -l` \
                --machinefile $OAR_NODE_FILE \
                -mca plm_rsh_agent "oarsh" \
                --prefix $HOME/.nix-profile \
                apptainer exec instance://multifast \
                /bin/bash -c "cd $TMPDIR; $TMPDIR/DNS_EXEC $DNS_NAME $DNS_TIME $PROW $PCOL 0 > $DNS_OUTPUT 2>temp_arrays.log"

	cd -
}

create_results_dir()
{
	# Create subfolders for the main channel
	cd $TMPDIR_host/$DNS_NAME/Results/3D	

	for numero in `seq $IT1 $EVERY $IT2`
    		do mkdir ./field"$numero"
	done

	cd -

	# Create subfolders for the embedded channel
	cd $TMPDIR_host/$DNS_NAME/Results/Embedded

	for numero in `seq $IT1 $EVERY $IT2`
		do mkdir ./field"$numero"
	done

	cd -

	# Create subfolders for the following channel
	cd $TMPDIR_host/$DNS_NAME/Results/following

	for numero in `seq $IT1 $EVERY $IT2`
		do mkdir ./field"$numero"
	done

	cd -
}

build_running_env(){

	echo $TEXT_LEVEL1 Building running environnment...

	rm -rf $TMPDIR_host
	mkdir $TMPDIR_host
	cp -r $SIMULATION_DIR_host/.recovery/arborescence $TMPDIR_host/$DNS_NAME

	cp -r $SIMULATION_DIR_host/Input $TMPDIR_host/$DNS_NAME
	cp -r $SIMULATION_DIR_host/Log $TMPDIR_host/$DNS_NAME
	cp -r $SIMULATION_DIR_host/.recovery $TMPDIR_host/$DNS_NAME

	cp $DNS_CODE_host/DNS_EXEC $TMPDIR_host

}

TEXT_LEVEL1="---------"
TEXT_LEVEL2="-------------"

DNS_OUTPUT=OUT_DNS
DNS_TIME=700000
PROW=0
PCOL=0

DNS_NAME=$1

DNS_CODE_host=/home/arrondeb/WORKSPACE/Codes/DNS/multifast
TMPDIR_host=/bettik/arrondeb/sim_data/$DNS_NAME
TMPDIR=/sim_data/$DNS_NAME
SIMULATION_DIR_host=/home/arrondeb/WORKSPACE/Codes/DNS/multifast/Simulations/$DNS_NAME
IT1=$2
EVERY=$4
IT2=$3


build_running_env
create_results_dir
launch_simulation
