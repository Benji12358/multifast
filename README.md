# MULTIFAST++

MULTIFAST++ is a parallel finite difference code that can perform a large number of flow simulations. The code has been parallelized thanks to openMP (Open Multi-Processing) and an elegant 2D domain decomposition strategy. Different temporal and numerical schemes have been implemented, covering a range from 2nd order accuracy with an explicit numerical scheme to 6th order accuracy with a compact scheme and a explicit optimized scheme (replicating quasi spectral accuracy). MULTIFAST++ contains a lot of modules allowing the use of passive scalar transport, the use of IBM (Immersed Boundary Method), the use of MHD (MagnetoHydrodynamic), the use of bubbles in the flow field or the use of the phase average. Examples of publications using this code can be found [here](https://doi.org/10.1016/j.compfluid.2014.10.009).

## Use Apptainer

To use MULTIFAST++, here are the steps to follow to setup the local environment with the help of Apptainer (v1.3.2 recommended). The image you will built reproduce all the steps written in the README of the `multifast++` branch.

#### STEP 1 : Build the image

After cloning the repo and switching to the `apptainer` branch, run the following command:

```bash
apptainer build multifast.sif multifast.def
```

This should create the file `multifast.sif` that contains the software environment to run MULTIFAST.

#### STEP 2 : Compile the code

To compile the code, you should run the Apptainer image. To do so, run the following command:

```bash
apptainer shell --containall --bind $PWD:/workspace multifast.sif
```

> The `bind` option is binding the `$PWD` folder from the host machine (where the source code is), to the `/workspace` folder of the container. 

Once this is done, you should see that the bash prompt changed from `username@host` to `Apptainer` meaning you have entered the container. Now, you can run the following command:

```bash
cd /workspace # where the source code of MULTIFAST is
make clean_all
make
```

> To exit the container, simply run `exit`

#### STEP 3 : Use the code

To use the code, you should run the Apptainer image and bind a folder to store the data:

```bash
apptainer shell --containall --bind $PWD:/workspace --bind /data/sim_data:/sim_data multifast.sif
```

> A second `bind` option is binding the `/data/sim_data` folder from the host machine (where the simulation data will be stored), to the `/sim_data` folder of the container. 

Once you have entered the container, you can run the minimal example:

```bash
cd /workspace
./job_laptop_singleproc.sh test 0 100 100
```

> For more complex example, you can look at the files `job_multiprocs_apptainer.sh` or `dahu_example.oar`.

## Usage

```bash
to be continued ...
```

## Contributing

```bash
to be continued ...
```

## License

```bash
to be continued ...
```
