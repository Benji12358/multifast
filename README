# MULTIFAST++

MULTIFAST++ is a parallel finite difference code that can perform a large number of flow simulations. The code has been parallelized thanks to openMP (Open Multi-Processing) and an elegant 2D domain decomposition strategy. Different temporal and numerical schemes have been implemented, covering a range from 2nd order accuracy with an explicit numerical scheme to 6th order accuracy with a compact scheme and a explicit optimized scheme (replicating quasi spectral accuracy). MULTIFAST++ contains a lot of modules allowing the use of passive scalar transport, the use of IBM (Immersed Boundary Method), the use of MHD (MagnetoHydrodynamic), the use of bubbles in the flow field or the use of the phase average. Examples of publications using this code can be found [here](https://doi.org/10.1016/j.compfluid.2014.10.009).

## Installation

To use MULTIFAST++, here are the steps to follow to setup the local environment.

#### STEP 1

First download all libraries listed in a folder libraries:

- gfortran/gcc/g++
- openmpi
- blas
- lapack
- fftw
- 2decomp_fft
- zlib
- hdf5

#### STEP 2 : Install openmpi

```bash
tar -xvf openmpi.tar.gz
cd openmpi/
mkdir -p /path/to/libraries/openmpi/openmpi
./configure --prefix=/path/to/libraries/openmpi/openmpi
make
make install
```

Then, openmpi needs to be added to the path. Proceed as follows:

```bash
vi ~/.bashrc
```

Add those lines at the end of the file:

```bash
export PATH="/path/to/libraries/openmpi/openmpi/bin:$PATH"
export LD_LIBRARY_PATH="/path/to/libraries/openmpi/openmpi/lib:$LD_LIBRARY_PATH"
```

Finally, update the path by sourcing the .bashrc file.

```bash
source ~/.bashrc
```

#### STEP 3 : Install blas

```bash
tar -xvf blas.tar.gz
cd BLAS
make
cp blas_LINUX.a libblas.a
mkdir -p /path/to/libraries/blas/blas
cp libblas.a /path/to/libraries/blas/blas/
```

#### STEP 4 : Install lapack

```bash
tar -xvf lapack.tgz
cd lapack
cp make.inc.example make.inc
vi make.inc
```

Edit the following line

```bash
BLASLIB=../../librefblas.a
```
to

```bash
BLASLIB=/path/to/libraries/blas/blas/libblas.a
```

Then,

```bash
make
mkdir -p /path/to/libraries/lapack/lapack
cp liblapack.a /path/to/libraries/lapack/lapack/
```

#### STEP 5 : Install fftw

```bash
tar -xvf fftw.tar.gz
cd fftw
mkdir -p /path/to/libraries/fftw/fftw
./configure --prefix=/path/to/libraries/fftw/fftw
make
make install
```

#### STEP 6 : Install 2decomp_fft

```bash
tar -xvf 2decomp_fft.tar.gz
cd 2decomp_fft
cd src
cp Makefile.inc.x86 Makefile.inc
vi Makefile.inc
```

Edit the following line

```bash
FFTW_PATH=/opt/fftw-3.3
```
to

```bash
FFTW_PATH=/path/to/libraries/fftw/fftw
```

Then,

```bash
cd ../
make
cd ../
cp -r 2decomp_fft /path/to/libraries/
```

#### STEP 7 : Install zlib

```bash
tar -xvf zlib.tar.gz
cd zlib
mkdir -p /path/to/libraries/zlib/zlib
./configure --prefix=/path/to/libraries/zlib/zlib
make
make check
make install
```

Then, zlib needs to be added to the path. Proceed as follows:

```bash
vi ~/.bashrc
```

Add those lines at the end of the file:

```bash
export PATH="/path/to/libraries/zlib/zlib/bin:$PATH"
export LD_LIBRARY_PATH="/path/to/libraries/zlib/zlib/lib:$LD_LIBRARY_PATH"
```

Finally, update the path by sourcing the .bashrc file.

```bash
source ~/.bashrc
```

#### STEP 8 : Install hdf5

```bash
tar -xvf hdf5.tar.gz
cd hdf5
mkdir -p /path/to/libraries/hdf5/hdf5
export CC=mpicc
export CCP="mpicc -E"
export CFLAGS="-O3"
export FC="mpif90"
export FCFLAGS="-O3"
export LIBS="-lz"
export LDFLAGS="-L/path/to/libraries/zlib/zlib/lib"
export CPPFLAGS="-I/path/to/libraries/zlib/zlib/include"
./configure --prefix=/path/to/libraries/hdf5/hdf5 --enable-fortran --enable-parallel --with-zlib=/path/to/libraries/zlib/zlib
make
make install
```

Then, hdf5 needs to be added to the path. Proceed as follows:

```bash
vi ~/.bashrc
```

Add those lines at the end of the file:

```bash
export PATH="/path/to/libraries/hdf5/hdf5/bin:$PATH"
export LD_LIBRARY_PATH="/path/to/libraries/hdf5/hdf5/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib64:/lib:/usr/lib:/usr/lib64
```

Finally, update the path by sourcing the .bashrc file.

```bash
source ~/.bashrc
```

#### STEP 9 : Compiling the code

Create a folder named WORKSPACE/Codes. Here, create a file libs.gfortran. This file will contain all the path to the local libraries which are necessary for the code to work successfully. Modify this file as follow:


```bash
LAPACK_LIB=/path/to/libraries/lapack/lapack-3.4.2/liblapack.a /home/benj/libraries/blas/blas-3.8.0/libblas.a

FLAG_FFT = -I /path/to/libraries/fftw/fftw/include/
LIBS_FFT = -L /path/to/libraries/fftw/fftw/lib/ -lfftw3 -lm

DECOMP2D_INC=-I /path/to/libraries/2decomp_fft/include
DECOMP2D_LIB=-L /path/to/libraries/2decomp_fft/lib -l2decomp_fft

HDF5_INC=-I /path/to/libraries/hdf5/hdf5/include
HDF5_LIB=/path/to/libraries/hdf5/hdf5/lib/libhdf5_fortran.a /path/to/libraries/hdf5/hdf5/lib/libhdf5.a -lz -ldl
```

Then, from here, create a folder DNS/MULTIFAST++. This folder will contain the code, so put the code here. Edit the following line of Makefile.inc as follows:

```bash
ROOT =/bettik/umairm/WORKSPACE
```
to

```bash
ROOT =/path/to/WORKSPACE
```

Finally, make the code:

```bash
make clean_all
make
```

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
