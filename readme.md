# File descriptions
## Main function: 
There are 4 four files for each porpuse, but only hybrid version has complete comments.
* Serial: sheep_v7_serial.cpp
* OpenMP: sheep_v7_omp.cpp
* MPI: sheep_v7_mpi.cpp
* Hybrid: sheep_v7_hyb.cpp
## Visualization: 
* makevideo.py: generate videos for plot and Visualization( only for serial)
* plotNum.py: plot the static population dynamic( only for serial)
* plotNum_txt.py: plot the static population dynamic ( can be used for all file, needed turn on **outputPopulation** in definition)
## Testing: 
* run_test.py: generate scripts and sbatch to slurm
* parse_slurm.py: collecting output information and print on the screen

# Required package
* install hdf5 library for c as the steps in http://micro.stanford.edu/wiki/Install_HDF5#Build_and_Installation_from_Sources
* install python3 package in requirements.txt

# Run program
This program is tested and run on Moss and Getfix. 

1. Set parameter in each sheep_*.cpp file if you want to try different parameters. 
N: size of landscape, should be square
T: time
initSheepNum: init sheep number
sheepGainFromFood: sheep gain energe from consumption of grass
sheepReproduce: probability of reproduction(%)
initWolveNum: init wolves number
wolveGainFromFood:Wolves gain energe from consumption of sheep
wolveReproduce: probability of reproduction(%)
Grass: (1/0)grass flag for version 1 or version 2 model.
initGrass: init grass number
grassRegrowth: time period for gass regrowth

    ex: <code>time ./executefile N T initSheepNum initWolveNum initGrass </code>

2. make
3. serial: ./sheep, OpenMp: ./sheep_omp, MPI: ./sheep_mpi, Hybrid: ./sheep_hyb
4. python3 makevideo.py (turn on **visualization** definition)
5. python3 plotNum.py (turn on **visualization** definition)
6. python3 plotNum_txt.py (turn on **outputPopulation** definition)

# Note
1. when T>100, plot.mp4 video will be subsample and visualization video will be divided into several videos with maximum 100 frames per videos. The initial setting for output number of visualization video is 1.
2. Only serial version can do visualization.
3. outputPopulation can be used in all files to check the correctness



