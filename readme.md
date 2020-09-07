# File descriptions
## Main function: 
sheep.c
## Visualization: 
makevideo.py: generate videos for plot and Visualization
plotNum.py: plot the static population dynamic
# Required package
install python3 package in requirements.txt

# Run program
1. Set parameter in sheep.c if you want to try different parameters. 
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

2. make
3. ./sheep
4. python3 makevideo.py 
5. python3 plotNum.py

# Note
when T>100, plot.mp4 video will be subsample and visualization video will be divided into several videos with maximum 100 frames per videos. The initial setting for output number of visualization video is 1.



