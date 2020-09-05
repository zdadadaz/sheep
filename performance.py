import subprocess
import numpy as np

def run_cmd(cmd):
    time = ''
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell = True)
    count = 0
    for line in iter(p.stdout.readline, b''):
        if line and count == 1: # Don't print blank lines
            time = line.decode("utf-8").split('\t')[1]
            break
        count+=1
    return float(time[2:-2])

Ns = [100,400, 900, 1600, 2500, 3600, 4900, 6400, 8100, 10000]
T = 1000
initSheeps = [4,16,36,64,100,144,192,256,324,400]
initWolfs = [2,8,18,32,50,72,96,128,162,200]
initGrasses = [48,192,432,768,1200,1728,2352,3072,3888,4800]

mean_time = []
for idx, (N, initSheep, initWolf, initGrass) in enumerate(zip(Ns, initSheeps, initWolfs, initGrasses)):
    # N = 2500
    # initSheep= 100
    # initWolf = 50
    # initGrass =1200
    cmd = "time ./sheep " + str(N) +" " + str(T)+" " + str(initSheep)+" " + str(initWolf)+" " + str(initGrass)
    times= []
    for i in range(10):
        times.append(run_cmd(cmd))
    mean_time.append(np.array(times).mean()*1000)
    print("time :",idx, mean_time[-1])
print(mean_time)
