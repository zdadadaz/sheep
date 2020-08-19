import numpy as np
import h5py
import matplotlib.pyplot as plt
import math

FILE = "wolve.h5"
DATASET = "DS1"


hf = h5py.File(FILE, 'r')
sheep = hf[DATASET]
sheep = np.array(sheep)

N, T = sheep.shape

grid_size = int(np.sqrt(N))
marker = [0]*N

t=1
half = math.sqrt(N)
for t in range(0,100):
    qq = {}
    count = 0
    for i in range(N):
        if sheep[i][t] >0:
            count += 1
            qq[i] = sheep[i][t]
            y = i//half
            x = i - y*half
            marker[i] = plt.scatter(x, y, color='blue', alpha=1, s=10**2)
            plt.xlim(-1, grid_size)
            plt.ylim(-1, grid_size)
            plt.annotate(round(sheep[i][t]), (x, y))
    plt.title('t: {}, sheep num: {}'.format(str(t),str(count)))
    plt.pause(1)
    plt.clf()
    
plt.show()