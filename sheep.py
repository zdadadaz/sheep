import numpy as np
import h5py
import matplotlib.pyplot as plt


FILE = "sheep.h5"
DATASET = "DS1"


hf = h5py.File(FILE, 'r')
sheep = hf[DATASET]
sheep = np.array(sheep)

N, T = sheep.shape

grid_size = int(np.sqrt(N)*2)
# random non-repetitive coordinates
coords = np.array([(i, j) for i in range(grid_size) for j in range(grid_size)])
idx = np.random.choice(grid_size*grid_size, size=N, replace=False)
coords = coords[idx,:]
print(coords)

# plot the first frame and save the markers
marker = [0]*N
for i, coord in enumerate(coords):
    marker[i] = plt.scatter(coord[0], coord[1], color='blue', alpha=1, s=10**2)

# for t in range(1, T):
#     for i in range(N):
#         if time_series[i][t] > 0.5:
#             alpha = 1
#         else: # past half-life
#             alpha = 0.5
#             marker[i].set_alpha(alpha)
#     plt.title('t={}'format(str(t)))
#     plt.pause(0.01)
    
plt.show()
