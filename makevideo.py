from util import *
import numpy
import math
import cv2
import matplotlib.pyplot as plt


N = 2500
T = 100
outWidth = 1200
fps =10
path = "./animal.h5"
sheep_path = 'sheep.jpg'
wolf_path = 'wolf.jpg'

sheep, wolve, grass, count = read_hdf(path)

# setting
realW = int(math.sqrt(N))
step = int((outWidth)//(realW))

# preprocssing for images
sheep_im = preprocessing_img(sheep_path, step, 10)
wolf_im = preprocessing_img(wolf_path, step, 2)

# save number of plot
save_number_trend_plot(count,N,T,fps)

# save visualization
canvas_v = np.zeros((T,outWidth,outWidth,3))
for t in range(T):
    canvas = np.zeros((outWidth,outWidth,3))
    draw_one_canvas(canvas, sheep[(t*N):((t+1)*N)], sheep_im,  wolve[(t*N):((t+1)*N)], wolf_im, grass[(t*N):((t+1)*N)], outWidth, step,realW)
    canvas = cv2.cvtColor(canvas.astype(np.uint8), cv2.COLOR_BGR2RGB)
    canvas_v[t,...] = canvas

canvas_v = canvas_v.astype(np.uint8)
write_video(canvas_v.transpose((3,0,1,2)),"tmp.avi",outWidth, fps)


