from util import *
import numpy
import math
import cv2
import matplotlib.pyplot as plt



outWidth = 1200
fps =10
plot_path = 'plot.mp4'
visual_path = 'visual'
path = "./animal.h5"
sheep_path = 'sheep.jpg'
wolf_path = 'wolf.jpg'

sheep, wolve, grass, count, setting = read_hdf(path)
N = setting[0]
T = setting[1]
# setting
realW = int(math.sqrt(N))
step = int((outWidth)//(realW))

# preprocssing for images
sheep_im = preprocessing_img(sheep_path, step, 10)
wolf_im = preprocessing_img(wolf_path, step, 2)


# save number of plot
# print(T)
# print("T, count ",T, len(count))
# subsample = 1
# cc = []
# for t in range(0,T,subsample):
#     cc.append(count[t*3])
#     cc.append(count[1+t*3])
#     cc.append(count[2+t*3])
# T = len(cc)//3
# count = cc
save_number_trend_plot(count,N,T,fps, plot_path)

# save visualization
cnt = 0
count = T
while count > 0:
    if count > 100:
        tt = 100
        count -= 100
    else:
        tt = count
        count = 0
    # print(count)
    canvas_v = np.zeros((tt,outWidth,outWidth,3))
    for t in range(tt):
        # pool.map(draw_one_canvas, zip(video_paths, means, stds, periods))
        canvas = np.zeros((outWidth,outWidth,3))
        draw_one_canvas(canvas, sheep[(t*N):((t+1)*N)], sheep_im,  wolve[(t*N):((t+1)*N)], wolf_im, grass[(t*N):((t+1)*N)], outWidth, step,realW)
        canvas = cv2.cvtColor(canvas.astype(np.uint8), cv2.COLOR_BGR2RGB)
        canvas_v[t,...] = canvas
    canvas_v = canvas_v.astype(np.uint8)
    write_video(canvas_v.transpose((3,0,1,2)),visual_path+"_"+str(cnt)+".avi",outWidth, fps)
    print("write out {}", visual_path+"_"+str(cnt)+".avi")
    cnt+= 1

