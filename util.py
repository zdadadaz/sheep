import numpy as np
import h5py
import cv2
from matplotlib import animation
import matplotlib.pyplot as plt

def read_hdf(path):
    DATASETs = "Ds_sheep"
    DATASETw = "Ds_wolve"
    DATASETg = "Ds_grass"
    DATASETc = "Ds_count"
    DATASETset = "Ds_set"
    hf = h5py.File(path, 'r')
    sheep = hf[DATASETs]
    wolve = hf[DATASETw]
    grass = hf[DATASETg]
    count = hf[DATASETc]
    setting = hf[DATASETset]
    sheep = np.array(sheep)
    wolve = np.array(wolve)
    grass = np.array(grass)
    count = np.array(count)
    setting = np.array(setting)
    return (sheep, wolve, grass, count, setting)

def readimg(path):
    sheep_im = cv2.imread(path)
    sheep_im = cv2.cvtColor(sheep_im, cv2.COLOR_BGR2RGB).astype(int)
    return sheep_im

def saveimg(path,sheep_im):
    sheep_im[sheep_im<0] = 0
    sheep_im[sheep_im>255] = 255
    sheep_im = sheep_im.astype(np.uint8)
    cv2.imwrite(path, cv2.cvtColor(sheep_im, cv2.COLOR_RGB2BGR)) 

def remove_background(sheep_im, diff):
    # diff = 10
    for i in range(sheep_im.shape[0]):
        for j in range(sheep_im.shape[1]):
            if (abs(sheep_im[i,j,0]-sheep_im[i,j,1])<diff) and (abs(sheep_im[i,j,0]-sheep_im[i,j,2])<diff) and (abs(sheep_im[i,j,2]-sheep_im[i,j,1])<diff):
                sheep_im[i,j,0] = sheep_im[i,j,1] = sheep_im[i,j,2] = 0
    return sheep_im

def preprocessing_img(path, width, diff):
    sheep_im = readimg(path)
    sheep_im = remove_background(sheep_im, diff)
    sheep_im[sheep_im<0] = 0
    sheep_im[sheep_im>255] = 255
    sheep_im = sheep_im.astype(np.uint8)
    sheep_im = cv2.resize(sheep_im, (width,width),interpolation = cv2.INTER_AREA)
    return sheep_im

def draw_img(img, animal, x, y):
    for i in range(animal.shape[0]):
        for j in range(animal.shape[1]):
            if animal[i,j,0]!=0:
                img[x+i,y+j,:] = animal[i,j,:]
                
def draw_one_canvas(canvas, sheep, sheep_im,  wolf, wolf_im, grass, outWidth, step, half):
    green = np.zeros((step,step,3))
    green[:,:,1] = np.ones((step,step))*200
    brown = np.ones((step,step,3))
    brown[:,:,0] *= 150
    brown[:,:,1] *= 75
    brown[:,:,2] *= 0
    for idx, i in enumerate(range(0, outWidth-step//2, step)):
        for idj, j in enumerate(range(0, outWidth-step//2, step)):
            # print(idx,idj)
            if grass[idx*half + idj]==1:
                # print(i+step, j+step, canvas.shape, green.shape)
                canvas[i:(i+step),j:(j+step),:] = green
            else:
                canvas[i:(i+step),j:(j+step),:] = brown
            if sheep[idx*half + idj]>0:
                draw_img(canvas, sheep_im, i, j)
            if wolf[idx*half + idj]>0:
                draw_img(canvas, wolf_im, i, j)
def write_video(array,filename,width,fps):
    c, f, height, width = array.shape
    fourcc = cv2.VideoWriter_fourcc('M', 'J', 'P', 'G')
    out = cv2.VideoWriter(filename, fourcc, fps, (width, height))
    for i in range(f):
        out.write(array[:, i, :, :].transpose((1, 2, 0)))

def save_number_trend_plot(count,N, T,fps, filename):
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
    plt.rcParams['animation.html'] = 'jshtml'
    fig2 = plt.figure(2)
    ax=fig2.add_subplot(111)

    xx = [i for i in range(T)]
    yys = []
    yyw = []
    yyg = []
    for t in range(T):
        yys.append(count[0 + t*3])
        yyw.append(count[1 + t*3])
        yyg.append(count[2 + t*3])
    def animate(i):
        ax.clear()
        ax.set_title('sheep: {}, wolf: {}, grass/4: {}'.format(yys[i],yyw[i],yyg[i]))
        ax.set_xlim(-1, T)
        ax.set_ylim(-1, N//4)
        ax.plot(xx[:i], yys[:i],color='b')
        ax.plot(xx[:i], yyw[:i],color='r')
        ax.plot(xx[:i], np.array(yyg[:i])/4,color='g')
        ax.grid(b=True, which='major', color='k', linestyle='--',alpha=0.1)
        plt.grid()
    im_ani = animation.FuncAnimation(fig2, animate, interval=T)
    im_ani.save(filename, writer=writer)
