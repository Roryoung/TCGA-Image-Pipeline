import numpy as np
from skimage import color
from PIL import Image

class Normalizer:
    def __init__(self):
        self.means = np.empty((0, 3))        
        self.stds = np.empty((0,3))


    def fit(self, tile):
        tile_array = np.array(tile)
        lab = color.rgb2lab(tile_array)

        self.means = np.append(self.means, [[np.mean(lab[:,:,i]) for i in range(3)]], axis=0)
        self.stds = np.append(self.stds, [[np.std(lab[:,:,i]) for i in range(3)]], axis=0)


    def normalize(self, tile):
        avg_mean = np.mean(self.means, axis=0)
        avg_stds = np.mean(self.stds, axis=0)

        tile_array = np.array(tile)
        lab = color.rgb2lab(tile_array)

        t_mean = [0,0,0]
        t_std  = [1,1,1]

        # Each channel 
        for i in range(3):
            t_mean[i] = np.mean(lab[:,:,i])
            t_std[i]  = np.std(lab[:,:,i])
            tmp = ( (lab[:,:,i] - t_mean[i]) * (avg_stds[i] / t_std[i]) ) + avg_mean[i]
            if i == 0:
                tmp[tmp<0] = 0
                tmp[tmp>100] = 100
                lab[:,:,i] = tmp
            else:
                tmp[tmp<-127] = -127
                tmp[tmp>127] = 127
                lab[:,:,i] = tmp

        return Image.fromarray((color.lab2rgb(lab)*255).astype(np.uint8))