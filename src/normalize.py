import os
from optparse import OptionParser

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


    def fit_dir(self, current_path):
        for filename in os.listdir(current_path):
            if filename != "rejected":
                new_path = os.path.join(current_path, filename)
                if os.path.isfile(new_path):
                    tile = Image.open(new_path)
                    tile = self.fit(tile)
                else:
                    self.fit_dir(new_path)



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


    def normalize_dir(self, current_path):
        for filename in os.listdir(current_path):
            if filename != "rejected":
                new_path = os.path.join(current_path, filename)
                if os.path.isfile(new_path):
                    tile = Image.open(new_path)
                    tile = self.normalize(tile)
                    tile.save(new_path)
                else:
                    self.normalize_dir(new_path)

if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <tile_dir>')
    (opts, args) = parser.parse_args()

    try:
        tile_dir = args[0]
    except IndexError:
        parser.error('Missing tile directory argument')

    normalizer = Normalizer()
    normalizer.fit_dir(tile_dir)
    normalizer.normalize_dir(tile_dir)