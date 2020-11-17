import os
from optparse import OptionParser

import numpy as np
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator
import staintools
from PIL import Image
import cv2
import shutil

class Tile:
    """
        This class will save tiles of the given H&E stained slide at different zoom levels.
    """

    def __init__(self, slide_loc, output_dir, background=0.2, threshold=225, size=255, reject_rate=0.1, ignore_repeat=False):
        """
            Args:
                - slide_loc: A .svs file of the H&E stained slides
                - output_dir: The location where the tiles will be saved
                - background: The maximum precentage of background allowed for a saved tile 
                - threshold: The maximum greyscale value for a pixel to be considered background
                - size: The width and hight of the tiles at each zoom level
                - reject_rate: The precentage of rejected tiles to save
                - ignore_repeat: Automatically overwrte repeated files in the dataset
            Returns:
                - None
        """
        self.background = background
        self.threshold = threshold
        self.size = size
        self.reject_rate = reject_rate

        self.slide = open_slide(slide_loc)
        self.dz = DeepZoomGenerator(self.slide, size, 0)

        file_name = os.path.basename(slide_loc).split(".")[0]
        self.slide_output_dir = os.path.join(output_dir, file_name)

        proceed = "y"

        if os.path.exists(self.slide_output_dir):
            if not ignore_repeat:
                print(f"{file_name} is already in the dataset. Do you wish to overwrite these tiles? [y/n]")
                proceed = input()
            if proceed == "y":
                shutil.rmtree(self.slide_output_dir)

        if proceed == "y":
            os.makedirs(self.slide_output_dir)

            maxZoom = float(self.slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER]) / self.slide.level_downsamples[0]

            for level in range(1, self.dz.level_count):
                thisMag = maxZoom/pow(2,self.dz.level_count-(level+1))
                tile_dir = os.path.join(self.slide_output_dir, str(thisMag))
                print("\rTiling slide {}, Zoom level {:.2f}".format(file_name, thisMag), end="")

                self.tile_level(level, tile_dir)

            print()


    def tile_level(self, level, tile_dir):
        """
            This function will save all the relevant tiles for a given zoom level.

            Args:
                - level: The level of zoom to save
                - tile_dir: The output directory for this zoom level

            Returns:
                - None
        """
        if not os.path.exists(tile_dir):
            os.makedirs(tile_dir)
        
        reject_dir = os.path.join(tile_dir, "rejected")
        if not os.path.exists(reject_dir):
            os.makedirs(reject_dir)

        cols, rows = self.dz.level_tiles[level]
        for row in range(rows):
            for col in range(cols):
                tile = self.dz.get_tile(level, (col, row))
                tile_name = os.path.join(tile_dir, f"{col}_{row}")

                tile_background = self.get_background_percentage(tile)
                #tile = normalize_tile(tile)

                if tile_background < self.background:
                    tile.save(f"{tile_name}.jpeg")

                else:
                    if np.random.uniform() < self.reject_rate:
                        reject_tile_name = os.path.join(reject_dir, f"{col}_{row}_{tile_background}")
                        tile.save(f"{reject_tile_name}.jpeg")

    def get_background_percentage(self, tile):
        """
            This function calculates the amount of white background in a given tile.

            Args:
                - tile: A PIL Image object of the slide tile

            Returns:
                - Float: The precentage of background in the tile
        """

        #convert tile to strict black and white (not greyscale)    
        grey = tile.convert("L")
        bw = grey.point(lambda x: 0 if x< self.threshold else 255)

        bw_arr = np.array(bw)
        foreground_cells = np.count_nonzero(bw_arr == 0)
        background_cells = np.count_nonzero(bw_arr == 255)
        return background_cells/(foreground_cells+background_cells)

    def normalize_tile(self, tile):
        """
            This function normalizes the staining of H&E histology slides using the macenko method.

            Args:
                - tile: A PIL Image object of the slide tile
            
            Returns:
                - A PIL Image object of the normalized slide tile
        """
        tile_array = np.array(tile)
        tile_array = staintools.LuminosityStandardizer.standardize(tile_array)

        normalizer = staintools.StainNormalizer(method='macenko')
        normalizer.fit(tile_array)
        normalized_tile = normalizer.transform(tile_array)

        return Image.fromarray(normalized_tile)



if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide> <output_folder> [options]')
    parser.add_option('-o', '--output', metavar='NAME', dest='output_dir', help='base name of output file')
    parser.add_option('-b', '--background', metavar='PIXELS', dest='background', type='float', default=0.2, help='Percentage of background allowed [0.2]')
    parser.add_option('-t', '--threshold', metavar='', dest='threshold', type='int', default=225, help='Backgorund threshold [225]')
    parser.add_option('-s', '--size', metavar='PIXELS', dest='tile_size', type='int', default=255, help='tile size [255]')
    parser.add_option('-r', '--reject', dest='reject', type='float', default=0.1, help='Precentage of rejected background tiles to save [0.1]')
    parser.add_option('-i', '--ignore_repeat', dest='ignore_repeat', action="store_true", help='Automatically overwrte repeated files in the dataset [False]')

    (opts, args) = parser.parse_args()

    try:
        slide_path = args[0]
    except IndexError:
        parser.error('Missing slide argument')

    if opts.output_dir is None:
        parser.error("Missing output directory argument")

    tile = Tile(
        slide_loc=slide_path,
        output_dir=opts.output_dir,
        background=opts.background,
        size=opts.tile_size,
        threshold=opts.threshold,
        reject_rate=opts.reject,
        ignore_repeat=opts.ignore_repeat
    )
