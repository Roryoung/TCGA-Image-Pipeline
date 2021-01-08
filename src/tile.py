import os
from optparse import OptionParser

import numpy as np
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator
from PIL import Image
import shutil


from scipy.ndimage.morphology import binary_fill_holes
from skimage.color import rgb2gray
from skimage.feature import canny
from skimage.morphology import binary_closing, binary_dilation, disk

class Tile:
    """
        This class will save tiles of the given H&E stained slide at different zoom levels.
    """

    def __init__(self, slide_loc, output_dir, normalizer=None, background=0.2,
                 size=255, reject_rate=0.1, ignore_repeat=False):
        """
            Args:
                - slide_loc: A .svs file of the H&E stained slides
                - output_dir: The location where the tiles will be saved
                - normalizer: A tile normalizer object
                - background: The maximum precentage of background allowed for a saved tile 
                - size: The width and hight of the tiles at each zoom level
                - reject_rate: The precentage of rejected tiles to save
                - ignore_repeat: Automatically overwrte repeated files in the dataset
        """
        self.normalizer = normalizer
        self.background = background
        self.size = size
        self.reject_rate = reject_rate

        self.slide = open_slide(slide_loc)
        self.dz = DeepZoomGenerator(self.slide, size, 0)

        self.file_name = ".".join(os.path.basename(slide_loc).split(".")[:-1])
        self.slide_output_dir = os.path.join(output_dir, self.file_name)
        self.tiles = {}
        self.reject_tiles = {}

        proceed = "y"

        if os.path.exists(self.slide_output_dir):
            if not ignore_repeat:
                print(f"{self.file_name} is already in the dataset. Do you wish to overwrite these tiles? [y/n]")
                proceed = input()
            if proceed == "y":
                shutil.rmtree(self.slide_output_dir)

        if proceed == "y":
            os.makedirs(self.slide_output_dir)
            self._save_tiles()
            print()


    def _save_tiles(self):
        """
            This function will save all the relevant tiles for a given zoom level.

            Args:
                - level: The level of zoom to save
                - tile_dir: The output directory for this zoom level

            Returns:
                - None
        """
        max_zoom = float(self.slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER]) / self.slide.level_downsamples[0]

        for level in range(1, self.dz.level_count):
            this_mag = max_zoom/pow(2,self.dz.level_count-(level+1))
            cols, rows = self.dz.level_tiles[level]

            tile_dir = os.path.join(self.slide_output_dir, str(this_mag))
            if not os.path.exists(tile_dir):
                os.makedirs(tile_dir)

            print(f"\rCreating {self.file_name} | zoom: x{this_mag:.2f}", end="")
            for row in range(rows):
                for col in range(cols):
                    tile = self.dz.get_tile(level, (col, row))

                    if self._keep_tile(tile, self.size, 1 - self.background):
                        if self.normalizer is not None:
                            self.normalizer.fit(tile)
                        
                        tile_name = os.path.join(tile_dir, f"{col}_{row}")
                        tile.save(f"{tile_name}.jpeg")

                    else:
                        if np.random.uniform() < self.reject_rate and tile.size == (self.size, self.size):
                            reject_dir = os.path.join(tile_dir, "rejected")
                            if not os.path.exists(reject_dir):
                                os.makedirs(reject_dir)

                            reject_tile_name = os.path.join(reject_dir, f"{col}_{row}")
                            tile.save(f"{reject_tile_name}.jpeg")


    def _keep_tile(self, tile, tile_size, tissue_threshold):
        """
        Determine if a tile should be kept.
        
        Args:
            - tile: A PIL Image object of the slide tile
            - tile_size: The width and height of a square tile to be generated.
            - tissue_threshold: Tissue percentage threshold.
        Returns:
            A Boolean indicating whether or not a tile should be kept.

        Check 0:
            The tile must be the specified height and width

        Check 1:
            - Convert image to greyscale with 0 = background, 1 = tissue
            - Canny edge detect
            - Binary dilation followed by erosion
            - Binary dilation to fill gaps in tissue
            - Calcualte tissue precentage and test against given minumum tissue value

        Check 2:
            - Convert tile to optical density space
            - Threshold values
            - Binary dilation followed by erosion
            - Binary dilation to fill gaps in tissue
            - Calcualte tissue precentage and test against given minumum tissue value
        """

        tile = np.array(tile)

        if tile.shape[0:2] == (tile_size, tile_size):
            tile_orig = tile
            tile = rgb2gray(tile)
            tile = 1 - tile
            tile = canny(tile)
            tile = binary_closing(tile, disk(10))
            tile = binary_dilation(tile, disk(10))
            tile = binary_fill_holes(tile)

            percentage_1 = tile.mean()
            check_1 = percentage_1 >= tissue_threshold

            
            tile = tile_orig.astype(np.float64)
            tile = -np.log((tile+1)/240) #convert to optical density
            beta = 0.15
            tile = np.min(tile, axis=2) >= beta
            tile = binary_closing(tile, disk(2))
            tile = binary_dilation(tile, disk(2))
            tile = binary_fill_holes(tile)
            percentage_2 = tile.mean()
            check_2 = percentage_2 >= tissue_threshold

            return check_1 and check_2
        else:
            return False

if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide> <output_folder> [options]')
    parser.add_option('-b', '--background', dest='background', type='float', default=0.2, help='Percentage of background allowed, default=0.2')
    parser.add_option('-s', '--size', dest='tile_size', type='int', default=255, help='Size of the output tiles, default=255')
    parser.add_option('-r', '--reject', dest='reject', type='float', default=0.1, help='Precentage of rejected background tiles to save, default=0.1')
    parser.add_option('-i', '--ignore_repeat', dest='ignore_repeat', action="store_true", help='Automatically overwrte repeated files in the dataset, default=False')

    (opts, args) = parser.parse_args()

    try:
        slide_path = args[0]
    except IndexError:
        parser.error('Missing slide argument')

    try:
        output_dir = args[1]
    except IndexError:
        parser.error('Missing output directory argument')

    tile = Tile(
        slide_loc=slide_path,
        output_dir=output_dir,
        background=opts.background,
        size=opts.tile_size,
        reject_rate=opts.reject,
        ignore_repeat=opts.ignore_repeat
    )
