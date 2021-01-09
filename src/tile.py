import os
from optparse import OptionParser

import h5py
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

    def __init__(self, slide_loc, set_hdf5_file, normalizer=None, background=0.2,
                 size=255, reject_rate=0.1, ignore_repeat=False):
        """
            Args:
                - slide_loc: A .svs file of the H&E stained slides
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
        self.tiles = {}
        self.reject_tiles = {}

        proceed = "y"

        if self.file_name in set_hdf5_file:
            if not ignore_repeat:
                print(f"{self.file_name} is already in the dataset. Do you wish to overwrite these tiles? [y/n]")
                proceed = input()
            if proceed == "y":
                del set_hdf5_file[self.file_name]

        if proceed == "y":
            self.h5_group = set_hdf5_file.create_group(self.file_name)
            self._save_tiles()
            print()


    def _create_name_dataset(self, hdf5_file, name):
        names_db_shape = (0, 1)
        max_names_db_shape = (None, 1)
        dt = h5py.special_dtype(vlen=str)

        return hdf5_file.create_dataset(name=name, shape=names_db_shape, maxshape=max_names_db_shape, dtype=dt)


    def _create_image_dataset(self, hdf5_file, name, size, n_ch=3):
        img_db_shape = (0, size, size, n_ch)
        max_img_db_shape = (None, size, size, n_ch)
        dt = np.float32

        return hdf5_file.create_dataset(name=name, shape=img_db_shape, maxshape=max_img_db_shape, dtype=np.float32)


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

            zoom_hdf5 = self.h5_group.create_group(str(this_mag))
        
            img_storage = self._create_image_dataset(zoom_hdf5, 'images', self.size)
            name_storage = self._create_name_dataset(zoom_hdf5, 'file_name')

            reject_img_storage = self._create_image_dataset(zoom_hdf5, "reject_images", self.size)
            reject_name_storage = self._create_name_dataset(zoom_hdf5, "reject_file_name")

            print(f"\rCreating {self.file_name} | zoom: x{this_mag:.2f}", end="")
            for row in range(rows):
                for col in range(cols):
                    tile = np.array(self.dz.get_tile(level, (col, row)))
                    tile_name = f"{col}_{row}"

                    if self._keep_tile(tile, self.size, 1 - self.background):
                        if self.normalizer is not None:
                            self.normalizer.fit(tile)
                        
                        max_index = img_storage.shape[0] + 1
                        img_storage.resize(max_index, axis=0)
                        name_storage.resize(max_index, axis=0)

                        img_storage[max_index-1, :, :, :] = tile
                        name_storage[max_index-1, :] = tile_name

                    else:
                        if np.random.uniform() < self.reject_rate and tile.shape == (self.size, self.size, 3):
                            max_index = reject_img_storage.shape[0] + 1
                            reject_img_storage.resize(max_index, axis=0)
                            reject_name_storage.resize(max_index, axis=0)

                            reject_img_storage[max_index-1, :, :, :] = tile
                            reject_name_storage[max_index-1, :] = tile_name


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

# if __name__ == "__main__":
#     parser = OptionParser(usage='Usage: %prog <slide> <output_folder> [options]')
#     parser.add_option('-b', '--background', dest='background', type='float', default=0.2, help='Percentage of background allowed, default=0.2')
#     parser.add_option('-s', '--size', dest='tile_size', type='int', default=255, help='Size of the output tiles, default=255')
#     parser.add_option('-r', '--reject', dest='reject', type='float', default=0.1, help='Precentage of rejected background tiles to save, default=0.1')
#     parser.add_option('-i', '--ignore_repeat', dest='ignore_repeat', action="store_true", help='Automatically overwrte repeated files in the dataset, default=False')

#     (opts, args) = parser.parse_args()

#     try:
#         slide_path = args[0]
#     except IndexError:
#         parser.error('Missing slide argument')


#     tile = Tile(
#         slide_loc=slide_path,
#         background=opts.background,
#         size=opts.tile_size,
#         reject_rate=opts.reject,
#         ignore_repeat=opts.ignore_repeat
#     )
