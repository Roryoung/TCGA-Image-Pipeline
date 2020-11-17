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


def get_background_percentage(tile, threshold=225):
    """
        This function calculates the amount of white background in a given tile.

        Args:
            - tile: A PIL Image object of the slide tile
            - threshold: The maximum greyscale value for a pixel to be considered background

        Returns:
            - Float: The precentage of background in the tile
    """

    #convert tile to strict black and white (not greyscale)    
    grey = tile.convert("L")
    bw = grey.point(lambda x: 0 if x< 225 else 255)

    bw_arr = np.array(bw)
    foreground_cells = np.count_nonzero(bw_arr == 0)
    background_cells = np.count_nonzero(bw_arr == 255)
    return background_cells/(foreground_cells+background_cells)


def normalize_tile(tile):
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


def tile_level(dz, level, tile_dir, background, threshold):
    """
        This function will save all the relevant tiles for a given zoom level.

        Args:
            - dz: A DeepZoomGenerated object for the slide
            - level: The level of zoom to save
            - tile_dir: The output directory for this zoom level
            - background: The maximum precentage of background allowed for a saved tile 
            - threshold: The maximum greyscale value for a pixel to be considered background

        Returns:
            - None
    """
    if not os.path.exists(tile_dir):
        os.makedirs(tile_dir)
    
    reject_dir = os.path.join(tile_dir, "rejected")
    if not os.path.exists(reject_dir):
        os.makedirs(reject_dir)

    cols, rows = dz.level_tiles[level]
    for row in range(rows):
        for col in range(cols):
            tile = dz.get_tile(level, (col, row))
            tile_name = os.path.join(tile_dir, f"{col}_{row}")

            tile_background = get_background_percentage(tile, threshold)
            #tile = normalize_tile(tile)

            if tile_background < background:
                tile.save(f"{tile_name}.jpeg")

            else:
                if np.random.uniform() < 0.1:

                    reject_tile_name = os.path.join(reject_dir, f"{col}_{row}_{tile_background}")
                    tile.save(f"{reject_tile_name}.jpeg")


def tile(slide_loc, output_dir, background=0.2, threshold=225, size=255):
    """
        This function will save tiles of the given H&E stained slide at different zoom levels.

        Args:
            - slide_loc: A .svs file of the H&E stained slides
            - output_dir: The location where the tiles will be saved
            - background: The maximum precentage of background allowed for a saved tile 
            - threshold: The maximum greyscale value for a pixel to be considered background
            - size: The width and hight of the tiles at each zoom level

        Returns:
            - None
    """
    slide = open_slide(slide_loc)
    dz = DeepZoomGenerator(slide, size, 0)

    file_name = os.path.basename(slide_loc).split(".")[0]
    slide_output_dir = os.path.join(output_dir, file_name)

    if os.path.exists(slide_output_dir):
        shutil.rmtree(slide_output_dir)
        
    os.makedirs(slide_output_dir)

    maxZoom = float(slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER]) / slide.level_downsamples[0]

    for level in range(1, dz.level_count):
        thisMag = maxZoom/pow(2,dz.level_count-(level+1))
        tile_dir = os.path.join(slide_output_dir, str(thisMag))

        tile_level(dz, level, tile_dir, background, threshold)
                        


if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide> <output_folder> [options]')
    parser.add_option('-o', '--output', metavar='NAME', dest='output_dir', help='base name of output file')
    parser.add_option('-b', '--background', metavar='PIXELS', dest='background', type='float', default=0.2, help='Percentage of background allowed [0.2]')
    parser.add_option('-t', '--threshold', metavar='', dest='threshold', type='int', default=225, help='Backgorund threshold [225 ]')
    parser.add_option('-s', '--size', metavar='PIXELS', dest='tile_size', type='int', default=255, help='tile size [255 ]')

    (opts, args) = parser.parse_args()

    try:
        slide_path = args[0]
    except IndexError:
        parser.error('Missing slide argument')

    if opts.output_dir is None:
        parser.error("Missing output directory argument")

    tile(slide_path, output_dir=opts.output_dir, background=opts.background, size=opts.tile_size)
