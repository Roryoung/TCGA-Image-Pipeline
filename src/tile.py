import os
from optparse import OptionParser

import numpy as np
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator
import staintools
from PIL import Image
import cv2


def get_background_percentage(tile):
    """
        This function calculates the amount of white background in a given tile.

        Args:
            - tile: A PIL Image object of the slide tile

        Returns:
            - Float: The precentage of background in the tile
    """

    #convert tile to strict black and white (not greyscale)
    grey = tile.convert("L")
    bw = grey.point(lambda x: 0 if x>= 225 else 255)
    
    #floodfill black and white image
    des = np.array(bw)
    contour, hier = cv2.findContours(des,cv2.RETR_CCOMP,cv2.CHAIN_APPROX_SIMPLE)

    for cnt in contour:
        cv2.drawContours(des,[cnt],0,255,-1)

    #calculate background %
    bw_arr = np.array(des)
    white_cells = np.count_nonzero(bw_arr == 0)
    black_cells = np.count_nonzero(bw_arr == 255)
    return white_cells/(white_cells+black_cells)


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


def tile(slide_loc, output_dir, background=0.2, size=255):
    """
        This function will save tiles of the given H&E stained slide at different zoom levels.

        Args:
            - slide_loc: A .svs file of the H&E stained slides
            - output_dir: The location where the tiles will be saved
            - background: The maximum precentage of background allowed for a saved tile 
            - size: The width and hight of the tiles at each zoom level

        Returns:
            - None
    """
    slide = open_slide(slide_loc)
    dz = DeepZoomGenerator(slide, size, 0)


    file_name = os.path.basename(slide_loc).split(".")[0]
    slide_output_dir = os.path.join(output_dir, file_name)
    if not os.path.exists(slide_output_dir):
        os.makedirs(slide_output_dir)

    maxZoom = float(slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER]) / slide.level_downsamples[0]

    for level in range(1, dz.level_count):
        thisMag = maxZoom/pow(2,dz.level_count-(level+1))
        tile_dir = os.path.join(slide_output_dir, str(thisMag))
        if not os.path.exists(tile_dir):
            os.makedirs(tile_dir)

        cols, rows = dz.level_tiles[level]
        for row in range(rows):
            for col in range(cols):
                tile = dz.get_tile(level, (col, row))
                tile_name = os.path.join(tile_dir, f"{col}_{row}")

                avg_brightness = get_background_percentage(tile)
                #tile = normalize_tile(tile)

                if avg_brightness < background:
                    tile.save(f"{tile_name}.jpeg")

                else:
                    if np.random.uniform() < 0.1:
                        reject_dir = os.path.join(slide_output_dir, str(thisMag), "rejected")
                        if not os.path.exists(reject_dir):
                            os.makedirs(reject_dir)

                        reject_tile_name = os.path.join(reject_dir, f"{col}_{row}_{avg_brightness}")
                        tile.save(f"{reject_tile_name}.jpeg")
                        


if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide> <output_folder> [options]')
    parser.add_option('-o', '--output', metavar='NAME', dest='output_dir', help='base name of output file')
    parser.add_option('-b', '--background', metavar='PIXELS', dest='background', type='float', default=0.2, help='Max background threshold [0.2]; percentage of background allowed')
    parser.add_option('-s', '--size', metavar='PIXELS', dest='tile_size', type='int', default=255, help='tile size [255 ]')
    
    (opts, args) = parser.parse_args()

    try:
        slide_path = args[0]
    except IndexError:
        parser.error('Missing slide argument')

    if opts.output_dir is None:
        parser.error("Missing output directory argument")

    tile(slide_path, output_dir=opts.output_dir, background=opts.background, size=opts.tile_size)
