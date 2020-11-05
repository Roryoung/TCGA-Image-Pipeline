import os
from optparse import OptionParser

import numpy as np
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator


def get_avg_brightness(tile):
    grey = tile.convert("L")
    bw = grey.point(lambda x: 0 if x< 200 else 1, "F")
    return np.average(bw)


def tile(slide_loc, output_dir, background=0.1, size=299):
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

                avg_brightness = get_avg_brightness(tile)

                if avg_brightness < background:
                    tile.save(tile_name + ".jpeg")
                


if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide> <output_folder> [options]')
    parser.add_option('-o', '--output', metavar='NAME', dest='output_dir', help='base name of output file')
    parser.add_option('-b', '--background', metavar='PIXELS', dest='background', type='float', default=0.1, help='Max background threshold [0.1]; percentage of background allowed')
    parser.add_option('-s', '--size', metavar='PIXELS', dest='tile_size', type='int', default=299, help='tile size [299]')
    
    (opts, args) = parser.parse_args()

    try:
        slide_path = args[0]
    except IndexError:
        parser.error('Missing slide argument')

    if opts.output_dir is None:
        parser.error("Missing output directory argument")

    tile(slide_path, output_dir=opts.output_dir, background=opts.background, size=opts.tile_size)
