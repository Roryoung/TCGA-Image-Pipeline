import os
from optparse import OptionParser
import numpy as np

from tile import Tile
from normalize import Normalizer

def build_dataset(slide_dir, output_dir, background=0.2, threshold=225, size=255, reject_rate=0.1, ignore_repeat=False):
    normalizer = Normalizer()
    
    slide_tiles = []
    for filename in os.listdir(slide_dir):
        slide_path = os.path.join(slide_dir, filename)
        tile = Tile(
            slide_loc=slide_path,
            output_dir=output_dir,
            normalizer=normalizer,
            background=background,
            threshold=threshold,
            size=size,
            reject_rate=reject_rate,
            ignore_repeat=ignore_repeat
        )

    normalizer.normalize_dir(output_dir)

if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide_folder> <output_folder> [options]')
    parser.add_option('-b', '--background', dest='background', type='float', default=0.2, help='Percentage of background allowed, defualt=0.2')
    parser.add_option('-t', '--threshold', dest='threshold', type='int', default=225, help='Backgorund threshold, defualt=225')
    parser.add_option('-s', '--size', dest='tile_size', type='int', default=255, help='tile size, defualt=255')
    parser.add_option('-r', '--reject', dest='reject', type='float', default=0.1, help='Precentage of rejected background tiles to save, defualt=0.1')
    parser.add_option('-i', '--ignore_repeat', dest='ignore_repeat', action="store_true", help='Automatically overwrte repeated files in the dataset, defualt=False')

    (opts, args) = parser.parse_args()

    try:
        slide_dir = args[0]
    except IndexError:
        parser.error('Missing input directory argument')

    try:
        output_dir = args[1]
    except IndexError:
        parser.error('Missing input directory argument')

    build_dataset(
        slide_dir=slide_dir,
        output_dir=output_dir,
        background=opts.background,
        size=opts.tile_size,
        threshold=opts.threshold,
        reject_rate=opts.reject,
        ignore_repeat=opts.ignore_repeat
    )