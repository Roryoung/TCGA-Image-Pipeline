import os
from optparse import OptionParser
import numpy as np

from tile import Tile
from normalize import Normalizer

def build_dataset(slide_dir, output_dir, background=0.2, size=255, reject_rate=0.1, ignore_repeat=False):
    """
        General function used to build the image dataset. This crops each image at different zoom levels then normalizes the output directory.

        Args:
            - slide_dir: The location where the svs images are stored.
            - output_dir: The location where the set of tiles will be saved.
            - background: Percentage of background allowed. (default = 0.2)
            - size: The width and height of the output tiles. (default = 255)
            - reject_rate: The percentage of rejected tiles kept for sanity check (default = 0.1)
            - ignore_repeat: Automatically overwrite repeated files in the dataset. (default = False)
    """
    normalizer = Normalizer()
    
    for filename in os.listdir(slide_dir):
        slide_path = os.path.join(slide_dir, filename)
        Tile(
            slide_loc=slide_path,
            output_dir=output_dir,
            normalizer=normalizer,
            background=background,
            size=size,
            reject_rate=reject_rate,
            ignore_repeat=ignore_repeat
        )

    normalizer.normalize_dir(output_dir)

if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide_folder> <output_folder> [options]')
    parser.add_option('-b', '--background', dest='background', type='float', default=0.2, help='Percentage of background allowed, default=0.2')
    parser.add_option('-s', '--size', dest='tile_size', type='int', default=255, help='tile size, default=255')
    parser.add_option('-r', '--reject', dest='reject', type='float', default=0.1, help='Precentage of rejected background tiles to save, default=0.1')
    parser.add_option('-i', '--ignore_repeat', dest='ignore_repeat', action="store_true", help='Automatically overwrite repeated files in the dataset, default=False')

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
        reject_rate=opts.reject,
        ignore_repeat=opts.ignore_repeat
    )