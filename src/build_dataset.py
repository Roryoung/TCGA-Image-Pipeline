import os
from optparse import OptionParser
import numpy as np
from random import shuffle
import h5py

from tile import Tile
from normalize import Normalizer
from labeling_util import *
from get_set_data import split_to_sets

def build_dataset(slide_dir, output_dir, background=0.2, size=255, reject_rate=0.1, ignore_repeat=False):
    data = get_projects_info(["TCGA-DLBC"])

    all_cases = list(data['case to images'].keys())
    shuffle(all_cases)

    #split to train and val+test
    train_len = int(0.8*len(all_cases))
    train_set = all_cases[:train_len]
    all_cases = all_cases[train_len:]

    #split val+test into val and test
    val_len = int(0.5*len(all_cases))
    val_set = all_cases[:val_len]
    test_set = all_cases[val_len:]

    train_data = split_to_sets(train_set, data, os.path.join("data", "train"))
    val_data = split_to_sets(val_set, data, os.path.join("data", "val"))
    test_data = split_to_sets(test_set, data, os.path.join("data", "test"))

    dataset = [
        (list(train_data["image to sample"].keys()), "train"),
        (list(val_data["image to sample"].keys()), "val"),
        (list(test_data["image to sample"].keys()), "test")
    ]

    # train_images = ["TCGA-44-7671-01A-01-BS1.914604a2-de9c-404d-9fa5-23fbd0b76da3.svs"]
    # val_images = ["TCGA-FF-8041-01A-01-TS1.b8b69ce3-a325-4864-a5b0-43c450347bc9.svs"]
    # test_images = ["TCGA-G8-6326-01A-01-TS1.e0eb24da-6293-4ecb-8345-b70149c84d1e.svs"]

    # # # train_images = []
    # # val_images = []
    # # test_images = []

    # dataset = [
    #     (train_images, "train"),
    #     (val_images, "val"),
    #     (test_images, "test")
    # ]

    normalizer = Normalizer()
    for images, label in dataset:
        set_path = os.path.join(output_dir, label)

        set_hdf5_path = os.path.join(output_dir, f"{label}.h5")
        set_hdf5_file = h5py.File(set_hdf5_path, 'a')

        for filename in images:
            download_image(filename, slide_dir)
            
            Tile(
                slide_loc=os.path.join(slide_dir, filename),
                set_hdf5_file=set_hdf5_file,
                normalizer=normalizer,
                background=background,
                size=size,
                reject_rate=reject_rate,
                ignore_repeat=ignore_repeat
            )
        
        set_hdf5_file.close()

    normalizer.normalize_dir(output_dir)


if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide_folder> <output_folder> [options]')
    parser.add_option('-b', '--background', dest='background', type='float', default=0.2, help='Percentage of background allowed, defualt=0.2')
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
        reject_rate=opts.reject,
        ignore_repeat=opts.ignore_repeat
    )