import os
from optparse import OptionParser
import numpy as np
from random import shuffle
import h5py

from tile import Tile
from normalize import Normalizer
from labeling_util import *
from get_set_data import split_to_sets, load_set_data

def build_dataset(slide_dir, output_dir, projects, background=0.2, size=255, reject_rate=0.1, ignore_repeat=False):
    proceed = None
    train_path = os.path.join(output_dir, "train.h5")
    val_path = os.path.join(output_dir, "val.h5")
    test_path = os.path.join(output_dir, "test.h5")

    if (os.path.isfile(train_path) and os.path.isfile(val_path) and os.path.isfile(test_path)):
        while not (proceed == "C" or proceed == "A" or proceed == "R" or proceed == "Q"):
            print(
                """A dataset already exists in this directory. Do you want to \n
                    - Continue to build the datset [C] \n
                    - Reset the dataset [R] \n
                    - Quit [Q]
                """
            )
            proceed = input().upper()
            if proceed == "R":
                os.remove(train_path)
                os.remove(val_path)
                os.remove(test_path)

    if proceed == "C":
        train_data = load_set_data(train_path)
        val_data = load_set_data(val_path)
        test_data = load_set_data(test_path)

        train_h5 = h5py.File(train_path, 'a')
        val_h5 = h5py.File(val_path, 'a')
        test_h5 = h5py.File(test_path, 'a')
    elif proceed == "R" or proceed == None:
        if projects is None:
            raise ValueError("Missing list of projects to download.")
        data = get_projects_info(projects)

        train_h5 = h5py.File(train_path, 'a')
        val_h5 = h5py.File(val_path, 'a')
        test_h5 = h5py.File(test_path, 'a')

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

        train_data = split_to_sets(train_set, data, train_path)
        val_data = split_to_sets(val_set, data, val_path)
        test_data = split_to_sets(test_set, data, test_path)

    if proceed != "Q":
        dataset = [
            (list(train_data["image to sample"].keys()), train_h5),
            (list(val_data["image to sample"].keys()), val_h5),
            (list(test_data["image to sample"].keys()), test_h5)
        ]

        # train_images = ["TCGA-44-7671-01A-01-BS1.914604a2-de9c-404d-9fa5-23fbd0b76da3.svs"]
        # val_images = ["TCGA-FF-8041-01A-01-TS1.b8b69ce3-a325-4864-a5b0-43c450347bc9.svs"]
        # test_images = ["TCGA-G8-6326-01A-01-TS1.e0eb24da-6293-4ecb-8345-b70149c84d1e.svs"]

        # # # train_images = []
        # val_images = []
        # test_images = []

        # dataset = [
        #     (train_images, train_h5),
        #     (val_images, val_h5),
        #     (test_images, test_h5)
        # ]

        normalizer = Normalizer()
        for images, h5_file in dataset:
            image_h5_file = h5_file.require_group("images")

            for filename in images:
                if proceed != "C" or ".".join(filename.split(".")[:-1]) not in image_h5_file:
                    download_image(filename, slide_dir)
                    
                    Tile(
                        slide_loc=os.path.join(slide_dir, filename),
                        set_hdf5_file=image_h5_file,
                        normalizer=normalizer,
                        background=background,
                        size=size,
                        reject_rate=reject_rate,
                        ignore_repeat=ignore_repeat
                    )
            
            h5_file.close()
            
        normalizer.normalize_dir(output_dir)

def list_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

if __name__ == "__main__":
    parser = OptionParser(usage='Usage: %prog <slide_folder> <output_folder> -p <projects> [options]')
    parser.add_option('-p', '--projects', dest="projects", type='string', action='callback', help="List of TCGA cases e.g. TCGA-LUAD,TCGA-BRCA", callback=list_callback)
    parser.add_option('-b', '--background', dest='background', type='float', default=0.2, help='Percentage of background allowed, defualt=0.2')
    parser.add_option('-s', '--size', dest='tile_size', type='int', default=255, help='tile size, defualt=255')
    parser.add_option('-r', '--reject', dest='reject', type='float', default=0.1, help='Precentage of rejected background tiles to save, defualt=0.1')
    parser.add_option('-i', '--ignore_repeat', dest='ignore_repeat', action="store_true", help='Automatically overwrte repeated files in the dataset, defualt=False')

    (opts, args) = parser.parse_args()

    try:
        slide_dir = args[0]
    except IndexError:
        parser.error('Missing slide image directory argument.')

    try:
        output_dir = args[1]
    except IndexError:
        parser.error('Missing output directory argument.')

    # if opts.projects is None:
    #     raise parser.error("Missing list of projects to download.")

    build_dataset(
        slide_dir=slide_dir,
        output_dir=output_dir,
        projects=opts.projects,
        background=opts.background,
        size=opts.tile_size,
        reject_rate=opts.reject,
        ignore_repeat=opts.ignore_repeat
    )