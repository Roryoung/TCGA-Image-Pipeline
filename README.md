# TCGA-Image-Pipeline #

## Setup ##
  - System Packages
    - `openslide`


  - Python packages
    - `pip3 install -r requirements.txt`

## Usage ##

### Tile Creation ###
    Usage: tile.py <slide> <output_folder> [options]

    Options:
    -h, --help            show this help message and exit
    -b BACKGROUND, --background=BACKGROUND
                            Percentage of background allowed, default=0.2
    -s TILE_SIZE, --size=TILE_SIZE
                            Size of the output tiles, default=255
    -r REJECT, --reject=REJECT
                            Precentage of rejected background tiles to save,
                            default=0.1
    -i, --ignore_repeat   Automatically overwrte repeated files in the dataset,
                            default=False

### Stain Normalization ###
    Usage: normalize.py <tile_dir>

    Options:
    -h, --help  show this help message and exit

### Build Dataset ###
    Usage: build_dataset.py <slide_folder> <output_folder> [options]

    Options:
    -h, --help            show this help message and exit
    -b BACKGROUND, --background=BACKGROUND
                            Percentage of background allowed, defualt=0.2
    -s TILE_SIZE, --size=TILE_SIZE
                            tile size, defualt=255
    -r REJECT, --reject=REJECT
                            Precentage of rejected background tiles to save,
                            defualt=0.1
    -i, --ignore_repeat   Automatically overwrte repeated files in the dataset,
                            defualt=False