import os
import json
from json import JSONEncoder

import numpy as np


class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        return JSONEncoder.default(self, obj)


def get_data_dict(case_set, data, output):
    set_data_dict = {}
    for case in case_set:
        for key, value in data["data dict"].items():
            if case in value:
                if key not in set_data_dict:
                    set_data_dict[key] = {}

                set_data_dict[key][case] = value[case]

    with open(os.path.join(output, "data_dict.json"), "w") as fp:
        json.dump(set_data_dict, fp, cls=NumpyArrayEncoder)

    return set_data_dict


def get_case_to_images(case_set, data, output=None):
    set_case_to_images = {}
    for case in case_set:
        for key, value in data["case to images"].items():
            if case == key:
                set_case_to_images[key] = value

    if output is not None:
        with open(os.path.join(output, "case_to_images.json"), "w") as fp:
            json.dump(set_case_to_images, fp, cls=NumpyArrayEncoder)

    return set_case_to_images


def get_image_to_sample(case_set, data, output=None):
    set_image_to_sample = {}
    for key, value in get_case_to_images(case_set, data).items():
        for image in value:
            set_image_to_sample[image] = key
    
    if output is not None:
        with open(os.path.join(output, "image_to_sample.json"), "w") as fp:
            json.dump(set_image_to_sample, fp, cls=NumpyArrayEncoder)

    return set_image_to_sample
        

def get_mutational_signatures(case_set, data, output):
    mutational_signatures = data["mutational signatures"][data["mutational signatures"]["case_id"].isin(case_set)]
    mutational_signatures.to_csv(os.path.join(output, "mutational_signatures.csv"))

    return mutational_signatures


def get_labels(case_set, data, output):
    sample_ids = list(map( lambda x : data['image to sample'][x], list(get_image_to_sample(case_set, data).keys())))
    labels = data['labels'][data['labels']['sample.barcode'].isin(sample_ids)]
    labels.to_csv(os.path.join(output, "labels.csv"))

    return labels


def get_hugo_symbols(case_set, data, output):
    hugo_symbols = data["hugo symbols"][data["hugo symbols"]["case_barcode"].isin(case_set)]
    hugo_symbols.to_csv(os.path.join(output, "hugo_symbols.csv"))

    return hugo_symbols


def split_to_sets(case_set, data, output):
    if not os.path.exists(output):
        os.makedirs(output)

    return {
        "data dict": get_data_dict(case_set, data, output),
        "image to sample": get_image_to_sample(case_set, data, output),
        "case to images":  get_case_to_images(case_set, data, output),
        "labels": get_labels(case_set, data, output),
        "mutational signatures": get_mutational_signatures(case_set, data, output),
        "hugo symbols": get_hugo_symbols(case_set, data, output)
    }
