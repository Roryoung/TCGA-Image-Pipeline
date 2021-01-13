import os

import h5py
import pandas as pd

from labeling_util import store_hugo, load_hugo

def recursive_save_to_h5(h5_file, path, item):
    if isinstance(item, dict):
        for key, value in item.items():
            recursive_save_to_h5(h5_file, f"{path}{key}/", value)
    elif isinstance(item, list):
        if len(item) > 0:
            for key, value in enumerate(item):
                recursive_save_to_h5(h5_file, f"{path}{key}/", value)
        else:
            h5_file.create_group(f"{path}/")
    elif item is None:
        h5_file[path] = h5py.Empty("f")
    else:
        h5_file[path] = item

    h5_file[path].attrs["type"] = type(item).__name__


def get_data_dict(case_set, data, h5_file_name):
    set_data_dict = {}
    for case in case_set:
        for key, value in data["data dict"].items():
            if case in value:
                if key not in set_data_dict:
                    set_data_dict[key] = {}

                set_data_dict[key][case] = value[case]

    with h5py.File(h5_file_name, "a") as h5_file:
        recursive_save_to_h5(h5_file, "data_dict/", set_data_dict)

    return set_data_dict


def get_case_to_images(case_set, data, h5_file_name=None):
    set_case_to_images = {}
    for case in case_set:
        for key, value in data["case to images"].items():
            if case == key:
                set_case_to_images[key] = value
    
    if h5_file_name is not None:
        with h5py.File(h5_file_name, "a") as h5_file:
            recursive_save_to_h5(h5_file, "case_to_images/", set_case_to_images)

    return set_case_to_images


def get_image_to_sample(case_set, data, h5_file_name=None):
    set_image_to_sample = {}
    for key, value in get_case_to_images(case_set, data).items():
        for image in value:
            set_image_to_sample[image] = key
    
    if h5_file_name is not None:
        with h5py.File(h5_file_name, "a") as h5_file:
            recursive_save_to_h5(h5_file, "image_to_sample/", set_image_to_sample)

    return set_image_to_sample
        

def get_mutational_signatures(case_set, data, h5_file_name):
    mutational_signatures = data["mutational signatures"][data["mutational signatures"]["case_id"].isin(case_set)]
    mutational_signatures.to_hdf(h5_file_name, "/mutational_signatures", format="table")

    return mutational_signatures


def get_labels(case_set, data, h5_file_name):
    sample_ids = list(map( lambda x : data['image to sample'][x], list(get_image_to_sample(case_set, data).keys())))
    labels = data['labels'][data['labels']['sample.barcode'].isin(sample_ids)]
    labels.to_hdf(h5_file_name, "/labels", format="table")

    return labels


def get_hugo_symbols(case_set, data, h5_file_name):
    hugo_symbols = data["hugo symbols"][data["hugo symbols"]["case_barcode"].isin(case_set)]
    with h5py.File(h5_file_name, "a") as h5_file:
        store_hugo(h5_file,hugo_symbols,overwrite=True)

    return hugo_symbols


def split_to_sets(case_set, data, h5_file_name):
    return {
        "data dict": get_data_dict(case_set, data, h5_file_name),
        "image to sample": get_image_to_sample(case_set, data, h5_file_name),
        "case to images":  get_case_to_images(case_set, data, h5_file_name),
        "labels": get_labels(case_set, data, h5_file_name),
        "mutational signatures": get_mutational_signatures(case_set, data, h5_file_name),
        "hugo symbols": get_hugo_symbols(case_set, data, h5_file_name)
    }

def recursive_load_from_h5(h5_file, path):
    if h5_file[path].attrs["type"] == dict.__name__:
        return_dict = {}
        for key, value in h5_file[path].items():
            return_dict[key] = recursive_load_from_h5(h5_file, f"{path}/{key}")
        return return_dict

    elif h5_file[path].attrs["type"] == list.__name__:
        return_list = []
        for key, value in h5_file[path].items():
            return_list.append(recursive_load_from_h5(h5_file, f"{path}/{key}"))
        
        return return_list
    elif h5_file[path] is h5py.Empty("f"):
        return None
    else:
        return h5_file[path][()]

def load_set_data(h5_file_loc):
    with h5py.File(h5_file_loc, "r") as h5_file:
        return {
            "data dict": recursive_load_from_h5(h5_file, "data_dict"),
            "image to sample": recursive_load_from_h5(h5_file, "image_to_sample/"),
            "case to images":  recursive_load_from_h5(h5_file, "case_to_images/"),
            "labels": pd.read_hdf(h5_file_loc, key="labels"),
            "mutational signatures": pd.read_hdf(h5_file_loc, key="mutational_signatures"),
            "hugo symbols": load_hugo(h5_file)
        }
