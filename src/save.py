import os

import h5py
from PIL import Image
import numpy as np

def create_h5_dataset(images_path, zoom, hdf5_file):
    images_path = os.path.join(images_path, zoom)

    num_samples = len(next(os.walk(images_path))[2])
    patch_h = 255
    patch_w = 255
    n_ch = 3

    img_db_shape = (num_samples, patch_h, patch_w, n_ch)
    names_db_shape = (num_samples, 1)

    dt = h5py.special_dtype(vlen=str)
    img_storage = hdf5_file.create_dataset(name='images', shape=img_db_shape, dtype=np.float32)
    name_storage = hdf5_file.create_dataset(name='file_name', shape=names_db_shape, dtype=dt)

    reject_images_path = os.path.join(images_path, "rejected")
    if os.path.exists(reject_images_path):
        reject_num_samples = len(next(os.walk(reject_images_path))[2])
    else:
        reject_num_samples = 0


    reject_img_db_shape = (reject_num_samples, patch_h, patch_w, n_ch)
    reject_names_db_shape = (reject_num_samples, 1)

    reject_img_storage = hdf5_file.create_dataset(name="reject_images", shape=reject_img_db_shape, dtype=np.float32)
    reject_name_storage = hdf5_file.create_dataset(name="reject_file_name", shape=reject_names_db_shape, dtype=dt)

    for index_img, entry in enumerate(os.scandir(images_path)):
        if not entry.is_dir():
            img = Image.open(images_path +'/' + entry.name)
            tissue_img = np.array(img)

            img_storage[index_img, :, :, :] = tissue_img
            name_storage[index_img, :] = entry.name
        elif entry.name == "rejected":
            for reject_index_img, reject_entry in enumerate(os.scandir(entry)):
                img = Image.open(os.path.join(entry, reject_entry.name))
                tissue_img = np.array(img)

                reject_img_storage[reject_index_img, :, :, :] = tissue_img
                reject_name_storage[reject_index_img, :] = entry.name


def crate_h5_set(output, set_name):
    set_path = os.path.join(output, set_name)

    set_hdf5_path = os.path.join(output, f"{set_name}.h5")
    set_hdf5_file = h5py.File(set_hdf5_path, 'w')

    for case in os.scandir(set_path):
        patient_hdf5 = set_hdf5_file.create_group(case.name)
        for zoom in os.scandir(case):
            zoom_hdf5 = patient_hdf5.create_group(zoom.name)
            create_h5_dataset(case, zoom.name, zoom_hdf5)

    set_hdf5_file.close()


def h5_output(output):
    for set_name in os.scandir(output):
        crate_h5_set(output, set_name.name)

