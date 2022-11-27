__all__ = ['bin_frame', 'central_crop']

import numpy as np


def bin_frame(data_frame, binning):

    binning = int(binning)

    if binning <= 1:
        return data_frame

    new_frame = []
    for xx in range(len(data_frame))[::binning]:
        new_frame.append(list(np.sum(data_frame[xx:xx + binning], 0)))

    new_frame = np.swapaxes(np.array(new_frame), 0, 1)

    new_frame2 = []
    for xx in range(len(new_frame))[::binning]:
        new_frame2.append(list(np.sum(new_frame[xx:xx + binning], 0)))

    data_frame = np.swapaxes(np.array(new_frame2), 0, 1)

    return data_frame


def central_crop(original_array, destination_fits):

    crop1 = len(original_array) / 2 - len(destination_fits[1].data) / 2
    crop2 = len(original_array) / 2 + len(destination_fits[1].data) / 2

    return original_array[crop1:crop2, crop1:crop2]
