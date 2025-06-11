import os.path

import numpy as np
import numpy.typing as npt

from datetime import datetime

import SimpleITK as sitk

from tqdm import tqdm


def get_acq_datetime(dicom_path: str) -> datetime:
    """Get an image acquisition datetime from its dicom header.
    Dicom images store the acquisition date and time in tags in the images
    dicom header. This function reads the relevant tags and turns it into a
    datetime object.

    Arguments:
    dicom_path  --  The path to the dicom file.

    Return value:
    A datetime object representing the date and time of the acquisition.
    """

    # Read the dicom image into Simple ITK
    img = sitk.ReadImage(dicom_path)

    # Read the relevant header tags as strings
    img_time = img.GetMetaData('0008|0032')
    img_date = img.GetMetaData('0008|0022')

    # Format the strings into ISO 8601 format [ YYYY-MM-DD hh:mm:ss.ffffff ]
    sd = img_date[:4] + "-" + img_date[4:6] + "-" + img_date[6:]
    sd = sd + " " + img_time[:2] + ":" + img_time[2:4] + ":" + img_time[4:6]
    sd = sd + "." + img_time[-1].ljust(6, "0")
    return datetime.fromisoformat(sd)


def get_tac_from_paths(series_path: str,
                       roi_paths: list[str],
                       progress: bool = True) \
        -> dict[str, npt.NDArray[np.float64]]:

    # Prepare series reader
    reader = sitk.ImageSeriesReader()

    # Get dicom file names in folder sorted according to acquisition time.
    dcm_names = reader.GetGDCMSeriesFileNames(series_path)

    # Read in all rois and filenames
    rois = [(sitk.ReadImage(roi), os.path.basename(roi)) for roi in roi_paths]

    # Prepare label statistics filter
    label_stats_filter = sitk.LabelStatisticsImageFilter()

    # Get acquisition time of first image
    acq0 = get_acq_datetime(dcm_names[0])

    # Initialise container for results
    res: dict[str, npt.NDArray[np.float64]] = {
        'tacq': np.zeros(len(dcm_names))
    }
    for roi in rois:
        res[roi[1]] = np.zeros(len(dcm_names))

    for i, name in enumerate(tqdm(dcm_names, disable=(not progress))):

        # Load images in order
        img = sitk.ReadImage(name)

        # Find acquisition time and store in list
        res['tacq'][i] = (get_acq_datetime(name) - acq0).total_seconds()

        for roi in rois:

            # Apply label stats filter on original img and read ROI means
            label_stats_filter.Execute(img, roi[0])

            # Append the roi sum value to the list for each label.
            res[roi[1]][i] = label_stats_filter.GetSum(1)

    return res
