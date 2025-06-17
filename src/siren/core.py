import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import numpy.typing as npt

from datetime import datetime

import SimpleITK as sitk
from matplotlib.gridspec import GridSpec

from tqdm import tqdm


def get_acq_datetime(dicom_path: str) -> datetime:
    '''
    Get the acquisition time and date form the dicom metadata of an image
    :param dicom_path: The path to the dicom file
    :return: A datetime object reflecting the acquisition time
    '''

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
    '''
    Extract TACs from an image series and a list of roi files.
    :param series_path: The path to the dynamic image data (dicom).
    :param roi_paths: A list of all roi files to extract
    :param progress: Whether or not to show a progress bar.
    :return: Returns a dictionary object. The keys will be the filenames
    passed to this function, and the respective values till be the TAC for that
    ROI in a numpy array. Furthermore the key 'tacq' will contain the time data
    (in seconds).
    '''

    # Prepare series reader
    reader = sitk.ImageSeriesReader()

    # Get dicom file names in folder sorted according to acquisition time.
    dcm_names = reader.GetGDCMSeriesFileNames(series_path)

    # Read in all rois and filenames
    rois = [(sitk.ReadImage(roi), roi) for roi in roi_paths]

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


def make_renogram(tac: dict[str, npt.NDArray[np.float64]],
                  left_kidney: str,
                  right_kidney: str,
                  aorta: str,
                  res_dict: dict[str, float]):
    '''
    Draw a renogram on the screen
    '''

    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(2, 1, height_ratios=[3, 1])

    a_max = np.max([np.max(tac[left_kidney]), np.max(tac[right_kidney])])
    aorta_max = np.max(tac[aorta])

    axs = fig.add_subplot(gs[0])
    axs.plot(tac['tacq'] / 60, tac[left_kidney]/a_max, 'g-',
             label='kidney_left')
    axs.plot(tac['tacq'] / 60, tac[right_kidney]/a_max, 'r-',
             label='kidney_right')
    axs.plot(tac['tacq'] / 60, 0.1 * tac[aorta] / aorta_max, 'k-',
             label='aorta')
    axs.axvline(res_dict['t_peak'] / 60, linestyle='--', color='b',
                label='Aorta peak')
    axs.axvline(res_dict['tmax_left'] / 60, linestyle='--', color='g')
    axs.axvline(res_dict['tmax_right'] / 60, linestyle='--', color='r')
    # if res_dict['thalf_left'] > 0:
    #     axs.axvline(res_dict['thalf_left'] / 60, linestyle='--', color='g')
    # if res_dict['thalf_right'] > 0:
    #     axs.axvline(res_dict['thalf_right'] / 60, linestyle='--', color='r')
    axs.axvline(res_dict['tfunc_min'] / 60, linestyle='--', color='k',
                label='Split (perf/func)')
    axs.axvline(res_dict['tfunc_max'] / 60, linestyle='--', color='k')
    axs.axvline(0.0, linestyle='--', color='k')
    axs.axvline(res_dict['tperf_max'] / 60, linestyle='--', color='k')

    axs.set_xlabel('Time [minutes]')
    axs.set_ylabel('Activity (normalised)')

    axs.legend()
    axs.set_xticks(np.arange(0.0, float(np.max(tac['tacq']))/60, 1.0))
    axs.set_yticks(np.arange(0.0, 1.0, 0.1))
    axs.grid()

    df = pd.DataFrame({
        "Left": [f'{res_dict["vol_left"]:.1f}',
                 f'{res_dict["split_left"]:.1f}',
                 f'{res_dict["perf_left"]:.1f}',
                 f'{(res_dict["tmax_left"] - res_dict["t_peak"]) / 60:.1f}',
                 f'{res_dict["ret20_left"]:.1f}'],
        "Right": [f'{res_dict["vol_right"]:.1f}',
                  f'{res_dict["split_right"]:.1f}',
                  f'{res_dict["perf_right"]:.1f}',
                  f'{(res_dict["tmax_right"] - res_dict["t_peak"]) / 60:.1f}',
                  f'{res_dict["ret20_right"]:.1f}']
    })

    axs = fig.add_subplot(gs[1])
    axs.table(cellText=df.values, loc='center',
              colLabels=df.columns,
              rowLabels=['Volume [cm^3]',
                         'Split function [%]',
                         'Perfusion split [%]',
                         'Tmax [min]',
                         'Retention @ 20 min [%]'],
              bbox=[0.3, 0, 0.6, 1])  # type: ignore
    axs.axis('off')

    plt.show()
