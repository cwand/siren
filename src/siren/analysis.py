import numpy as np
import numpy.typing as npt
import scipy
import SimpleITK as sitk


def find_peak(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float = 0.0) -> tuple[float, float]:
    '''
    Find the maximum value in the given TAC.
    :param tac: TAC dictionary
    :param label: Label of interest
    :param start: Start search from the given time
    :return: The time and value of the maximum data point in the TAC.
    '''

    # Filter the data to start at the given start value
    t = tac['tacq'][tac['tacq'] >= start]
    a = tac[label][tac['tacq'] >= start]

    # Find index of largest element
    idx = np.argmax(a)

    return float(t[idx]), float(a[idx])


def find_first_under(tac: dict[str, npt.NDArray[np.float64]],
                     label: str,
                     value: float,
                     start: float = 0.0) -> float:
    '''
    Find the first time point where a label TAC is under a certain value
    :param tac: The TAC dictionary
    :param label: The label of interest
    :param value: The trigger value
    :param start: Start the search from a given start time, ignoring earlier
    times.
    :return: The value of tac['tacq'] at the first point after start where
    tac[label] < value. If that time is not observed in the TAC, -1 is
    returned.
    '''

    a = tac[label]
    t = tac['tacq']
    t_true = t[np.nonzero(np.logical_and(a < value, t >= start))]
    if t_true.size == 0:
        return -1
    else:
        return float(np.min(t_true))


def integrate(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float,
              end: float) -> float:
    '''
    Trapezoid integration of the TAC data
    :param tac: The TAC dictionary
    :param label: The label of interest
    :param start: Start the integration from this point
    :param end: End the integration at this point
    :return: The trapezoid integral of the TAC data. Only data between start
    and end are included.
    '''

    # Filter the data
    t = tac['tacq'][np.logical_and(tac['tacq'] >= start, tac['tacq'] <= end)]
    a = tac[label][np.logical_and(tac['tacq'] >= start, tac['tacq'] <= end)]

    return float(scipy.integrate.trapezoid(a, t))


def avg(tac: dict[str, npt.NDArray[np.float64]],
        label: str,
        start: float) -> float:

    a = tac[label][tac['tacq'] >= start]
    return float(np.mean(a))


def get_volume(roi_path: str) -> float:
    roi = sitk.ReadImage(roi_path)
    binary_mask = sitk.BinaryThreshold(roi,
                                       lowerThreshold=1,
                                       upperThreshold=1)
    n_vox = sitk.GetArrayFromImage(binary_mask).sum()
    spacing = roi.GetSpacing()
    return float(n_vox * spacing[0] * spacing[1] * spacing[2] / 1000)
