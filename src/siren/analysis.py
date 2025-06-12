import numpy as np
import numpy.typing as npt
import scipy


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


def find_peak_half(tac: dict[str, npt.NDArray[np.float64]],
                   label: str,
                   start: float = 0.0) -> float:
    '''
    Find the time for the TAC to decay to half of the maximum value
    :param tac: The TAC dictionary
    :param label: The label of interest
    :param start: The max value is taken as the largest value after "start"
    :return: Return the time point where the TAC first dives under half the
    maximum value. If the data never goes under the half value a negative value
    is returned.
    '''

    # Filter out elements before start
    t = tac['tacq'][tac['tacq'] >= start]
    a = tac[label][tac['tacq'] >= start]

    # Find index of the largest value
    idx = np.argmax(a)
    # Save the maximum value
    amax = a[idx]

    # Iterate through the data and find the first element under the half value
    while a[idx] > 0.5 * amax:
        idx = idx + 1
        if idx == len(a):
            # We have iterated to the end, return an error code.
            return -1.0

    return float(t[idx])


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
