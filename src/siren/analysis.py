import numpy as np
import numpy.typing as npt
import scipy


def find_peak(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float = 0.0) -> tuple[float, float]:

    t = tac['tacq'][tac['tacq'] >= start]
    a = tac[label][tac['tacq'] >= start]

    idx = np.argmax(a)

    return float(t[idx]), float(a[idx])


def find_peak_half(tac: dict[str, npt.NDArray[np.float64]],
                   label: str,
                   start: float = 0.0) -> float:

    t = tac['tacq'][tac['tacq'] >= start]
    a = tac[label][tac['tacq'] >= start]

    idx = np.argmax(a)
    amax = a[idx]
    while a[idx] > 0.5 * amax:
        idx = idx + 1
        if idx == len(a):
            return -1.0

    return float(t[idx])


def integrate(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float,
              end: float) -> float:

    t = tac['tacq'][np.logical_and(tac['tacq'] >= start, tac['tacq'] <= end)]
    a = tac[label][np.logical_and(tac['tacq'] >= start, tac['tacq'] <= end)]

    return float(scipy.integrate.trapezoid(a, t))
