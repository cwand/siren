import numpy as np
import numpy.typing as npt

from datetime import datetime

def get_acq_datetime(dicom_path: str) -> datetime: ...

def get_tac_from_paths(series_path: str,
                       roi_paths: list[str],
                       progress: bool = ...) \
        -> dict[str, npt.NDArray[np.float64]]: ...

def make_renogram(tac: dict[str, npt.NDArray[np.float64]],
                  left_kidney: str,
                  right_kidney: str,
                  aorta: str,
                  res_dict: dict[str, float]): ...

def find_peak(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float = ...) -> tuple[float, float]: ...

def find_first_under(tac: dict[str, npt.NDArray[np.float64]],
                     label: str,
                     value: float,
                     start: float = 0.0) -> float: ...

def integrate(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float,
              end: float) -> float: ...

def avg(tac: dict[str, npt.NDArray[np.float64]],
        label: str,
        start: float) -> float: ...

def get_volume(roi_path: str) -> float: ...
