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
                  t_peak: float,
                  t_max_left: float,
                  t_max_right: float,
                  t_half_left: float,
                  t_half_right: float,
                  t_func_min: float,
                  t_func_max: float,
                  split_function_left_kidney: float,
                  split_function_right_kidney: float,
                  retention20_left: float,
                  retention20_right: float): ...

def find_peak(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float = ...) -> tuple[float, float]: ...

def find_peak_half(tac: dict[str, npt.NDArray[np.float64]],
                   label: str,
                   start: float = 0.0) -> float: ...

def integrate(tac: dict[str, npt.NDArray[np.float64]],
              label: str,
              start: float,
              end: float) -> float: ...
