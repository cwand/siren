import numpy as np
import numpy.typing as npt

from datetime import datetime

def get_acq_datetime(dicom_path: str) -> datetime: ...

def get_tac_from_paths(series_path: str,
                       roi_paths: list[str],
                       progress: bool = ...) \
        -> dict[str, npt.NDArray[np.float64]]: ...
