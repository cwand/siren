import sys
import importlib.metadata
import time
import argparse


def main(sys_args: list[str]):

    # Get version number from pyproject.toml
    __version__ = importlib.metadata.version("tictac")
    start_time = time.time_ns()

    print("Starting SIREN", __version__)
    print()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", metavar="IMG_DIR", required=True,
                        help="Path to dynamic dicom images.")
    parser.add_argument("-r", metavar="ROI", required=True, nargs="*",
                        help="ROI-image files to analyse.")

    # args = parser.parse_args(sys_args)

    # Report successful end of program
    run_time = (time.time_ns() - start_time) * 1e-9
    print(f'SIREN finished successfully in {run_time:.1f} seconds.')
    print()


if __name__ == "__main__":
    main(sys.argv[1:])
