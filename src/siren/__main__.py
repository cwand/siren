import os
import sys
import importlib.metadata
import argparse

import numpy as np

import siren


def main(sys_args: list[str]):

    # Get version number from pyproject.toml
    __version__ = importlib.metadata.version("siren")

    print("Starting SIREN", __version__)
    print()
    print(f'Current working directory: {os.getcwd()}')
    print()

    # Parse system arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", metavar="IMG_DIR", required=True,
                        help="Path to dynamic dicom images.")
    parser.add_argument("--left_kidney", metavar="ROI_LEFT_KIDNEY",
                        help="ROI-image file for the left kidney.")
    parser.add_argument("--right_kidney", metavar="ROI_RIGHT_KIDNEY",
                        help="ROI-image file for the right kidney.")
    parser.add_argument("--aorta", metavar="ROI_AORTA",
                        help="ROI-image file for the aorta input function.")

    args = parser.parse_args(sys_args)

    # Load data and extract TACs
    print(f'Image data: {args.i}')
    print('ROI data:')
    print(f' Aorta: {args.aorta}')
    print(f' Left kidney: {args.left_kidney}')
    print(f' Right kidney: {args.right_kidney}')

    roi_paths = [
        args.aorta,
        args.left_kidney,
        args.right_kidney
    ]

    print('Extracting TACs from images...')
    tac = siren.get_tac_from_paths(
        series_path=args.i,
        roi_paths=roi_paths
    )
    print()

    # Start analysis
    print('Finding aorta peak...')
    t_peak, a_peak = siren.find_peak(tac, args.aorta)
    print(f'Input peak time: {t_peak:.1f} seconds')
    print()

    print('Finding Tmax...')
    t_max_left, a_max_left = siren.find_peak(tac,
                                             args.left_kidney,
                                             start=t_peak + 60)
    t_max_right, a_max_right = siren.find_peak(tac,
                                               args.right_kidney,
                                               start=t_peak + 60)
    print('Tmax:')
    print(f' Left kidney:  {(t_max_left-t_peak)/60:.1f} min')
    print(f' Right kidney: {(t_max_right-t_peak)/60:.1f} min')
    print()

    print('Computing split function...')
    t_func_min = t_peak + 60
    t_func_max = t_peak + 150
    print(f'Split function will be evaluated between {t_func_min/60:.1f} and '
          f'{t_func_max/60:.1f} mins.')
    func_left = siren.integrate(tac, args.left_kidney,
                                t_func_min, t_func_max)
    func_right = siren.integrate(tac, args.right_kidney,
                                 t_func_min, t_func_max)
    sf_left = 100 * func_left/(func_left + func_right)
    sf_right = 100 * func_right/(func_left + func_right)
    print(f' Left kidney: {sf_left:.1f} %')
    print(f' Right kidney: {sf_right:.1f} %')
    print()

    print('Retention at T=20 min (% of uptake at Tmax)')
    idx = np.searchsorted(tac['tacq'], 20 * 60 + t_peak, side='right')
    if idx == len(tac['tacq']):
        print("Renogram ended before T = 20 min. Using final time point.")
        idx = idx - 1
    ret20_left = 100 * tac[args.left_kidney][idx] / a_max_left
    ret20_right = 100 * tac[args.right_kidney][idx] / a_max_right
    print(f' Left: {ret20_left:.1f} %')
    print(f' Right: {ret20_right:.1f} %')
    print()

    print('Computing T½:')
    t_half_left = siren.find_peak_half(tac,
                                       args.left_kidney,
                                       start=t_peak + 60)
    if t_half_left < 0:
        print("Could not find T½ for left kidney!")
    else:
        print(f' Left: {(t_half_left - t_peak)/60:.1f} min')
    t_half_right = siren.find_peak_half(tac,
                                        args.right_kidney,
                                        start=t_peak + 60)
    if t_half_right < 0:
        print("Could not find T½ for right kidney!")
    else:
        print(f' Right: {(t_half_right - t_peak)/60:.1f} min')
    print()

    siren.make_renogram(tac=tac,
                        left_kidney=args.left_kidney,
                        right_kidney=args.right_kidney,
                        t_peak=t_peak,
                        t_max_left=t_max_left,
                        t_max_right=t_max_right,
                        t_half_left=t_half_left,
                        t_half_right=t_half_right,
                        t_func_min=t_func_min,
                        t_func_max=t_func_max,
                        split_function_left_kidney=sf_left,
                        split_function_right_kidney=sf_right,
                        retention20_left=ret20_left,
                        retention20_right=ret20_right)

    # Report successful end of program
    print('SIREN finished successfully.')
    print()


if __name__ == "__main__":
    main(sys.argv[1:])
