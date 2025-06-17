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

    res_dict: dict[str, float] = {}

    print('Finding kidney volumes...')
    vol_left = siren.get_volume(args.left_kidney)
    vol_right = siren.get_volume(args.right_kidney)
    print('Volumes:')
    print(f' Left:  {vol_left:.1f} cm^3')
    print(f' right: {vol_right:.1f} cm^3')
    res_dict['vol_left'] = vol_left
    res_dict['vol_right'] = vol_right
    print()

    print('Finding aorta peak...')
    t_peak, a_peak = siren.find_peak(tac, args.aorta)
    print(f'Input peak time: {t_peak:.1f} seconds')
    res_dict['t_peak'] = t_peak
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
    res_dict['tmax_left'] = t_max_left
    res_dict['tmax_right'] = t_max_right
    print()

    print('Computing perfusion split...')
    # t_perf_max = siren.find_first_under(tac,
    #                                     label=args.aorta,
    #                                     value=siren.avg(tac,
    #                                                     label=args.aorta,
    #                                                     start=t_peak + 60),
    #                                     start=t_peak)
    t_perf_max = 2 * t_peak
    print(f'Perfusion split will be evaluated between 0 and'
          f' {t_perf_max:.0f} seconds.')
    res_dict['tperf_max'] = t_perf_max
    perf_left = siren.integrate(tac, args.left_kidney, 0, t_perf_max)
    perf_right = siren.integrate(tac, args.right_kidney, 0, t_perf_max)
    ps_left = 100 * perf_left / (perf_left + perf_right)
    ps_right = 100 * perf_right / (perf_left + perf_right)
    print(f' Left kidney: {ps_left:.1f} %')
    print(f' Right kidney: {ps_right:.1f} %')
    res_dict['perf_left'] = ps_left
    res_dict['perf_right'] = ps_right
    print()

    print('Computing split function...')
    t_func_min = t_peak + 60
    t_func_max = t_peak + 150
    print(f'Split function will be evaluated between {t_func_min/60:.1f} and '
          f'{t_func_max/60:.1f} mins.')
    res_dict['tfunc_min'] = t_func_min
    res_dict['tfunc_max'] = t_func_max
    func_left = siren.integrate(tac, args.left_kidney,
                                t_func_min, t_func_max)
    func_right = siren.integrate(tac, args.right_kidney,
                                 t_func_min, t_func_max)
    sf_left = 100 * func_left/(func_left + func_right)
    sf_right = 100 * func_right/(func_left + func_right)
    print(f' Left kidney: {sf_left:.1f} %')
    print(f' Right kidney: {sf_right:.1f} %')
    res_dict['split_left'] = sf_left
    res_dict['split_right'] = sf_right
    print()

    print('Retention at T=20 min (% of uptake at Tmax)')
    idx = np.searchsorted(tac['tacq'], 20 * 60 + t_peak, side='right')
    if idx == len(tac['tacq']):
        print("Renogram ended before T = 20 min. Using final time point.")
        idx = idx - 1
    ret20_left = float(100 * tac[args.left_kidney][idx] / a_max_left)
    ret20_right = float(100 * tac[args.right_kidney][idx] / a_max_right)
    print(f' Left: {ret20_left:.1f} %')
    print(f' Right: {ret20_right:.1f} %')
    res_dict['ret20_left'] = ret20_left
    res_dict['ret20_right'] = ret20_right
    print()

    print('Computing T½:')
    t_half_left = siren.find_first_under(tac,
                                         label=args.left_kidney,
                                         value=a_max_left * 0.5,
                                         start=t_max_left)
    if t_half_left < 0:
        print("Could not find T½ for left kidney!")
    else:
        print(f' Left: {(t_half_left - t_peak)/60:.1f} min')
    t_half_right = siren.find_first_under(tac,
                                          label=args.right_kidney,
                                          value=a_max_right * 0.5,
                                          start=t_max_right)
    if t_half_right < 0:
        print("Could not find T½ for right kidney!")
    else:
        print(f' Right: {(t_half_right - t_peak)/60:.1f} min')
    res_dict['thalf_left'] = t_half_left
    res_dict['thalf_right'] = t_half_right
    print()

    siren.make_renogram(tac=tac,
                        left_kidney=args.left_kidney,
                        right_kidney=args.right_kidney,
                        aorta=args.aorta,
                        res_dict=res_dict)

    # Report successful end of program
    print('SIREN finished successfully.')
    print()


if __name__ == "__main__":
    main(sys.argv[1:])
