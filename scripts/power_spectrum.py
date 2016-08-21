#!/home/meteoboy4/anaconda/envs/Meteorology/bin/python

import sys
import argparse
import numpy as np
import xarray as xr

from detection import fft


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Power spectrum analysis using fft along both latitude and longitude")
    parser.add_argument("input", help="the path to the netCDF file to be analyzed")
    parser.add_argument("output", help="the output directory to which the plots go")

    args = parser.parse_args()

    target_file = xr.open_dataset(args.input)
    variable = target_file[target_file.keys()[-1]]
    variable_slice = variable.sel(longitude=slice(120, 180), latitude=slice(60, 30), time=slice('2005-01-08', '2005-01-11'))

    latitude_range = range(30, 60, 1)
    for i, lat in enumerate(latitude_range):
        fft2 = fft.FFT2(variable_slice.sel(latitude=lat), latitude=lat, folder=args.output)
        if i == 0:
            power = fft2.shiftized
        else:
            power += fft2.shiftized
        fft2.draw()
    power /= len(latitude_range)
    fft.draw(target=power, folder=args.output, latitude=lat, plot_name='{0}-{1}mean'.format(latitude_range[0], latitude_range[-1]), xcood=fft2.xcood, ycood=fft2.ycood)

    longitude_range = range(120, 160, 1)
    for i, lon in enumerate(longitude_range):
        fft2 = fft.FFT2(variable_slice.sel(longitude=lon), longitude=lon, folder=args.output)
        if i == 0:
            power = fft2.shiftized
        else:
            power += fft2.shiftized
        fft2.draw()
    power /= len(longitude_range)
    fft.draw(target=power, folder=args.output, longitude=lon, plot_name='{0}-{1}mean'.format(longitude_range[0], longitude_range[-1]), xcood=fft2.xcood, ycood=fft2.ycood)

if __name__ == '__main__':
    main()
