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
    parser.add_argument("--latitude", nargs=2, type=float, default=[60., 30.],
                        help="the latitude range for analysis, must be two numbers indicting the maximum and minimum")
    parser.add_argument("--longitude", nargs=2, type=float, default=[120., 180.],
                        help="the longitude range for analysis, must be two numbers indicting the minimum and maximum")
    parser.add_argument("--level", nargs=3, type=float,
                        help="the contour levels you wish to draw, indicating the minimum, maximum and spacing respectively")

    args = parser.parse_args()
    north_latitude = args.latitude[0]
    south_latitude = args.latitude[1]
    west_longitude = args.longitude[0]
    east_longitude = args.longitude[1]
    levels = np.linspace(args.level[0], args.level[1], args.level[2])
    folder_name = args.output + '_lat{0}-{1}_lon{2}-{3}'.format(south_latitude, north_latitude, west_longitude,
                                                                east_longitude)

    target_file = xr.open_dataset(args.input)
    variable = target_file[target_file.keys()[-1]]
    variable_slice = variable.sel(longitude=slice(west_longitude, east_longitude),
                                  latitude=slice(north_latitude, south_latitude),
                                  time=slice('2005-01-08', '2005-01-11'))

    latitude_range = np.arange(south_latitude, north_latitude, 1)
    for i, lat in enumerate(latitude_range):
        fft2 = fft.FFT2(variable_slice.sel(latitude=lat), latitude=lat, folder=folder_name, levels=levels)
        if i == 0:
            power = fft2.shiftized
        else:
            power += fft2.shiftized
        fft2.draw()
    power /= len(latitude_range)
    fft.draw(target=power, folder=folder_name, latitude=lat, levels=levels,
             plot_name='{0}-{1}mean'.format(latitude_range[0], latitude_range[-1]), xcood=fft2.xcood, ycood=fft2.ycood)

    longitude_range = np.arange(west_longitude, east_longitude, 1)
    for i, lon in enumerate(longitude_range):
        fft2 = fft.FFT2(variable_slice.sel(longitude=lon), longitude=lon, folder=folder_name, levels=levels)
        if i == 0:
            power = fft2.shiftized
        else:
            power += fft2.shiftized
        fft2.draw()
    power /= len(longitude_range)
    fft.draw(target=power, folder=folder_name, longitude=lon, levels=levels,
             plot_name='{0}-{1}mean'.format(longitude_range[0], longitude_range[-1]), xcood=fft2.xcood,
             ycood=fft2.ycood)


if __name__ == '__main__':
    main()
