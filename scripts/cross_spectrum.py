#!/home/meteoboy4/anaconda/envs/Meteorology/bin/python

import sys
import argparse
import numpy as np
import datetime
import xarray as xr

from detection import fft


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Cross spectrum analysis using fft along both latitude and longitude")
    parser.add_argument("diversion", help="the path to the netCDF file of 'diversion'")
    parser.add_argument("vorticity", help="the path to the netCDF file of 'vorticity'")
    parser.add_argument("time", type=int, default=41, help="the time to analyze (in indexing order)")
    parser.add_argument("output", help="the output directory to which the plots go")

    args = parser.parse_args()

    time = args.time
    div_file = xr.open_dataset(args.diversion)
    vort_file = xr.open_dataset(args.vorticity)
    div = div_file[div_file.keys()[-1]]
    vort = vort_file[vort_file.keys()[-1]]
    div_slice = div.sel(longitude=slice(120, 160), latitude=slice(60, 30)).isel(time=time)
    vort_slice = vort.sel(longitude=slice(120, 160), latitude=slice(60, 30)).isel(time=time)
    timestring = datetime.datetime.utcfromtimestamp(div_slice.time.data.tolist()/1e9).strftime('%Y-%m-%d_%H:%M')
    folder_name = args.output + '_{0}'.format(timestring)

    latitude_range = range(30, 60, 1)
    for lat in latitude_range:
        cross = fft.CrossSpectrum(div_slice.sel(latitude=lat), vort_slice.sel(latitude=lat), latitude=lat, folder=folder_name)
        cross.draw()

if __name__ == '__main__':
    main()
