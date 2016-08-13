#!/home/meteoboy4/anaconda/envs/Meteorology/bin/python

import xarray as xr
import numpy as np
import math
import os
import matplotlib.pyplot as plt


def draw_2d_fft(ax, fig, xcood, ycood, var, param_dict={}):
    dmap = ax.contourf(xcood, ycood, var, **param_dict)
    ax.set_ylim(0, max(ycood))
    fig.colorbar(dmap, orientation='horizontal')
    return dmap


class FFT(object):
    def __init__(self, var, latitude=None, longitude=None, levels=None, folder=None):

        if len(var.shape) != 2:
            raise ValueError("The input value is not two dimensional!")

        self.var = var
        self.latitude = latitude
        self.longitude = longitude
        self.levels = levels
        self.folder = folder
        self.height = var.shape[0]
        self.width = var.shape[1]

    @property
    def fft2(self):
        return np.fft.fft2(self.var)

    def get_sample_spacing(self, dimth):
        dim = self.var.dims[dimth]
        if dim == 'time':
            return 1.  # The sample spacing is one day
        elif dim == 'latitude':
            return 110.  # The sample spacing between latitude is 110km
        elif dim == 'longitude':
            try:
                return 110. * math.cos(math.radians(self.latitude))  # The sample spacing between longitude
            except TypeError:
                print 'You must provide latitude when calculate zonally'

    @property
    def xcood(self):
        spacing = self.get_sample_spacing(1)
        return np.fft.fftshift(np.fft.fftfreq(self.width, d=spacing))

    @property
    def ycood(self):
        spacing = self.get_sample_spacing(0)
        return np.fft.fftshift(np.fft.fftfreq(self.height, d=spacing))

    @property
    def shiftized(self):
        shift = np.fft.fftshift(self.fft2)
        power = np.abs(shift) * np.abs(shift)
        return np.log10(power)

    def draw(self):
        if self.latitude:
            ppath = './FFT_spectrum/{0}/fft-D-latitude{1}.png'.format(self.folder, self.latitude)
        if self.longitude:
            ppath = './FFT_spectrum/{0}/fft-D-longitude{1}.png'.format(self.folder, self.longitude)
        pfolder = './FFT_spectrum/{0}'.format(self.folder)

        fig, ax = plt.subplots(1)
        if np.any(self.levels):
            dmap = draw_2d_fft(ax, fig, self.xcood, self.ycood, self.shiftized, {'levels': self.levels})
        else:
            dmap = draw_2d_fft(ax, fig, self.xcood, self.ycood, self.shiftized)

        if not os.path.exists(pfolder):
            os.makedirs(pfolder)
        plt.savefig(ppath)
        plt.close()


PLOT_FOLDER = '2013_summer_w_300hpa'
FILENAME = '/run/media/MeteoBoy4/Data/MData/ERA-Interim/2013/Wind/daily/w300daily.nc'
ncfile = xr.open_dataset(FILENAME)
variable = ncfile.w
vslice = variable.sel(longitude=np.arange(60, 120), latitude=np.arange(-10, 60), time=slice('2013-06-01', '2013-08-31'))
latrange = range(20, 41, 1)
for i, lat in enumerate(latrange):
    fft2 = FFT(vslice.sel(latitude=lat), latitude=lat, folder=PLOT_FOLDER)  # , levels = np.linspace(-12, -2, num=30))
    if i == 0:
        Dave = fft2.shiftized
    else:
        Dave += fft2.shiftized
    fft2.draw()
Dave = Dave / len(latrange)
Xcood = fft2.xcood
Ycood = fft2.ycood
fig, ax = plt.subplots(1)
draw_2d_fft(ax, fig, Xcood, Ycood, Dave)
plt.savefig('./FFT_spectrum/{0}/fft-D-latitude-20-40mean.ps'.format(PLOT_FOLDER))
plt.close()
lonrange = range(80, 101, 1)
for i, lon in enumerate(lonrange):
    fft2 = FFT(vslice.sel(longitude=lon), longitude=lon, folder=PLOT_FOLDER)  # , levels = np.linspace(-12, -2, num=30))
    if i == 0:
        Dave = fft2.shiftized
    else:
        Dave += fft2.shiftized
    fft2.draw()
Dave = Dave / len(lonrange)
Xcood = fft2.xcood
Ycood = fft2.ycood
fig, ax = plt.subplots(1)
draw_2d_fft(ax, fig, Xcood, Ycood, Dave)
plt.savefig('./FFT_spectrum/{0}/fft-D-longitude-80-100mean.ps'.format(PLOT_FOLDER))
plt.close()
# plt.show()
