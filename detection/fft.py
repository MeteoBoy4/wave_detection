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


def draw(target, folder=None, plot_name=None, latitude=None, longitude=None, xcood=None, ycood=None, levels=None):
    if latitude and plot_name:
        ppath = '{0}/fft-D-latitude{1}.ps'.format(folder, plot_name)
    elif latitude:
        ppath = '{0}/fft-D-latitude{1}.png'.format(folder, latitude)
    if longitude and plot_name:
        ppath = '{0}/fft-D-longitude{1}.ps'.format(folder, plot_name)
    elif longitude:
        ppath = '{0}/fft-D-longitude{1}.png'.format(folder, longitude)

    fig, ax = plt.subplots(1)
    if np.any(levels):
        dmap = draw_2d_fft(ax, fig, xcood, ycood, target, {'levels': levels})
    else:
        dmap = draw_2d_fft(ax, fig, xcood, ycood, target)

    if not os.path.exists(folder):
        os.makedirs(folder)
    plt.savefig(ppath)
    plt.close()

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
                print('You must provide latitude when calculate zonally')

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
        draw(self.shiftized, folder=self.folder, latitude=self.latitude, longitude=self.longitude, xcood=self.xcood,
             ycood=self.ycood, levels=self.levels)


