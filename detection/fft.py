#!/home/meteoboy4/anaconda/envs/Meteorology/bin/python

import xarray as xr
import numpy as np
import math
import os
import matplotlib.pyplot as plt


def draw_2d_fft(ax, fig, xcood, ycood, var, param_dict={}):
    dmap = ax.contourf(xcood, ycood, var, **param_dict)
    ax.set_ylim(0, max(ycood))
    ax.set_xlim(-max(xcood)/6., max(xcood)/6.)
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


def draw_cross(target, folder=None, latitude=None, longitude=None, xcood=None):
    if not os.path.exists(folder):
        os.makedirs(folder)
    path_amplitude = '{0}/amplitude-latitude{1}.png'.format(folder, latitude)
    path_angle = '{0}/angle-latitude{1}.png'.format(folder, latitude)
    path_condense = '{0}/condense-latitude{1}.png'.format(folder, latitude)
    path = {'amplitude':path_amplitude, 'angle':path_angle, 'condense':path_condense}

    for keys in target.keys():
        fig, ax = plt.subplots()
        ax.plot(xcood, target[keys])
        plt.savefig(path[keys])
        plt.close()

    # plt.plot(xcood, target['amplitude'])
    # plt.axis([-max(xcood) / 6., max(xcood) / 6., 0, 1.2 * max(target['amplitude'])])
    # plt.savefig(path['amplitude'])
    # plt.close()

    # plt.plot(xcood, target['angle'])
    # plt.axis(-max(xcood) / 6., max(xcood) / 6., 0, max(target['angle']))
    # plt.savefig(path['angle'])
    # plt.close()
    #
    # plt.plot(xcood, target['condense'])
    # plt.axis(-max(xcood) / 6., max(xcood) / 6., 0, max(target['condense']))
    # plt.savefig(path['condense'])
    # plt.close()

class FFT2(object):
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
            return 110. * abs(self.var.latitude.data[0] - self.var.latitude.data[1])  # The sample spacing between one latitude is 110km
        elif dim == 'longitude':
            try:
                return 110. * math.cos(math.radians(self.latitude)) * abs(self.var.longitude.data[0] - self.var.longitude.data[1])  # The sample spacing between longitude
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


class CrossSpectrum(object):
    def __init__(self, var1, var2, latitude, folder=None):
        self.var1 = var1
        self.var2 = var2
        self.width = var1.shape[0]
        self.latitude = latitude
        self.folder = folder

    @property
    def fft_first(self):
        return np.fft.fft(self.var1)

    @property
    def fft_second(self):
        return np.fft.fft(self.var2)

    def get_sample_spacing(self, dimth):
        dim = self.var1.dims[dimth]
        if dim == 'time':
            return 1.  # The sample spacing is one day
        elif dim == 'latitude':
            return 110. * abs(self.var1.latitude.data[0] - self.var1.latitude.data[1])  # The sample spacing between one latitude is 110km
        elif dim == 'longitude':
            try:
                return 110. * math.cos(math.radians(self.latitude)) * abs(self.var1.longitude.data[0] - self.var1.longitude.data[1])  # The sample spacing between longitude
            except TypeError:
                print('You must provide latitude when calculate zonally')

    @property
    def xcood(self):
        spacing = self.get_sample_spacing(0)
        return np.fft.fftfreq(self.width, d=spacing)

    @property
    def part1(self):
        x = self.fft_first
        y = self.fft_second
        return x.real * y.real + x.imag * y.imag

    @property
    def part2(self):
        x = self.fft_first
        y = self.fft_second
        return x.imag * y.real - x.real * y.imag

    @property
    def amplitude(self):
        return np.sqrt(self.part1 ** 2 + self.part2 ** 2)

    @property
    def angle(self):
        return np.degrees(np.arctan(self.part2 / self.part1))

    @property
    def condense(self):
        return self.amplitude ** 2 / np.abs(self.fft_first) ** 2 / np.abs(self.fft_second) ** 2

    def draw(self):
        target = {'amplitude': self.amplitude, 'angle': self.angle, 'condense': self.condense}
        draw_cross(target, folder=self.folder, latitude=self.latitude, xcood=self.xcood)
