#!/usr/bin/env python

from __future__ import division

import os
from bisect import bisect_left
from distutils.util import strtobool

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import lmfit

class Fitter(object):
    '''
    Class that can fit functions to data
    '''

    def __init__(self, fit_func='gauss'):
        if fit_func == 'gauss' or fit_func == 'gaus':
            self.func = self._gauss
        elif fit_func == 'lorentz':
            self.func = self._lorentz
        else:
            print 'Function unknown:', fit_func

    def _gauss(self, x, params):
        '''
        Returns shifted gaussian
        '''
        x0 = params['x0'].value
        sigma = params['sigma'].value
        norm = params['norm'].value
        offset = params['offset'].value
        return offset + norm*np.exp(-.5*((x-x0)/sigma)**2)

    def _lorentz(self, x, params):
        '''
        Returns shifted lorentzian
        '''
        x0 = params['x0'].value
        hwhm = params['hwhm'].value
        norm = params['norm'].value
        offset = params['offset'].value
        return offset + norm/(1 + ((x-x0)/hwhm)**2)


    def residuals(self, params, x_data, y_data, y_sigma):
        '''
        Calculates residuals
        '''
        if y_sigma is None:
            y_sigma = np.ones(y_data.shape)
        return ((y_data - self.func(x_data, params)) / y_sigma )**2


    def fit(self, x_data, y_data, params, y_sigma=None):
        '''
        Performs fit using lmfit
        '''

        result = lmfit.minimize(self.residuals, params, args=(x_data, y_data, y_sigma))
        # result on more data points
        x = np.linspace(x_data[0], x_data[-1], 1E4)
        fit_result = self.func(x, params)
        return result, params, x, fit_result


class FitsFile(object):
    '''
    Loads and holds a LISA FITS file.
    '''

    def __init__(self, filename, dtype=None):
        # define filename
        self.filename = filename
        # set some params
        self.flat_processed = False
        self.bias_processed = False
        self.dark_processed = False
        self.cal_processed = False
        # load fits data
        self.load_file()
        if not dtype:
            # ask user for dtype (cal/bias/dark/flat/light)
            self.get_dtype()
        else:
            self.dtype = dtype

    def load_file(self):
        '''
        Loads the FITS file and its header
        '''
        print 'Loading {0}'.format(self.filename)
        with fits.open(self.filename) as f:
            self.flux = f[0].data
            self.header = f[0].header
            self.exptime = self.header['exptime']

    def get_dtype(self):
        '''
        Asks user for data type.
        '''
        dtypes = ('cal', 'light', 'dark', 'bias', 'flat')
        ans = raw_input('Please enter data type {0}\n'.format(dtypes))
        while True:
            if ans not in dtypes:
                ans = raw_input('Data type \'{0}\' not recognised, try again\n'.format(ans))
            else:
                break

    def process_flat(self, flat):
        '''
        Divide frame by flat.
        '''
        if not self.dtype in ('light', 'cal'):
            print 'Cannot divide {0} frame by flat frame, not doing anything'.format(self.dtype)
        elif self.flat_processed:
            print 'Already processed flat frame, not doing anything'
        else:
            # scale flat so max flux is 1
            flat_flux_scaled = flat.flux / np.amax(flat.flux)
            # find where pixels are ~zero
            near_zero = np.where(self.flux <= 500)
            # divide frame by flat
            self.flux /= flat_flux_scaled
            # set near zero values to zero again
            self.flux[near_zero] = 0
            self.flat_processed = True

    def process_dark(self, dark):
        '''
        Subtract dark from frame
        '''
        if not self.dtype in ('light'):
            print 'Cannot subtract dark frame from {0} frame, not doing anything'.format(self.dtype)
        elif self.dark_processed:
            print 'Already processed flat frame, not doing anything'
        else:
            # subtract dark from frame
            self.flux -= dark.flux
            self.dark_processed = True

    def process_bias(self, bias):
        '''
        Subtract bias from frame
        '''
        if self.dtype == 'bias':
            print 'Cannot subtract bias frame from {0} frame, not doing anything'.format(self.dtype)
        elif self.bias_processed:
            print 'Already processed bias, not doing anything'
        else:
            # subtract bias from frame
            self.flux -= bias.flux
            self.bias_processed = True

    def collapse(self):
        '''
        Converts to 1D by summing along y axis.
        '''
        self.flux = np.sum(self.flux, axis=0)


def shorten(index_data, val_min, val_max, *outdata):
    '''
    Select interval from data (1D or 2D)

    Arguments:
        index_data    <array>    data got get indices from (e.g. wavelengths)
        val_min    <float>    minimum value in index_data to be used
        val_max    <float>    maximum value in index_data to be used
        *outdata    <array>    one or more lists of data that are shortened (e.g. flux)
    Returns:
        outlist    <array>    list with shortened index_data and outdata
    '''

    i_min = bisect_left(index_data, val_min)
    i_max = bisect_left(index_data, val_max)
    data_short = index_data[i_min:i_max]
    outlist = [data_short]
    for d in outdata:
        if len(d.shape) == 2:
            d = d.transpose()
        d = d[i_min:i_max]
        if len(d.shape) == 2:
            d = d.transpose()
        outlist.append(d)

    return outlist


def line_flux(lamb, flux):
    '''
    Estimates total line flux
    '''
    # neon lines in angstrom
    linelist = (5031.3504, 5400.5617, 5562.7662, 5656.5664, 5689.8163, 5719.2248, 5748.2985, 5764.4188, 5804.4496, 5820.1558, 5852.4878, 5881.895, 5944.8342, 5975.534, 6029.9971, 6074.3377, 6096.1631, 6143.0626, 6163.5939, 6217.2812, 6266.495, 6304.789, 6334.4278, 6382.9917, 6402.246, 6506.5281, 6532.8822, 6598.9529, 6678.2764, 6717.043, 6929.4673, 7032.4131, 7173.9381, 7245.1666)
    # estimate total line flux
    flux_tot = 0
    for line in linelist:
        px = bisect_left(lamb, line)
        flux_tot += np.sum(flux[px-1:px+2])
    return flux_tot


def wav_cal(cal, show_plots=False):
    '''
    Determines wavelength calibration from 1D calibration frame.
    '''
    # approximate range:
    # ~ 3800-7300 AA in 1391 pixels
    startx = 3770
    # neon lines in angstrom
    linelist = (5031.3504, 5400.5617, 5562.7662, 5656.5664, 5689.8163, 5719.2248, 5748.2985, 5764.4188, 5804.4496, 5820.1558, 5852.4878, 5881.895, 5944.8342, 5975.534, 6029.9971, 6074.3377, 6096.1631, 6143.0626, 6163.5939, 6217.2812, 6266.495, 6304.789, 6334.4278, 6382.9917, 6402.246, 6506.5281, 6532.8822, 6598.9529, 6678.2764, 6717.043, 6929.4673, 7032.4131, 7173.9381, 7245.1666)
    # estimate wavelength scale
    lamb = np.linspace(startx, 7350, len(cal.flux))
    deltax = lamb[1] - lamb[0]
    
    if show_plots:
        # make a plot of theoretical line positions
        fig, ax = plt.subplots()
        ax.plot(lamb, cal.flux)
        for line in linelist:
            ax.axvline(line, ls='--', c='k')
        fig.suptitle('Theoretical line positions')
        plt.show()

    # find line flux for different shifts in the wavelength scale
    shifts = np.linspace(-20, 20, 100)
    fluxes = np.asarray([line_flux(lamb+deltax*shift, cal.flux) for shift in shifts])
    # shorten data for fitting
    shift_min = -5
    shift_max = 5
    shifts_short, fluxes_short = shorten(shifts, shift_min, shift_max, fluxes)
    
    # fit maximum with a gaussian
    params = lmfit.Parameters()
    fluxmax_guess = np.amax(fluxes_short)
    shift_guess = shifts_short[np.argmax(fluxes_short)]
    params.add('x0', shift_guess, min=shift_guess-2, max=shift_guess+2, vary=True)
    params.add('sigma', 5, min=2, max=10, vary=True)
    params.add('norm', fluxmax_guess, min=.98*fluxmax_guess, max=1.02*fluxmax_guess, vary=True)
    params.add('offset', 0, vary=False)
    fitter = Fitter(fit_func='gauss')
    result, params, x, fit_result = fitter.fit(shifts_short, fluxes_short, params)

    # save best fit final shift
    shift = params['x0'].value

    if show_plots:
        # make a plot of best fit shift
        fig, ax = plt.subplots()
        # data
        ax.scatter(shifts, fluxes, c='k')
        # best fit
        ax.plot(x, fit_result, c='r')
        ax.set_xlim(shifts[0], shifts[-1])
        ax.set_xlabel('Pixel shift')
        ax.set_ylabel('Total line flux')
        fig.suptitle('Correlation')
        plt.show()

    # find best pixel position for all lines
    linelist_px = [(line - startx)/deltax for line in linelist] + shift
    positions = []
    for line in linelist_px:
        maxval = np.amax(cal.flux[line-2:line+3])
        pos = np.where(cal.flux == maxval)[0]
        positions.append(pos[0])
    # find best fit wavelength scale
    deltax, startx = np.polyfit(positions, linelist, deg=1)
    print startx, deltax
    lamb = startx + deltax * np.arange(len(cal.flux))
    # create O-C data
    o_c = np.asarray(linelist) - (startx + deltax * linelist_px)

    if show_plots:
        # plot best fit wavelength scale
        fig, ax = plt.subplots()
        ax.scatter(positions, linelist)
        # add best fit wavelength scale
        ax.plot(lamb, c='r')
        ax.set_xlim(positions[0], positions[-1])
        ax.set_xlabel('Pixel')
        ax.set_ylabel(r'Wavelength ($\AA$)')
        fig.suptitle('Pixel to wavelength conversion')

        # Create O-C plot
        fig, ax = plt.subplots()
        ax.scatter(linelist_px, o_c)
        ax.set_xlim(positions[0], positions[-1])
        ax.set_xlabel('Pixel')
        ax.set_ylabel('O-C')
        fig.suptitle('O-C wavelenghts')
        plt.show()

    # return wavelength scale and flux
    return lamb


def plot(frame):
    '''
    Plots 1D or 2D image
    '''
    flux = frame.flux
    ndims = len(flux.shape)
    fig, ax = plt.subplots()
    
    if ndims == 1:
        ax.plot(flux)
    elif ndims == 2:
        ax.imshow(flux, vmin=0, vmax=np.amax(flux), cmap='gray')
    else:
        print 'Number of dimensions is wrong: {0:.0f}'.format(ndims)
        return
    fig.suptitle(frame.dtype.title())
    plt.show()


def main():
    # linelist
    linelist = [ (r'H$\alpha$', 6562.79), (r'H$\beta$', 4861.35), (r'H$\gamma$', 4340.472), ('Ca II H', 3968.47), ('Ca II K', 3933.66) ]
    # define files
    data_dir = 'test_data'
    light_file = 'jupiter.fit'
    flat_file = 'flat.fit'
    cal_file = 'neon.fit'
    # load files
    light = FitsFile(os.path.join(data_dir, light_file), dtype='light')
    flat = FitsFile(os.path.join(data_dir, flat_file), dtype='flat')
    cal = FitsFile(os.path.join(data_dir, cal_file), dtype='cal')

    # apply flat correction
    light.process_flat(flat)
    cal.process_flat(flat)

    # convert light and calibration to 1D
    light.collapse()
    cal.collapse()

    # wavelength calibration
    lamb = wav_cal(cal, show_plots=True)
    # save flux of light frame
    flux = light.flux

    # plot science spectrum
    fig, ax = plt.subplots()
    ax.plot(lamb, flux)
    # add some lines
    for name, wavelength in linelist:
        wav_index = bisect_left(lamb, wavelength)
        fl_wav = np.average(flux[wav_index-5:wav_index+5])
        ax.axvline(wavelength, ls='--', color='k')
        ax.annotate(name, xy=(wavelength, 1.5*fl_wav), horizontalalignment='center', bbox=dict(facecolor='w', edgecolor='None'))
    ax.set_xlim(lamb[0], lamb[-1])
    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel('Counts')
    fig.suptitle('Science spectrum')
    plt.show()

    # save science spectrum
    try:
        ans = strtobool(raw_input('Save science spectrum? [Y/n]\n'))
    except ValueError:
        ans = 1
    if ans == 1:
        save_name = raw_input('Please input filename. (output.dat)')
        if save_name == '':
            save_name = 'output.dat'
        # save the file
        data = np.array([lamb, flux])
        np.savetxt(save_name, data.T, fmt='%.2f %.0f')
        print 'Saved data to {0}'.format(save_name)


if __name__ == '__main__':
    main()
