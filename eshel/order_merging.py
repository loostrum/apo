#!/usr/bin/env python
'''
Order merging for normalized eShel spectra.
Author: Leon Oostrum
E-Mail: l.c.oostrum@uva.nl
'''

from __future__ import division
import os
import sys
import glob
from distutils.util import strtobool
from bisect import bisect_left, bisect_right

import numpy as np
import matplotlib.pyplot as plt


def read_file(filename):
    '''
    Read an ascii file as produced by specnorm.py.
    '''
    wave, flux = np.loadtxt(filename).T
    return list(wave), list(flux)

if __name__ == '__main__':
    if not len(sys.argv) in (2, 3):
        print 'Usage: order_merging.py DATADIR SAVEDIR (optional)'
        sys.exit(1)
    # set data dir
    data_dir = sys.argv[1]
    if not os.path.isdir(data_dir):
        print 'Data directory does not exist.'
        sys.exit(1)

    # set save dir
    if len(sys.argv) == 2:
        # save dir was not specified
        try:
            ans = strtobool(raw_input('Save directory not specified. Use data directory? [Y/n]\n'))
        except ValueError:
            ans = 1
        if ans:
            save_dir = data_dir[:]
        else:
            sys.exit(1)
    else:
        save_dir = sys.argv[2]
        # check if save dir exists
        if not os.path.isdir(save_dir):
            print 'Save directory does not exist.'
            sys.exit(1)

    # get list of files
    filelist = glob.glob(os.path.join(data_dir, '*P_1B_[0-9][0-9]_norm.dat'))
    pre = filelist[0].split('P_1B')[0] + 'P_1B_'
    aft = '_norm.dat'
    # Get object name
    obj = filelist[0].split('-')[2]
    # create dict which returns filename as function of order
    files = dict([(int(f.split(pre)[1].split(aft)[0]), f) for f in filelist])
    # get list of orders. Reverse to get shorter wavelengths first
    orders = sorted(files.keys(), reverse=True)

    # load first order
    print 'Processing order {0} (1/{1}'.format(orders[0], len(orders))
    wave_full, flux_full = read_file(files[orders[0]])
    # array to save order merge locations
    merge_locs = []
    # loop over orders, but skip first one
    for i, order in enumerate(orders[1:]):
        print 'Processing order {0} ({1}/{2})'.format(order, i+1, len(orders))
        # load data
        wave, flux = read_file(files[order])
        # find overlap with previous order
        min_old = bisect_left(wave_full, wave[0])
        max_new = bisect_right(wave, wave_full[-1])
        # save merge locations
        merge_locs.append([wave[0], wave_full[-1]])
        # average the overlapping part
        part1 = np.array(flux_full[min_old:])
        part2 = np.array(flux[:max_new])
        flux_avg = list(np.mean(np.array([part1, part2]), axis=0))
        # add to final data
        wave_full.extend(wave[max_new:])
        flux_full[min_old:] = flux_avg
        flux_full.extend(flux[max_new:])

    # make a plot if needed
    try:
        ans = strtobool(raw_input('Show the spectrum [Y/n]?\n'))
    except ValueError:
        ans = 1
    if ans == 1:
        fig, ax = plt.subplots()
        ax.plot(wave_full, flux_full, c='k')
        for wav_min, wav_max in merge_locs:
            ax.axvspan(wav_min, wav_max, alpha=.1)
        ax.set_xlim(wave_full[0], wave_full[-1])
        ax.set_ylim(ymin=max(ax.set_ylim()[0], 0))
        ax.ticklabel_format(axis='both', useOffset=False)
        ax.set_xlabel(r'Wavelength ($\AA$)')
        ax.set_ylabel('Flux (norm.)')
        fig.suptitle(obj)
        plt.show()

    # save spectrum
    try:
        ans = strtobool(raw_input('Save the spectrum? [Y/n]?\n'))
    except ValueError:
        ans = 1
    if ans == 1:
        filename = os.path.join(save_dir, obj.lower()+'_merged.dat')
        if os.path.isfile(filename):
            try:
                ans = strtobool(raw_input('File already exists: {0}, overwrite? [Y/n]?\n'.format(filename)))
            except ValueError:
                ans == 1
            if not ans:
                exit()
        data = np.array([wave_full, flux_full]).T
        np.savetxt(os.path.join(save_dir, obj.lower()+'_merged.dat'), data)
