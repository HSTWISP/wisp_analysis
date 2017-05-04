#!/usr/bin/env python
##########################################################################
##########################################################################
# get_z_interactive.py By Nathaniel R. Ross, UCLA, nross@astro.ucla.edu
# usage: python get_z_interactive.py line_list
# Reads in list of emission lines from the WISP survey HST/WFC3/IR Grisms and
# plots the spectra, iterates through the detected emission lines, allowing the
# user to make line identifications or reject spurious lines quickly and
# efficiently.
#
# Version 1.0 updates the program to look for wavelength-calibrated 2d grism
# stamps first. Also, will now look for default line list name and make
# linelistfile an optional parameter. Furthermore, I have added the option
# to save the incomplete line list to file as you progress, but default is NOT to
# do this.
# Optional arguments for go() are 1. linelistfile (String) - path to line list file
#                                 2. save_temp (Bool) - Save progress in
#                                    linelist/Par???lines_with_redshifts.incomplete.dat
#                                 3. recover_temp (Bool) - Recover progress from previous
#                                    session using the .incomplete.dat files
#
# ** Major change in version 1.0:
# We now use xpa instead of IRAF to display the 2d grism stamps and full 2d
# direct images (instead of cutout stamp). Reason for this change was the desire
# to show the grism stamps that have a wavelength solution applied. When using
# the IRAF display command, ds9 would fail to recognize this coordinate system.
# The XPA system can be downloaded here: http://hea-www.harvard.edu/RD/xpa/index.html
#
##########################################################################
import os
from glob import glob
import distutils
#import numpy as np
import fileinput
import scipy
import pylab as plt
from scipy.interpolate import spline
from astropy.table import Table
from distutils.sysconfig import *  # question-- what is this for?
import sys
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
# Explicitly import readline to make the text entry process less tortuous
# on OSX
import readline
# SQLLite database support for data persistence
from WISPLFDatabaseManager import WISPLFDatabaseManager as WDBM

from wisp_analysis import *


#######################################################################
# define wavelengths of lines of interest
# this is super not the way to do this, but oh well
lam_Halpha = 6563.0
lam_Hbeta = 4861.0
lam_Hg = 4341.0
lam_Oiii_1 = 4959.0
lam_Oiii_2 = 5007.0
lam_Oii = 3727.0
lam_Sii = 6724.0
lam_Siii_1 = 9069.0
lam_Siii_2 = 9532.0
# lam_Lya=1216.0
lam_He = 10830.0
# lam_Fe=12600.0
# lam_Pag=10940.0
# lam_Pab=12810.0

suplines = [lam_Oii, lam_Hg, lam_Hbeta, lam_Oiii_2,
            lam_Halpha, lam_Sii, lam_Siii_1, lam_Siii_2, lam_He]
suplines_str = ['[OII]', r'H$\gamma$', r'H$\beta$', '[OIII]',
                r'H$\alpha$', '[SII]', '[SIII]', '[SIII]', 'HeI']
#######################################################################


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def getzeroorders(zeroorderpath, g='G141', magcut=23.5):  # MB: changed from 23.0
    """
    Changing to return a table
    """
    zop = open(zeroorderpath, 'r')
    zox = []
    zoy = []
    zoid = []
    zmag = []
    for line in zop:
        if len(line) > 60:
            linesplit = line.split()
            zox.append(float(linesplit[1][0:-1]))
            zoy.append(float(linesplit[2][0:-3]))
            zoid.append(int(linesplit[-2].split('{')[-1]))
            zmag.append(float(linesplit[-1][1:-2]))  # get mag from reg file

    zop.close()
    zoid = np.array(zoid)
    zoy = np.array(zoy)
    zox = np.array(zox)
    zmag = np.array(zmag)
    cond = (zmag <= magcut)
    t = Table([zox[cond], zoy[cond], zoid[cond]], names=('x', 'y', 'objid'))
    return t


def getfirstorders(firstorderpath):
    """
    Changing to return a table
    """
    fop = open(firstorderpath, 'r')
    fox = []
    foy = []
    folen = []
    fowid = []
#    foid=[] nothing is done with 1st order IDs anymore
    for line in fop:
        # if line[0]!='#':
        linesplit = line.split()
        fox.append(float(linesplit[1][0:-1]))  # [0:-1] strips off the comma.
        foy.append(float(linesplit[2][0:-1]))
        folen.append(float(linesplit[3][0:-1]))
        # python is weird.
        fowid.append(float(linesplit[-1].split('{')[-1].split('}')[0]))
    print foid
    t = Table([fox, foy, folen, fowid], names=(
        'x', 'y', 'len', 'width', 'objid'))
    return t


def get_remaining_objects(full_obj_list, objid_done):
    # array of flags
    wdone = np.in1d(np.array(full_obj_list), objid_done)
    remaining = np.copy(full_obj_list)
    mask = np.ones(full_obj_list.shape, dtype=bool)
    mask[wdone] = False
    remaining = remaining[mask]
    return remaining


def print_help_message():
    """
    Just putting this here to keep it out of the way.
    """
    msg = "Enter option (read carefully, options have changed): \n \
             \t a = accept object fit \n \
             \t ac = accept object fit, noting contamination\n  \
             \t r = reject object \n \
             \t z = enter a different z guess  \n \
             \t w = enter a different emission line wavelength guess  \n \
             \t ha,  or hb, o31, o32, o2, s2, s31, s32 = change redshift guess \n \
             \t n = skip to next brightest line found in this object \n \
             \t fw = change the fwhm guess in pixels \n \
             \t c = add comment \n \
             \t contam = specify contamination to line flux and/or continuum \n \
             \t f = set/unset flag(s) \n \
             \t t = change transition wavelength \n \
             \t m1, m2, or m3 =mask up to three discontinuous wavelength regions \n \
             \t nodes = change the wavelengths for the continuum spline nodes \n \
             \t bluecut = change the blue cutoff of the G102 spec \n \
             \t redcut  = change the red cutoff of the G141 spec \n \
             \t reset = reset interactive options back default for this object \n \
             \t lin = linear z-scale \n \
             \t log = logarithmic  \n \
             \t zs102 = z1,z2 comma-separated range for G102 zscale \n  \
             \t zs141 = z1,z2 comma-separated range for G141 zscale \n  \
             \t dc = recenter direct images \n \
             \t reload = reload direct images \n \
             \t dr = reload direct image reg files\n \
             \t h = print this message\n \
             \t q = quit\n"
    print msg


def write_object_summary(par, obj, fitresults, snr_meas_array, contamflags):
    """ """
    # string names for output 
    linenames = np.array(['[OII]', 'Hgamma', 'Hbeta', '[OIII]', \
                          'Halpha', '[SII]', '[SIII]', '[SIII]', 'HeI'])
    # string names for accessing fitresults
    fluxstrs = ['oii','hg','hb','oiii','hanii','sii','siii_9069','siii_9532']
    linefluxes = np.array([fitresults['%s_flux'%fs] for fs in fluxstrs])

    # initial message
    print '#'*72
    msg = '## Par{} Obj {}:\n##   Fit Redshift: z = {:.4f}\n'.format(par, obj, fitresults['redshift'])

    # lines with S/N > 3
    good_snr = np.where(snr_meas_array > 3)
    msg = msg + '##   Lines fit with S/N > 3:\n'
    for gsnr in good_snr[0]:
        msg = msg + '##\t%s: Flux = %.3e    S/N = %.2f\n'%(linenames[gsnr], 
                                linefluxes[gsnr], snr_meas_array[gsnr])

    cfout = ['%s:%i'%(cf,contamflags[cf]) for cf in contamflags if contamflags[cf] > 0]
    msg = msg + '##   Contamination flags set:\n##\t' + ', '.join(cfout)
    print(msg)
    print '#'*72


def plot_object(zguess, spdata, config_pars, snr_meas_array, full_fitmodel, full_contmodel, lamlines_found, index_of_strongest_line, contmodel, plottitle):
    """
    # save the figure for everything, junk objects and all
    # previous figures are overwritten
    """
    # the expected wavelengths of emission lines given the zguess
    lamobs = (1 + zguess) * np.array(suplines)

    plotfilename = os.path.join('figs', '%s_fit.png' % plottitle)

    spec_lam = spdata[0]
    spec_val = spdata[1]
    spec_unc = spdata[2]
    spec_con = spdata[3]
    spec_zer = spdata[4]
    # apply the mask to the wavelength array
    masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)

    plt.ion()
    fig = plt.figure(1, figsize=(11, 8))
    plt.clf()
    gs = gridspec.GridSpec(3, 4)
    ax1 = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2:, :])

    xmin = np.ma.min(spec_lam) - 200.0
    xmax = np.ma.max(spec_lam) + 200.0
    ymin = np.ma.min(spec_val)
    ymax = 1.5 * np.ma.max(spec_val)

    ax1.plot(spec_lam, spec_val, 'k', spec_lam, spec_con, 'r', ls='steps')
    ax1.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)

    # transforms for plotting in data and axes coordinates
    ax1trans = mtransforms.blended_transform_factory(
        ax1.transData, ax1.transAxes)
    ax2trans = mtransforms.blended_transform_factory(
        ax2.transData, ax2.transAxes)
    # plot observed wavelengths of all the possible lines.
    for li, lstring, sn_meas in zip(lamobs, suplines_str, snr_meas_array):
        if (li > xmin + 100) & (li < xmax - 100):
            for ax in [ax1, ax2]:
                ax.axvline(x=li, color='b')
            stringplot = lstring + '   (' + str(round(sn_meas, 2)) + ')'
            # use data coordinates for x-axis and axes coords for y-axis
            ax1.text(li, 0.85, stringplot, rotation='vertical',
                     ha='right', fontsize='16', transform=ax1trans)

    ax1.plot(spec_lam, full_fitmodel, color='r', lw=1.5)
    ax1.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=1.5)

    # plot 0th orders
    w = np.where(spec_zer == 3)
    spec_zero_bad = spec_zer * 0 - 1
    spec_zero_bad[w] = 1.
    # mild zeroth orders
    w = np.where(spec_zer == 2)
    spec_zero_mild = spec_zer * 0 - 1
    spec_zero_mild[w] = 1.
    for ax in [ax1, ax2]:
        # use data coordinates for x-axis and axes coords for y-axis
        trans = mtransforms.blended_transform_factory(
            ax.transData, ax.transAxes)
        ax.fill_between(spec_lam, 0, 1, where=spec_zero_bad == 1, color='red',
                        alpha=0.3, transform=trans, label='Major 0th order contam')
        ax.fill_between(spec_lam, 0, 1, where=spec_zero_mild == 1, color='orange',
                        alpha=0.3, transform=trans, label='Minor 0th order contam')

    # plot any masked regions
    for mr, label in zip(['mask_region1', 'mask_region2', 'mask_region3'], ['masked regions', None, None]):
        if (config_pars[mr][0] != 0.) & (config_pars[mr][1] != 0.):
            for ax in [ax1, ax2]:
                trans = mtransforms.blended_transform_factory(
                    ax.transData, ax.transAxes)
                ax.fill_between(
                    config_pars[mr], 0, 1, color='grey', alpha=0.3, transform=trans, label=label)
    ax1.legend(bbox_to_anchor=[1.05, 1.15])

    # find values of spec_lam nearest to the nodes
    nodelam = config_pars['node_wave']
    nl_arr = []
    cont_node = []
    for nl in nodelam:
        w = np.argmin(np.abs(spec_lam - nl))
        nl_arr.append(spec_lam[w])
        cont_node.append(full_contmodel[w])
    ax1.plot(nl_arr, cont_node, 'ko', ms=9)

    # repeat for line_candidates
    lf_lam = []
    lf_cont = []
    for lf in lamlines_found:
        w = np.argmin(np.abs(spec_lam - lf))
        lf_lam.append(spec_lam[w])
        lf_cont.append(full_contmodel[w])
    ax1.plot(lf_lam, lf_cont, 'bo', ms=9)

    # indicate "current" line
    current_lam = lamlines_found[index_of_strongest_line]
    current_cont = contmodel[
        np.argmin(np.abs(np.ma.compressed(masked_spec_lam) - current_lam))]
    ax1.plot(current_lam, current_cont, 'ro', ms=10)

    ax1.set_ylabel(
        r'F$_\lambda$ ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$', size='xx-large')
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    ax1.set_title(plottitle)

    # second panel for s/n
    s2n = (spec_val - full_contmodel) / spec_unc
    s2n_lam = spec_lam
    mask = np.logical_and(s2n > -10000., s2n < 10000.)
    s2n = s2n[mask]
    s2n_lam = s2n_lam[mask]
    ax2.plot(s2n_lam, s2n, 'k-', linestyle='steps')
    ymin = s2n.min()
    ymax = 1.5 * s2n.max()
    ax2.axhline(y=config_pars['n_sigma_above_cont'], c='r')
    for li in lamobs:
        ax2.axvline(x=li, color='b')
    ax2.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)
    ax2.set_xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
    ax2.set_ylabel(r'S/N', size='xx-large')
    ax2.set_xlim([xmin, xmax])
    ax2.set_ylim(ymin, ymax)
    # fig = plt.gcf() a
    fig.savefig(plotfilename)
    plt.draw()


def inspect_object(par, obj, objinfo, lamlines_found, ston_found, g102zeros, g141zeros, linelistoutfile, commentsfile, remaining, allobjects, show_dispersed=True):
    """An attempt to move all object-specific tasks
    """
    # set up and filenames
    specnameg102 = 'Par%i_G102_BEAM_%iA.dat' % (par, obj)
    specnameg141 = 'Par%i_G141_BEAM_%iA.dat' % (par, obj)
    plottitle = 'Par%i_BEAM_%i' % (par, obj)
    fitdatafilename = 'fitdata/%s_fitspec' % plottitle
    # read in 1D spectrum
    if os.path.exists(specnameg102):
        tab_blue = asciitable.read(
            specnameg102, names=['lambda', 'flux', 'ferror', 'contam', 'zero'])
    else:
        tab_blue = None
    # check for g141, too. there are maybe 4 G102-only fields
    if os.path.exists(specnameg141):
        tab_red = asciitable.read(
            specnameg141, names=['lambda', 'flux', 'ferror', 'contam', 'zero'])
    else:
        tab_red = None

    # display the object
    if g102zeros is not None:
        #show2dNEW('G102', par, obj, g102firsts, g102zeros, 'linear')
        show2dNEW('G102', par, obj, g102zeros, 'linear')
    if g141zeros is not None:
        show2dNEW('G141', par, obj, g141zeros, 'linear')
    # pan full images to the new object
    showDirectNEW(obj)
    if show_dispersed:
        showDispersed(obj)

    # define parameters for this object
    ra = objinfo['ra']
    dec = objinfo['dec']
    a_image = objinfo['a_image']
    b_image = objinfo['b_image']
    jmag = objinfo['jmag']
    jerr = objinfo['jerr']
    hmag = objinfo['hmag']
    herr = objinfo['herr']

    # start with a fresh set of config pars
    config_pars = read_config('default.config')

    # Data have been successfully loaded for this object. If it has been inspected
    # previously, the original results will have been stored in the SQLLite database
    # and a retrieval obtion should be offered.
    databaseManager = WDBM(dbFileNamePrefix='Par{}'.format(par))
    mostRecentObjectData = databaseManager.getMostRecentObject()
#    if mostRecentObjectData is not None :
#        print('Most recent object in database: Par {}, Obj {}, Date {}'.format(*mostRecentObjectData))
#    else :
#        print('Database is empty.')
    catalogueEntryData = databaseManager.loadCatalogueEntry(parNumber=par, objectId=obj)
    acceptPrevFit = False
    if catalogueEntryData is not None :
        nonFitResults, fitResults = catalogueEntryData
        (par_db, obj_db, ra_db, dec_db, jmag_db, hmag_db, a_image_db, b_image_db, contamflag_db, entrytime_db) = nonFitResults
#        print('Found previous fit results for Pointing {}, Object {}.\nEnter "y" to accept the earlier fit.'.format(par_db, obj_db))
        print('You have already fit Obj {}. Refit? [y/N]').format(obj)
        acceptPrevFit = raw_input('> ').strip().lower() == 'n'
#        print('Accepting previous fit.' if acceptPrevFit else 'Re-fitting this object.')
#    else :
#        print('No previous fit results found. Fitting this object now.')

    # get line, fwhm, z estimate
    # choose the lamline that has the highest S/N estimate
    s = np.argsort(ston_found)
    # reverse s/n order

    ston_found = ston_found[s[::-1]]
    lamlines_found = lamlines_found[s[::-1]]
    index_of_strongest_line = 0
    lamline = lamlines_found[index_of_strongest_line]
    zguess = lamline / lam_Halpha - 1
    # fwhm is defined for the red side, regardless of where line is
    fwhm_guess = 2.35 * a_image * config_pars['dispersion_red']

    # print object info to screen
    print
    print "=" * 72
    print "Par%i Obj %i:" % (par, obj)
    print "Initial redshift guess: z = %f" % (zguess)
    print "\nWhat would you like to do with this object?\nSee the README for options, or type 'h' to print them all to the screen."

    comment = ' '
    contamflags = {'o2':0, 'hg':0, 'hb':0, 'o3':0, 'ha':0, 's2':0, 's31':0, \
                   's32':0, 'he1':0}
    # Skip if previous fit is to be accepted
    done = 0 if not acceptPrevFit else 1
    while (done == 0):
        # sticking with the while loop to determine whether user is finished
        # with object

        # get spectrum for obj. do this every time because sometimes we
        # re-read with a mask or a different transition wavelength
        spdata = trim_spec(tab_blue, tab_red, config_pars,
                           mask_zeros=True, return_masks=True)
        spec_lam = spdata[0]
        spec_val = spdata[1]
        spec_unc = spdata[2]
        spec_con = spdata[3]
        spec_zer = spdata[4]
        mask_flg = spdata[5]

        # apply the mask to the wavelength array
        masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)
        # compress the masked arrays for fitting
        fit_inputs = [np.ma.compressed(masked_spec_lam), np.ma.compressed(
            spec_val), np.ma.compressed(spec_unc), config_pars, zguess, fwhm_guess, str(obj)]
        # parsing the input to facilitate parallel processing when fitting
        # is done in batch mode.
        fitresults = fit_obj(fit_inputs)
        zfit = fitresults['redshift']
        fitpars = fitresults['fit_parameters']
        fitpars_nolines = cp.deepcopy(fitpars)
        fitpars_nolines[9:19] = 0.
        fitpars_nolines[11] = 1.4  # can't kill this one or divide by zero.
        fitpars_nolines[12] = 0.1
        fitmodel = emissionline_model(fitpars, np.ma.compressed(
            masked_spec_lam)) * fitresults['fit_scale_factor']
        contmodel = emissionline_model(fitpars_nolines, np.ma.compressed(
            masked_spec_lam)) * fitresults['fit_scale_factor']
        # the fitting is done on compressed arrays, so we need to
        # create masked versions of the fit and continuum models
        full_fitmodel = np.zeros(spec_lam.shape, dtype=float)
        full_contmodel = np.zeros(spec_lam.shape, dtype=float)
        full_fitmodel[np.ma.nonzero(spec_val)] = fitmodel
        full_contmodel[np.ma.nonzero(spec_val)] = contmodel
        full_fitmodel = np.ma.masked_where(
            np.ma.getmask(spec_val), full_fitmodel)
        full_contmodel = np.ma.masked_where(
            np.ma.getmask(spec_val), full_contmodel)
        # measured S/N
        snr_meas_array = np.array([fitresults['oii_flux'] / fitresults['oii_error'],
                                   fitresults['hg_flux'] /
                                   fitresults['hg_error'],
                                   fitresults['hb_flux'] /
                                   fitresults['hb_error'],
                                   fitresults['oiii_flux'] /
                                   fitresults['oiii_error'],
                                   fitresults['hanii_flux'] /
                                   fitresults['hanii_error'],
                                   fitresults['sii_flux'] /
                                   fitresults['sii_error'],
                                   fitresults['siii_9069_flux'] /
                                   fitresults['siii_9069_error'],
                                   fitresults['siii_9532_flux'] /
                                   fitresults['siii_9532_error'],
                                   fitresults['he1_flux'] / fitresults['he1_error']])

        # plot the whole goddamn thing
        plot_object(zguess, spdata, config_pars, snr_meas_array, full_fitmodel,
                    full_contmodel, lamlines_found, index_of_strongest_line, contmodel, plottitle)
#        print "    Guess Redshift: z = %f" % (zguess)
        print "    Fit Redshift:   z = %f\n" % (zfit)
        # input
        option = raw_input("> ")

        # checking user's input. keeping this format the same as before
        # any time done is set to 1, the object is considered fit

        # reject object
        if option == 'r':
            done = 1
            zset = 0

        # accept object
        elif option == 'a':
            done = 1
            zset = 1
            flagcont = 1

        # accept object and note contamination
        elif option == 'ac':
            done = 1
            zset = 1
            flagcont = 2
            # add to contamination flags
            for k,v in contamflags.iteritems():
                contamflags[k] = contamflags[k] | 1

        # change redshift guess
        elif option == 'z':
            print "The current redshift guess is: %f\nEnter Redshift Guess:" % zguess
            zguess = float(raw_input("> "))

        # change wavelength guess
        elif option == 'w':
            print "The current emission line wavelength is: %f\nEnter Wavelength Guess in Angstroms:" % lamline
            zguess = float(raw_input("> ")) / lam_Halpha - 1.

        # change the fwhm guess
        elif option == 'fw':
            print "Enter a Guess for FWHM in pixels"
            print "The current fwhm_fit is:  " + str(fitresults['fwhm_g141'] / config_pars['dispersion_red']) + " and 2*A_image is: " + str(2 * a_image)
            fwhm_guess = config_pars['dispersion_red'] * float(raw_input(">"))

        # mask out 1, 2, or 3 regions of the spectrum
        elif option == 'm1':
            print "Enter wavelength window to mask out:  blue, red:"
            maskstr = raw_input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]),
                            float(maskstr.split(",")[1])]
            except ValueError:
                print 'Invalid entry. Enter wavelengths separated by commas'
            else:
                config_pars['mask_region1'] = maskwave
        elif option == 'm2':
            print "Enter wavelength window to mask out:  blue, red:"
            maskstr = raw_input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]),
                            float(maskstr.split(",")[1])]
            except ValueError:
                print 'Invalid entry. Enter wavelengths separated by commas'
            else:
                config_pars['mask_region2'] = maskwave
        elif option == 'm3':
            print "Enter wavelength window to mask out:  blue, red (Angstroms:"
            maskstr = raw_input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]),
                            float(maskstr.split(",")[1])]
            except ValueError:
                print 'Invalid entry. Enter wavelengths separated by commas'
            else:
                config_pars['mask_region3'] = maskwave

        # change the transition wavelength between the grisms
        elif option == 't':
            print "The current transition wavelength is: " + str(config_pars['transition_wave']) + "\nEnter the wavelength for the G102 to G141 transition:"
            try:
                config_pars['transition_wave'] = float(raw_input("> "))
            except ValueError:
                print 'Invalid entry. Enter wavelength of grism transition.'

        # change the nodes used for the continuum spline
        elif option == 'nodes':
            strnw = ','.join(str(nw) for nw in config_pars['node_wave'])
            print "Enter Wavelengths for Continuum Spline: w1, w2, w3, w4, ...."
            print "(current node wavelengths are: %s)" % strnw
            nodestr = raw_input("> ")
            nodesplit = nodestr.split(',')
            node_arr = []
            try:
                for nodelam in nodesplit:
                    node_arr.append(float(nodelam))
            except ValueError:
                print 'Invalid entry. Enter wavelengths separated by commas'
            else:
                node_arr = np.array(node_arr)
                # sort by wavelength
                node_arr = np.sort(node_arr)
                config_pars['node_wave'] = node_arr

        # reset all options
        elif option == 'reset':
            print "Reset configuration parameters, fwhm guess, and zguess to default values"
            config_pars = read_config('default.config')
            fwhm_guess = 2.35 * a_image * config_pars['dispersion_red']
            # reset strongest line, too
            index_of_strongest_line = 0
            lamline = lamlines_found[index_of_strongest_line]
            zguess = lamline / lam_Halpha - 1

        # change the blue cutoff of G102 (or whichever grism is present?)
        elif option == 'bluecut':
            print "The current blue cutoff is: " + str(config_pars['lambda_min']) + "\nChange the blue cutoff of G102:"
            try:
                config_pars['lambda_min'] = float(raw_input("> "))
            except ValueError:
                print 'Invalid entry. Enter wavelength of blue cutoff.'

        # change the red cutoff of G141
        elif option == 'redcut':
            print "The current red cutoff is: " + str(config_pars['lambda_max']) + "\nChage the red cutoff of G141:"
            try:
                config_pars['lambda_max'] = float(raw_input("> "))
            except ValueError:
                print 'Invalid entry. Enter wavelength of red cutoff.'

        # change to next brightest line
        elif option == 'n':
            nlines_found_cwt = np.size(lamlines_found)
            index_of_strongest_line = index_of_strongest_line + 1
            if index_of_strongest_line < (nlines_found_cwt):
                lamline = lamlines_found[index_of_strongest_line]
                zguess = lamline / 6564. - 1
            else:
                print 'There are no other automatically identified peaks. Select another option.'
                # stay at current line
                index_of_strongest_line -= 1

        # change to another line
        elif option == 'ha':
            zguess = (lamline / lam_Halpha) - 1
        elif option == 'hb':
            zguess = (lamline / lam_Hbeta) - 1
        elif option == 'o2':
            zguess = (lamline / lam_Oii) - 1
        elif option == 'o31':
            zguess = (lamline / lam_Oiii_1) - 1
        elif option == 'o32':
            zguess = (lamline / lam_Oiii_2) - 1
        elif option == 's2':
            zguess = (lamline / lam_Sii) - 1
        elif option == 's31':
            zguess = (lamline / lam_Siii_1) - 1
        elif option == 's32':
            zguess = (lamline / lam_Siii_2) - 1

        # note contamination
        elif option == 'contam':
            print "Specify contamination.\nEnter a comma-separated list of identifiers choosing from:\n  o2,hg,hb,o3,ha,s2,s31,s32,he1,c(ontinuum)"
            cf = raw_input("> ")
            cflags = [thing.strip() for thing in cf.split(',')]
            if 'c' in cflags:
                for k,v in contamflags.iteritems():
                    contamflags[k] = contamflags[k] | 2
            cflaglines = [thing for thing in cflags if thing != 'c']
            for contamflag in cflaglines:
                try:
                    contamflags[contamflag] = contamflags[contamflag] | 4
                except KeyError:
                    print '%s not known. Skipping'%contamflag

        # add a comment
        elif option == 'c':
            print "Enter your comment here:"
            comment = raw_input("> ")
            # sqlite3 database support - automatically creates and initializes DB if required
            databaseManager.saveAnnotation((par, obj, comment.decode('utf-8')))

        # set or unset one or more flags
        elif option == 'f':
            print('Enter a comma-separated list of flag, value pairs e.g. CONTAM, 1, CONTIN, 2:')
            print('Valid flags are {}'.format(WISPLFDatabaseManager.WISPLFDatabaseManager.validFlags))
            flagList = raw_input("> ")
            # sqlite3 database support - automatically creates and initializes DB if required
            databaseManager.setFlagsFromString(par, obj, flagList.decode('utf-8'))

        # print help message
        elif option == 'h':
            print_help_message()

        ### image/display options ###
        # change 2d stamp scaling to linear
        elif option == 'lin':
            if g102zeros is not None:
                show2dNEW('G102', par, obj, g102zeros, 'linear')
            if g141zeros is not None:
                show2dNEW('G141', par, obj, g141zeros, 'linear')

        # change 2d stamp scaling to log
        elif option == 'log':
            if g102zeros is not None:
                show2dNEW('G102', par, obj, g102zeros, 'log')
            if g141zeros is not None:
                show2dNEW('G141', par, obj, g141zeros, 'log')

        # change g102 2d stamp scaling to zscale
        elif option == 'zs102':
            print "Enter comma-separated range for G102 zscale: z1,z2"
            zscale = raw_input("> ")
            zs = zscale.split(',')
            try:
                z1 = float(zs[0])
                z2 = float(zs[1])
            except ValueError:
                print 'Invalid entry.'
            else:
                if g102zeros is not None:
                    show2dNEW('G102', par, obj, g102zeros, 'linear',
                              zran1=z1, zran2=z2)

        # change g141 2d stamp scaling to zscale
        elif option == 'zs141':
            print "Enter comma-separated range for G141 zscale: z1,z2"
            zscale = raw_input("> ")
            zs = zscale.split(',')
            try:
                z1 = float(zs[0])
                z2 = float(zs[1])
            except ValueError:
                print 'Invalid entry.'
            else:
                if g141zeros is not None:
                    show2dNEW('G141', par, obj, g141zeros, 'linear',
                              zran1=z1, zran2=z2)

        # recenter full images
        elif option == 'dc':
            showDirectNEW(obj)
            if show_dispersed:  # MB
                showDispersed(obj)

        # reload full iamges
        elif option == 'reload':
            showDirectNEW(obj, load_image=True)
            if show_dispersed:
                showDispersed(obj, load_image=True)

        # reload direct image region files
        elif option == 'dr':
            reloadReg()

        ### new options dealing with iterating objects ###
        # can't actually go back or choose another object now,
        # but allow for them now just in case
        elif option == 'b':
            print 'Please either reject or accept this object first.'
        elif 'obj' in option:
            print 'Please either reject or accept this object first.'
        # print remaining objects that have not yet been inspected
        elif option == 'left':
            print '    Remaining objects:'
            print remaining
        # print all objects in line list
        elif option == 'list':
            print '    All objects:'
            print allobjects

        # quit this object
        elif option == 'q':
            print 'Quitting Obj %i. Nothing saved to file' % (obj)
            print '-' * 72
            return 0

        # catch-all for everything else
        else:
            print "Invalid entry.  Try again."

        # print "OK"

    # only re-save data if the previous fit was discarded
    if not acceptPrevFit :
        # write to file if object was accepted
        if zset == 1 :
            # write object summary
            write_object_summary(par, obj, fitresults, snr_meas_array, 
                                 contamflags)
        
            # sqlite3 database support - automatically creates and initializes DB if required
            databaseManager.saveCatalogueEntry(databaseManager.layoutCatalogueData(par, obj, ra[0], dec[0], a_image[0],
                                                                                   b_image[0], jmag[0], hmag[0], fitresults, flagcont))

            writeToCatalog(linelistoutfile, par, obj, ra, dec, a_image,
                           b_image, jmag, hmag, fitresults, contamflags)

            writeFitdata(fitdatafilename, spec_lam, spec_val, spec_unc,
                         spec_con, spec_zer, full_fitmodel, full_contmodel, mask_flg)

            fitspec_pickle = open(fitdatafilename + '.pickle', 'wb')
            output_meta_data = [par, obj, ra, dec, a_image, b_image,
                                jmag, hmag, fitresults, flagcont, config_pars]
            pickle.dump(output_meta_data, fitspec_pickle)
            fitspec_pickle.close()
        else :
            # done == 1, but zset == 0 => rejected
            databaseManager.saveCatalogueEntry(databaseManager.layoutCatalogueData(par, obj, ra[0], dec[0], a_image[0],
                                                                                   b_image[0], jmag[0], hmag[0],
                                                                                   None,
                                                                                   None))
            databaseManager.setFlags(par, obj, [('REJECT', 1)])
            databaseManager.saveAnnotation((par, obj, 'REJECTED'))

        # write comments to file
        # if we go back to the previous objects, duplicate comments will still be
        # written
        writeComments(commentsfile, par, obj, comment)

        # write object to done file, incase process gets interrupted
        if not os.path.exists('linelist/done'):
            f = open('linelist/done', 'w')
        else:
            f = open('linelist/done', 'a')
        f.write('%i\n' % obj)
        f.close()
    # print closing line anyway
    print('-' * 72)


def check_input_objid(objlist, objid):
    """ """
    fulllist = ', '.join(['%i' % o for o in objlist])
    if objid not in objlist:
        print '   Obj %i is not in the line list' % (objid)
        print '   Full line list: \n%s' % (fulllist)
        o = raw_input('Please try again: > ')
        objid = int(re.search('\d+', o).group())
    return objid


def measure_z_interactive(linelistfile=" ", show_dispersed=True, use_stored_fit=False):
    #### STEP 0:   set ds9 window to tile mode ################################
    ###########################################################################
    # not the best way to do this, but matching the method in guis.py
    cmd = 'xpaset -p ds9 tile'
    os.system(cmd)
    cmd = 'xpaset -p ds9 tile grid layout 2 3'
    os.system(cmd)

    #### STEP 1:   get linelist ###############################################
    ###########################################################################
    if linelistfile == " ":
        files = glob('linelist/Par*lines.dat')
        if len(files) == 0:
            print 'No line list file found'
            return 0
        else:
            linelistfile = files[0]

    if not os.path.exists(linelistfile):
        print "Invalid path to line list file: %s" % (linelistfile)
        return 0
    else:
        print 'Found line list file: %s' % (linelistfile)

    #### STEP 2:   read the list of candidate lines  ####################
    ###########################################################################
    llin = asciitable.read(linelistfile, names=[
                           'parnos', 'grism', 'objid', 'wavelen', 'npix', 'ston'])
    parnos = llin['parnos']
    grism = llin['grism']
    objid = llin['objid']
    wavelen = llin['wavelen']
    npix = llin['npix']
    ston = llin['ston']
    objid_unique = np.unique(objid)

    #### STEP 3: define filenames and check for partially complete work #####
    #########################################################################
    if not os.path.exists('figs'):
        os.mkdir('figs')
    if not os.path.exists('fitdata'):
        os.mkdir('fitdata')
    parts = linelistfile.split('.')
    linelistoutfile = parts[0] + '_catalog.' + parts[1]
    commentsfile = parts[0] + '_comments.' + parts[1]
    # the file that will be used to determine which objects are "done"
    donefile = 'linelist/done'

    if os.path.isfile(linelistoutfile):
        print '\n Output file: \n  %s, \nalready exists\n' % linelistoutfile
        ask = raw_input('Append? [Y/n] ')
        if ask.lower() == 'n':
            os.unlink(linelistoutfile)
            os.unlink(commentsfile)
            # starting over, no objects have been done
            os.unlink(donefile)
            objid_done = np.array([])
        else:
            # an object may be written to the comment file before it has
            # actually been inspected, so use donefile for a list
            # of the "done" objects
            objid_done = np.atleast_1d(np.genfromtxt(donefile, dtype=int))
    else:
        if os.path.exists(donefile):
            os.unlink(donefile)
        objid_done = np.array([])

    #### STEP 4: create trace.reg files ############################
    #########################################################################
    trace102 = open('G102_trace.reg', 'w')
    trace102.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace102.write('wcs;\n')
    # sensitivity drops below 25% of max at wave < 8250 and wave > 11540
    # so box should be 3290 angstroms wide and be centered at 9895.
    trace102.write('box(9895,0,3290,1,1.62844e-12)\n')
    trace102.close()
    trace141 = open('G141_trace.reg', 'w')
    trace141.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace141.write('wcs;\n')
    # sensitivity drops below 25% of max at wave < 10917 and wave > 16904
    # so box should be 5897 angstroms wide and be centered at 13910.5
    trace141.write('box(13910.5,0,5897,1,0)\n')
    trace141.close()

    #### STEP 5:  Get zero and first order positions; unpack them ###########
    #########################################################################
    g102zeroordreg = "../DATA/DIRECT_GRISM/G102_0th.reg"
    g102firstordreg = "../DATA/DIRECT_GRISM/G102_1st.reg"
    g141zeroordreg = "../DATA/DIRECT_GRISM/G141_0th.reg"
    g141firstordreg = "../DATA/DIRECT_GRISM/G141_1st.reg"

    #### STEP 6:  Get object information from SExtractor catalog ############
    #########################################################################
    # a_image will give an order of magnitude estimate on the FWHM of the line
    #   this determines the initial guess and sets an upper limit on how broad
    #   it can be.
    # ra/dec, b_image, jmag, jerr, hmag, herr will be carried forward into
    #   the output linelist.
    # find all available cats
    secats = glob('../DATA/DIRECT_GRISM/fin_F*.cat')
    secats.sort()
    cat = asciitable.read(secats[0])
    beam = cat['col2']
    a_image = cat['col5']
    b_image = cat['col6']
    ra = cat['col8']
    dec = cat['col9']
    # which filter is this?
    if os.path.basename(secats[0]) == 'fin_F110.cat':
        jmag = cat['col13']
        jerr = cat['col14']
    else:
        jmag = np.ones(ra.shape, dtype=float) * 99.
        jerr = np.ones(ra.shape, dtype=float) * 99.
        hmag = cat['col13']
        herr = cat['col14']
    # read in second file if there are two
    if len(secats) == 2:
        cat2 = asciitable.read(secats[1])
        # second catalog should be hband
        hmag = cat2['col13']
        herr = cat2['col14']
    objtable = Table([beam, ra, dec, a_image, b_image, jmag, jerr, hmag, herr],
                     names=('obj', 'ra', 'dec', 'a_image', 'b_image', 'jmag',
                            'jerr', 'hmag', 'herr'))

    #### STEP 7:  Set up initial ds9 display ################################
    #########################################################################
    if os.path.exists(g102zeroordreg):
        g102zeroarr = getzeroorders(g102zeroordreg, g='G102')
        # nothing is done with the first orders anymore
        # g102firstarr=getfirstorders(g102firstordreg)
        show2dNEW('G102', parnos[0], objid_unique[0], g102zeroarr, 'linear')
    else:
        g102zeroarr = None
        g102firstarr = None
    if os.path.exists(g141zeroordreg):
        g141zeroarr = getzeroorders(g141zeroordreg, g='G102')
        # g141firstarr=getfirstorders(g141firstordreg)
        show2dNEW('G141', parnos[0], objid_unique[0], g102zeroarr, 'linear')
    else:
        g141zeroarr = None
        g141firstarr = None

    showDirectNEW(objid_unique[0], load_image=True)
    if show_dispersed:  # MB
        showDispersed(objid_unique[0], load_image=True)

    #### STEP 8:  Loop through objects ############
    #########################################################################
    remaining_objects = get_remaining_objects(objid_unique, objid_done)
    allobjects = [unique_obj for unique_obj in objid_unique]

    print '\nAs you loop through the objects, you can choose from the following\noptions at any time:\n\txxx = skip to object xxx\n\tb = revist the previous object\n\tleft = list all remaining objects that need review\n\tlist = list all objects in line list\n\tany other key = continue with the next object\n\tq = quit\n'

    while remaining_objects.shape[0] > 0:
        ndone = len(objid_done)
        progress = float(ndone) / float(len(objid_unique)) * 100.
        print "\nProgress: %.1f percent" % (progress)

        # do some things as long as there are still objects to inspect
        next_obj = remaining_objects[0]
        print 'Next up: Obj %i' % (next_obj)
        o = raw_input(
            "Enter 'obj xxx' to skip to Obj xxx or hit any key to continue. > ")

        if o == 'b':
            # need to figure out what object came before this one
            w = np.where(objid_unique == remaining_objects[0])
            # if on first object, this will roll around to previous object
            next_obj = objid_unique[w[0][0] - 1]
            print "Going back to previous object: Obj %i" % (next_obj)

        if o == 'left':
            #remaining_list = ', '.join(['%i'%i for i in remaining_objects])
            print '    Remaining objects:'
            print remaining_objects
            o = raw_input('> ')

        if o == 'list':
            print '    All objects:'
            print allobjects
            o = raw_input('> ')

        if o == 'q':
            print "Quitting; saved through previously completed object."
            return 0

        elif 'obj' in o:
            next_obj = int(re.search('\d+', o).group())
            # confirm that requested object is in line list
            next_obj = check_input_objid(objid_unique, next_obj)

        # pass the information for this object
        wlinelist = np.where(objid == next_obj)
        lamlines_found = wavelen[wlinelist]
        ston_found = ston[wlinelist]
        wcatalog = np.where(objtable['obj'] == next_obj)
        objinfo = objtable[wcatalog]
        inspect_object(parnos[0], next_obj, objinfo, lamlines_found, ston_found,
                       g102zeroarr, g141zeroarr, linelistoutfile, commentsfile,
                       remaining_objects, allobjects,
                       show_dispersed=show_dispersed)
        objid_done = np.append(objid_done, next_obj)
        remaining_objects = get_remaining_objects(objid_unique, objid_done)

    # outside the while loop, field is done
    redo = ' '
    while redo != 'q':
        print "You've finished this field.\nEnter an object ID below to revisit a particular object.\nOtherwise enter 'q' to quit the field."
        redo = raw_input("> ")
        try:
            next_obj = int(redo)
        except ValueError:
            print "Invalid entry. Enter an object ID or enter 'q' to quit"
        else:
            next_obj = check_input_objid(objid_unique, next_obj)
            # pass the information for this object
            wlinelist = np.where(objid == next_obj)
            lamlines_found = wavelen[wlinelist]
            ston_found = ston[wlinelist]
            wcatalog = np.where(objtable['Beam'] == next_obj)
            objinfo = objtable[wcatalog]
            inspect_object(parnos[0], next_obj, objinfo, lamlines_found,
                           ston_found, g102zeroarr, g141zeroarr,
                           linelistoutfile, commentsfile, remaining_objects,
                           allobjects, show_dispersed=show_dispersed)


# XXX this is never actually called????
# XXX        ### it may be desirable to overwrite the inital guesses, if we are trying to update the object.
# XXX        if (use_stored_fit ==True) & (os.path.exists('./fitdata/'+ fitdatafilename + '.pickle') == True):
# XXX            inpickle = './fitdata/' +fitdatafilename + '.pickle'
# XXX            fileObject = open(inpickle,'r')
# XXX            alldata = pickle.load(fileObject)
# XXX            fitresults_old = alldata[8]
# XXX            config_pars = alldata[10]
# XXX            zguess = fitresults_old['redshift']
# XXX            fwhm_guess = fitresults_old['fwhm_g141']

   # Clean up temp files
    if os.path.exists('./tempcoo.dat') == 1:
        os.unlink('./tempcoo.dat')
    if os.path.exists('./temp_zero_coords.coo') == 1:
        os.unlink('./temp_zero_coords.coo')
    if os.path.exists('./temp110.fits') == 1:
        os.unlink('./temp110.fits')
    if os.path.exists('./temp140.fits') == 1:
        os.unlink('./temp140.fits')
    if os.path.exists('./temp160.fits') == 1:
        os.unlink('./temp160.fits')
    if os.path.exists('./temp_zero_coords.reg') == 1:
        os.unlink('./temp_zero_coords.reg')
    if os.path.exists('G102_trace.reg') == True:
        os.unlink('G102_trace.reg')
    if os.path.exists('G141_trace.reg') == True:
        os.unlink('G141_trace.reg')


# parnos, objid are scalar not array.
def writeToCatalog(catalogname, parnos, objid, ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, contamflags):
    if not os.path.exists(catalogname):
        cat = open(catalogname, 'w')
        cat.write('#1  ParNo\n')
        cat.write('#2  ObjID\n')
        cat.write('#3 RA \n')
        cat.write('#4 Dec \n')
        cat.write('#5 Jmagnitude [99.0 denotes no detection] \n')
        cat.write('#6 Hmagnitude [99.0 denotes no detection] \n')
        cat.write('#7 A_IMAGE \n')
        cat.write('#8 B_IMAGE \n')
        cat.write('#9 redshift \n')
        cat.write('#10 redshift_err \n')
        cat.write('#11 dz_oiii \n')
        cat.write('#12 dz_oii \n')
        cat.write('#13 dz_siii_he1 \n')
        cat.write('#14 G141_FWHM_Obs  [Angs] \n')
        cat.write('#15 G141_FWHM_Obs_err  \n')
        cat.write('#16 oii_flux \n')
        cat.write('#17 oii_error \n')
        cat.write('#18 oii_EW_obs \n')
        cat.write('#19 oii_contam \n')
        cat.write('#20 hg_flux \n')
        cat.write('#21 hg_err \n')
        cat.write('#22 hg_EW_obs \n')
        cat.write('#22 hg_contam \n')
        cat.write('#23 hb_flux \n')
        cat.write('#24 hb_err \n')
        cat.write('#25 hb_EW_obs \n')
        cat.write('#26 hb_contam \n')
        cat.write('#27 oiii_flux [both lines] \n')
        cat.write('#28 oiii_err [both lines] \n')
        cat.write('#29 oiii_EW_obs [both lines] \n')
        cat.write('#30 oiii_contam [both lines] \n')
        cat.write('#31 hanii_flux \n')
        cat.write('#32 hanii_err \n')
        cat.write('#33 hanii_EW_obs \n')
        cat.write('#34 hanii_contam \n')
        cat.write('#35 sii_flux \n')
        cat.write('#36 sii_err \n')
        cat.write('#37 sii_EW_obs \n')
        cat.write('#38 sii_contam \n')
        cat.write('#39 siii_9069_flux \n')
        cat.write('#40 siii_9069_err \n')
        cat.write('#41 siii_9069_EW_obs \n')
        cat.write('#42 siii_9069_contam \n')
        cat.write('#43 siii_9532_flux \n')
        cat.write('#44 siii_9532_err \n')
        cat.write('#45 siii_9532_EW_obs \n')
        cat.write('#46 siii_9532_contam \n')
        cat.write('#47 he1_10830_flux \n')
        cat.write('#48 he1_10830_err \n')
        cat.write('#49 he1_10830_EW_obs \n')
        cat.write('#50 he1_10830_contam \n')

#        cat.write('#43 ContamFlag \n')
        cat.close()
       # cat.write('#45 Comment \n')

#    else:
#        cat = open(catalogname, 'a')

    outstr = '{:<8d}'.format(parnos) + \
        '{:<6d}'.format(objid) +\
        '{:<12.6f}'.format(ra_obj[0]) + \
        '{:<12.6f}'.format(dec_obj[0]) + \
        '{:<8.2f}'.format(jmag_obj[0]) + \
        '{:<8.2f}'.format(hmag_obj[0]) + \
        '{:<8.3f}'.format(a_image_obj[0]) + \
        '{:<8.3f}'.format(b_image_obj[0]) + \
        '{:>8.4f}'.format(fitresults['redshift']) + \
        '{:>10.4f}'.format(fitresults['redshift_err']) +\
        '{:>10.4f}'.format(fitresults['dz_oiii'])  + \
        '{:>10.4f}'.format(fitresults['dz_oii'])   + \
        '{:>10.4f}'.format(fitresults['dz_siii_he1']) +\
        '{:>10.2f}'.format(fitresults['fwhm_g141']) + \
        '{:>10.2f}'.format(fitresults['fwhm_g141_err'])  +  \
        '{:>13.2e}'.format(fitresults['oii_flux'])  + \
        '{:>13.2e}'.format(fitresults['oii_error']) +  \
        '{:>13.2e}'.format(fitresults['oii_ew_obs']) +\
        '{:>6d}'.format(contamflags['o2']) +\
        '{:>13.2e}'.format(fitresults['hg_flux']) +\
        '{:>13.2e}'.format(fitresults['hg_error']) + \
        '{:>13.2e}'.format(fitresults['hg_ew_obs']) +\
        '{:>6d}'.format(contamflags['hg']) +\
        '{:>13.2e}'.format(fitresults['hb_flux']) + \
        '{:>13.2e}'.format(fitresults['hb_error']) + \
        '{:>13.2e}'.format(fitresults['hb_ew_obs']) +\
        '{:>6d}'.format(contamflags['hb']) +\
        '{:>13.2e}'.format(fitresults['oiii_flux']) + \
        '{:>13.2e}'.format(fitresults['oiii_error']) + \
        '{:>13.2e}'.format(fitresults['oiii_ew_obs']) +\
        '{:>6d}'.format(contamflags['o3']) +\
        '{:>13.2e}'.format(fitresults['hanii_flux']) + \
        '{:>13.2e}'.format(fitresults['hanii_error']) + \
        '{:>13.2e}'.format(fitresults['hanii_ew_obs']) + \
        '{:>6d}'.format(contamflags['ha']) +\
        '{:>13.2e}'.format(fitresults['sii_flux']) + \
        '{:>13.2e}'.format(fitresults['sii_error']) + \
        '{:>13.2e}'.format(fitresults['sii_ew_obs']) +\
        '{:>6d}'.format(contamflags['s2']) +\
        '{:>13.2e}'.format(fitresults['siii_9069_flux']) + \
        '{:>13.2e}'.format(fitresults['siii_9069_error']) + \
        '{:>13.2e}'.format(fitresults['siii_9069_ew_obs']) +\
        '{:>6d}'.format(contamflags['s31']) +\
        '{:>13.2e}'.format(fitresults['siii_9532_flux']) + \
        '{:>13.2e}'.format(fitresults['siii_9532_error']) + \
        '{:>13.2e}'.format(fitresults['siii_9532_ew_obs']) +\
        '{:>6d}'.format(contamflags['s32']) +\
        '{:>13.2e}'.format(fitresults['he1_flux']) + \
        '{:>13.2e}'.format(fitresults['he1_error']) + \
        '{:>12.2e}'.format(fitresults['he1_ew_obs']) +\
        '{:>6d}'.format(contamflags['he1']) + '\n'
        #     '   ' + '{:<6d}'.format(flagcont) + ' \n'

#    """
    # if a row already exists for this object, comment it out
    objstr = '{:<8d}'.format(parnos) + '{:<6d}'.format(objid)
    for line in fileinput.input(catalogname, inplace=True):
        if objstr in line:
            print "#%s" % line,
        else:
            print '%s' % line,
#    """

    cat = open(catalogname, 'a')
    cat.write(outstr)
    cat.close()


def writeFitdata(filename, lam, flux, eflux, contam, zero, fit, continuum, masks):
    """ """
    fitspec_file = filename + '.dat'
    t = Table([lam, flux, eflux, contam, zero, fit, continuum, masks],
              names=('Lam', 'Flam', 'Flam_err', 'Contam', 'Zero', 'Fitmodel', 'Contmodel', 'Masked'))
    t['Lam'].format = '{:.1f}'
    t['Flam'].format = '{:.5e}'
    t['Flam_err'].format = '{:.5e}'
    t['Contam'].format = '{:.5e}'
    t['Zero'].format = '{:.0f}'
    t['Fitmodel'].format = '{:.5e}'
    t['Contmodel'].format = '{:.5e}'
    t['Masked'].format = '{:.0f}'
    asciitable.write(t, fitspec_file, fill_values=[(asciitable.masked, '--')],
                     overwrite=True, format='fixed_width_two_line',
                     position_char='#', delimiter_pad=' ')


def writeComments(filename, parnos, objid, comment):
    if os.path.exists(filename) == False:
        cat = open(filename, 'w')
    else:
        cat = open(filename, 'a')

    outstr = '{:<8d}'.format(parnos) + \
        '{:<6d}'.format(objid) +\
        comment + '\n'

    cat.write(outstr)
    cat.close()


def UpdateCatalog(linelistoutfile):

    allDirectoryFiles = os.listdir('./fitdata/')
    objid_list = []
    for obj in allDirectoryFiles:
        x = obj.split('_')[2]
        objid_list.append(int(x))
        Parno = obj.split('_')[0]   # this is inefficient, but I don't care.
    objid_list = np.sort(np.unique(np.array(objid_list)))
    for obj in objid_list:
        print obj
        inpickle = './fitdata/' + Parno + \
            '_BEAM_' + str(obj) + '_fitspec.pickle'
        fileObject = open(inpickle, 'r')
        alldata = pickle.load(fileObject)
        # definition from above
       #                      0          1                 2      3        4           5            6         7         8           9         10
       # output_meta_data = [parnos[0], objid_unique[i], ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont, config_pars]
        parnos = alldata[0]
        objid_unique = alldata[1]
        ra_obj = alldata[2]
        dec_obj = alldata[3]
        a_image_obj = alldata[4]
        b_image_obj = alldata[5]
        jmag_obj = alldata[6]
        hmag_obj = alldata[7]
        fitresults = alldata[8]
        flagcont = alldata[9]
        # config_pars = alldata[10] ## not used here.

        WriteToCatalog(linelistoutfile, parnos, objid_unique, ra_obj, dec_obj,
                       a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont)
