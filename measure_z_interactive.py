#!/usr/bin/env python
##################################################################################
##################################################################################
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
##################################################################################
import os
import distutils
#import numpy as np
import scipy
import pylab as plt
from scipy.interpolate import spline
from astropy.table import Table
from distutils.sysconfig import *    ### question-- what is this for? 
import sys
from matplotlib import gridspec
import matplotlib.transforms as mtransforms

from wisp_analysis import *



def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False



def getzeroorders (zeroorderpath,g='G141',magcut=23.5): # MB: changed from 23.0
    zop=open(zeroorderpath,'r')
    zox=[]
    zoy=[]
    zoid=[]
    zmag=[]
    for line in zop:
        if len(line)>60:
            linesplit = line.split()
            zox.append(float(linesplit[1][0:-1]))
            zoy.append(float(linesplit[2][0:-3])) 
            zoid.append(int(linesplit[-2].split('{')[-1]))
            zmag.append(float(linesplit[-1][1:-2]))  ### get mag from reg file

    zop.close()
    zoid=np.array(zoid)
    zoy=np.array(zoy)
    zox=np.array(zox)
    zmag=np.array(zmag)
    cond = (zmag <= magcut)
    return zox[cond],zoy[cond],zoid[cond]


def getfirstorders (firstorderpath):
    fop=open(firstorderpath,'r')
    fox=[]
    foy=[]
    folen=[]
    fowid=[]
    foid=[]
    for line in fop:
        #if line[0]!='#':
        linesplit = line.split() 
        fox.append(float( linesplit[1][0:-1] ))  ### [0:-1] strips off the comma. 
        foy.append(float(linesplit[2][0:-1]))
        folen.append(float(linesplit[3][0:-1])) 
        fowid.append(float(linesplit[-1].split('{')[-1].split('}')[0]))  ### python is weird.
    return (fox,foy,folen,fowid,foid)


def plot_broken_spec(wavelength, flux, ax, **kwargs):
    """ """
    diff = np.diff(wavelength)
    # define the gap size as a multiple of the dispersion, using G141
    gapsize = 2*50  # = ~2 pix in G141, ~4 pix in G102
    split_i = np.where(diff > gapsize)[0] + 1
    for _w,_f in zip(np.split(wavelength,split_i),np.split(flux,split_i)):
        ax.plot(_w, _f, **kwargs)


def write_fitdata(filename, lam, flux, eflux, contam, zero, fit, continuum, masks):
    """ """
    fitspec_file = filename+'.dat'
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
    asciitable.write(t, fitspec_file, fill_values=[(asciitable.masked,'--')], 
                     overwrite=True, format='fixed_width_two_line', 
                     position_char='#', delimiter_pad=' ')


def measure_z_interactive (linelistfile=" ", show_dispersed=True, use_stored_fit = False):
    #### STEP 0:   set ds9 window to tile mode ################################
    ###########################################################################
    ### not the best way to do this, but matching the method in guis.py
    cmd='xpaset -p ds9 tile'
    os.system(cmd)
    cmd='xpaset -p ds9 tile grid layout 2 3'
    os.system(cmd)


    #### STEP 1:   get linelist ###############################################
    ###########################################################################
    if linelistfile==" ":
        ### check the .dat files to find the ParXXX, which is in the name of the linelistfile 
        allDirectoryFiles=os.listdir('.')
        for files in allDirectoryFiles:
            if files[0:3]=='Par':
                llpts=files.split('_')
                linelistfile='linelist/'+llpts[0]+'lines.dat'
                break

    #### throw an error or print that the linelistfile is found         
    if linelistfile==" " or os.path.exists(linelistfile)==0:
        print "Invalid path to line list file: %s" % (linelistfile)
        return 0
    else:
        print "Found line list file %s" % (linelistfile)

    #### STEP 2:   read the list of candidate lines  ####################
    ###########################################################################
    llin = asciitable.read(linelistfile, names = ['parnos', 'grism', 'objid', 'wavelen', 'npix', 'ston']) 
    parnos = llin['parnos']
    grism = llin['grism']
    objid = llin['objid'] 
    wavelen = llin['wavelen'] 
    npix = llin['npix'] 
    ston = llin['ston'] 
    #flag  = llin['flag']   ## these flags are defined line-by-line and are not useful in the object-by-object context 
    nstart =0
    objid_unique = np.unique(objid) 

    #parnos,grism,objid,wavelen,npix,ston,flag,nstart,setzs,flagcont,comment=readll(linelistfile,recover_temp)
    #### STEP 3: define other filenames ############################
    #########################################################################
    ##### these are for recovering/storing partially complete work

    parts=linelistfile.split('.')
    linelistoutfile=parts[0]+'_catalog.'+parts[1]
    commentsfile = parts[0]+'_comments.'+parts[1]
    
    if os.path.isfile(linelistoutfile):
        print '\n Output file: \n  %s, \nalready exists\n' % linelistoutfile
        ask = raw_input('Append? [y/n] ')
        if ask.lower() == 'n':
           os.unlink(linelistoutfile) 
           os.unlink(commentsfile) 
        else : 
           ### because the comments are not machine readable, we have to do this the hard way.  
           f = open(commentsfile) 
           objid_done = [] 
           for line in f:
               x = line.split() 
               objid_done.append(float(x[1])) 
           objid_done = np.max(np.array(objid_done))
           w=np.where(objid_unique > objid_done) 
           #### note that objid_done only includes good objects. therefore 
           #w will return the last n objects that you decided were crap.  Sorry about that.
           try:
               nstart = w[0][0]
           except IndexError:
           #if np.size(nstart) == 0:  ## w is an empty list if you have already finished the field
               print "You've finished this field." 
               return 0 
           else:   
               ############################################################
               #### what the heck does this do?   
               #### this reloads the direct and grisms images if you are 
               #### restarting a partially completed field. 
               #### otherwise showDirectNEW only loads the images the first
               #### time; it adds time to reload the images for every object
               showDirectNEW(1,0)
               if show_dispersed:  # MB
                   showDispersed(1,0)
   
   
    #### STEP 4: create trace.reg files ############################
    #########################################################################
    trace102 = open('G102_trace.reg', 'w')
    trace102.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace102.write('wcs;\n')
    trace102.write('box(9950,0,3100,1,1.62844e-12)\n')
    trace102.close()
    trace141 = open('G141_trace.reg', 'w')
    trace141.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace141.write('wcs;\n')
    trace141.write('box(14142.49,0,5500,1,0)\n')
    trace141.close()



    ##### STEP 5:  defined some wavelengths; these will be used later for plotting.##### 
    ####################################################################################
    lam_Halpha=6563.0
    lam_Hbeta=4861.0
    lam_Hg = 4341.0 
    lam_Oiii_1=4959.0
    lam_Oiii_2=5007.0
    lam_Oii=3727.0
    lam_Sii=6724.0
    lam_Siii_1=9069.0
    lam_Siii_2=9532.0
    #lam_Lya=1216.0
    lam_He=10830.0
    #lam_Fe=12600.0
    #lam_Pag=10940.0
    #lam_Pab=12810.0
    
    suplines=[lam_Oii, lam_Hg, lam_Hbeta,lam_Oiii_2,lam_Halpha,lam_Sii,lam_Siii_1,lam_Siii_2,lam_He]
    suplines_str = ['[OII]', r'H$\gamma$', r'H$\beta$', '[OIII]', r'H$\alpha$', '[SII]', '[SIII]', '[SIII]', 'HeI']  
    
    #### STEP 6:  Get zero and first order positions; unpack them #############
    #########################################################################
    g102zeroarr=[]
    g102firstarr=[]
    g102zeroordreg="../DATA/DIRECT_GRISM/G102_0th.reg"
    g102firstordreg="../DATA/DIRECT_GRISM/G102_1st.reg"
    if os.path.exists(g102zeroordreg)==1:
        g102zeroarr=getzeroorders(g102zeroordreg,g='G102')
        g102firstarr=getfirstorders(g102firstordreg)

    g141zeroarr=[]
    g141firstarr=[]
    g141zeroordreg="../DATA/DIRECT_GRISM/G141_0th.reg"
    g141firstordreg="../DATA/DIRECT_GRISM/G141_1st.reg"
    if os.path.exists(g141zeroordreg)==1:
        g141zeroarr=getzeroorders(g141zeroordreg)
        g141firstarr=getfirstorders(g141firstordreg)
        
    if len(g102zeroarr)>0:
        g102zerox=g102zeroarr[0]
        g102zeroy=g102zeroarr[1]
        g102zeroid=g102zeroarr[2]
        g102firstx=g102firstarr[0]
        g102firsty=g102firstarr[1]
        g102firstlen=g102firstarr[2]
        g102firstwid=g102firstarr[3]
        g102firstid=g102firstarr[4]
    else:
        g102zerox=[]
        g102zeroy=[]
        g102zeroid=[]
        g102firstx=[]
        g102firsty=[]
        g102firstlen=[]
        g102firstwid=[]
        g102firstid=[]
        
    if len(g141zeroarr)>0:
        g141zerox=g141zeroarr[0]
        g141zeroy=g141zeroarr[1]
        g141zeroid=g141zeroarr[2]
        g141firstx=g141firstarr[0]
        g141firsty=g141firstarr[1]
        g141firstlen=g141firstarr[2]
        g141firstwid=g141firstarr[3]
        g141firstid=g141firstarr[4]
    else:
        g141zerox=[]
        g141zeroy=[]
        g141zeroid=[]
        g141firstx=[]
        g141firsty=[]
        g141firstlen=[]
        g141firstwid=[]
        g141firstid=[]
        


    ###### STEP 7 ################################################################## 
    ##### get the SEXtractor data on a_image, which gives an order of 
    ##### magnitude estimate on the FWHM of the line. this determines the initial guess and sets an upper limit on how broad it can be.
    
    #### ALSO GET ra/dec, b_image, jmag, jerr, hmag, herr because these will be carried forward into the output linelist.

    secat ="../DATA/DIRECT_GRISM/fin_F110.cat"
    setab = asciitable.read(secat) 
    beam_list = setab['col2']
    a_image_list = setab['col5'] 
    b_image_list = setab['col6'] 
    ra_list = setab['col8'] 
    dec_list =setab['col9'] 
    jmag_list = setab['col13']
    jerr_list = setab['col14'] 
    secat2 = "../DATA/DIRECT_GRISM/fin_F160.cat" 
    setab2 = asciitable.read(secat2) 
    hmag_list = setab2['col13'] 
    herr_list = setab2['col14'] 



    #### define some things... 
    #progress=0.0
    i=nstart   

    ###### STEP 8 ########################################################### 
    ### while loop through line/objects that need fitting####################
    while i<len(objid_unique):
        #### note progress. 
        progress=float(i)/float(len(objid_unique))*100.0
        print "Progress: %.1f percent" % (progress)
        
        ### get spectrum files for this particular object
       
        specname='Par' + str(parnos[0]) + '_BEAM_' + str(objid_unique[i]) + 'A.dat'
        specnameg102='Par' + str(parnos[0]) + '_G102_BEAM_' + str(objid_unique[i]) + 'A.dat'   ### may not exist.
        specnameg141='Par' + str(parnos[0]) + '_G141_BEAM_' + str(objid_unique[i]) + 'A.dat'
        plotTitle='Par' + str(parnos[0]) + '_BEAM_' + str(objid_unique[i])
        if os.path.exists('figs') == False: 
            os.mkdir('figs') 
        plotfilename = 'figs/'+plotTitle + '_fit.png'
        if os.path.exists('fitdata') == False: 
            os.mkdir('fitdata') 
        fitdatafilename = 'fitdata/'  +plotTitle + '_fitspec'

         
        ##### also  start with a fresh set of config pars... we may change these during the interactive fitting. 
        config_pars = read_config('default.config')

        ### get line, fwhm, and z estimate
        ### choose the lamline that has the highest s/n estimate 
        w=np.where(objid == objid_unique[i]) 
        w=w[0]

        lamlines_found = wavelen[w] 
        ston_found = ston[w] 

        s = np.argsort(ston_found) 
        ston_found = ston_found[s[::-1]]  ### reverse s/n order. 
        lamlines_found = lamlines_found[s[::-1]]  
        
        index_of_strongest_line  = 0 
        lamline = lamlines_found[index_of_strongest_line] 
        zguess = lamline /6564. - 1

        w=np.where(beam_list == objid_unique[i])
        w=w[0][0]
        fwhm_guess = 2.35 *  a_image_list[w] * config_pars['dispersion_red'] ### fwhm is defined for the red side, regardless of where line is.
        zguess = lamline /6564. - 1   ### since we are working with the strongest line only, assume it is Ha. 
        
        ### also set ra, dec, etc. for this object
        ra_obj = ra_list[w]
        dec_obj = dec_list[w] 
        a_image_obj = a_image_list[w] 
        b_image_obj = b_image_list[w] 
        jmag_obj = jmag_list[w] 
        jerr_obj = jerr_list[w] 
        hmag_obj= hmag_list[w] 
        herr_obj = herr_list[w] 



        ### it may be desirable to overwrite the inital guesses, if we are trying to update the object. 
        if (use_stored_fit ==True) & (os.path.exists('./fitdata/'+ fitdatafilename + '.pickle') == True): 
            inpickle = './fitdata/' +fitdatafilename + '.pickle' 
            fileObject = open(inpickle,'r') 
            alldata = pickle.load(fileObject)
            fitresults_old = alldata[8] 
            config_pars = alldata[10]
            zguess = fitresults_old['redshift'] 
            fwhm_guess = fitresults_old['fwhm_g141'] 
                


        #### show zero orders.   
        if len(g102zerox)>0:
            show2dNEW('G102',parnos[i],int(objid_unique[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'linear')

        if len(g141zerox)>0:
            show2dNEW('G141',parnos[i],int(objid_unique[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'linear')
            showDirectNEW(objid_unique[i],i)
            if show_dispersed:  # MB
               showDispersed(objid_unique[i],i)    

       


         ### flag to determine whether we've measured the z or not. 
        next=0
        comment = ' '  ### in case comment isn't filled in, we can still write it.

        while (next==0):

            ### do this every time because it is fast, and sometimes we re-read with mask or different transition wave.
            #### read spectra and unpack  
            if os.path.exists(specnameg102) == True: 
               tab_blue = asciitable.read(specnameg102, names = ['lambda',  'flux', 'ferror', 'contam', 'zero']) 
               tab_red =  asciitable.read(specnameg141,  names = ['lambda',  'flux', 'ferror', 'contam', 'zero'])
               spdata = trim_spec(tab_blue, tab_red, config_pars, mask_zeros=True, return_masks=True)    ### this is a from the wisp package.
            else: 
               tab_red =  asciitable.read(specnameg141,  names = ['lambda',  'flux', 'ferror', 'contam', 'zero'])
               spdata = trim_spec(None, tab_red, config_pars, mask_zeros=True, return_masks=True) 
    
            spec_lam = spdata[0]
            spec_val = spdata[1]
            spec_unc = spdata[2] 
            spec_con = spdata[3] 
            spec_zer = spdata[4]
            mask_flg = spdata[5]

            ### get the zeroth orders
            w=np.where(spec_zer == 3) 
            spec_zero_bad = spec_zer * 0 -1 
            spec_zero_bad[w] = 1.
            ### mild zeroth orders
            w=np.where( spec_zer == 2)  
            spec_zero_mild = spec_zer * 0 -1
            spec_zero_mild[w] = 1.

            # apply the mask to the wavelength array
            masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)
            # compress the masked arrays for fitting
            fit_inputs = [np.ma.compressed(masked_spec_lam), np.ma.compressed(spec_val), np.ma.compressed(spec_unc), config_pars, zguess, fwhm_guess, str(objid_unique[i])]
            fitresults = fit_obj(fit_inputs)  ### parsing the inputs this way facilitates parallel processing when fitting is done in batch mode.
            zfit = fitresults['redshift']
            print "Guess Redshift z =  %f" % (zguess) 
            print "Fit Redshift z =  %f" % (zfit)
            lamobs = (1 + zguess) * np.array(suplines) 
            fitpars = fitresults['fit_parameters'] 
            fitpars_nolines = cp.deepcopy(fitpars) 
            fitpars_nolines[9:19] = 0. 
            fitpars_nolines[11] = 1.4  ### can't kill this one or divide by zero. 
            fitpars_nolines[12] = 0.1
            fitmodel = emissionline_model(fitpars, np.ma.compressed(masked_spec_lam)) * fitresults['fit_scale_factor']  
            contmodel = emissionline_model(fitpars_nolines, np.ma.compressed(masked_spec_lam)) * fitresults['fit_scale_factor']  
            # the fitting is done on compressed arrays, so we need to 
            # create masked versions of the fit and continuum models
            full_fitmodel = np.zeros(spec_lam.shape, dtype=float)
            full_contmodel = np.zeros(spec_lam.shape, dtype=float)
            full_fitmodel[np.ma.nonzero(spec_val)] = fitmodel
            full_contmodel[np.ma.nonzero(spec_val)] = contmodel
            full_fitmodel = np.ma.masked_where(np.ma.getmask(spec_val), full_fitmodel)
            full_contmodel = np.ma.masked_where(np.ma.getmask(spec_val), full_contmodel)

            snr_meas_array = np.array( [ fitresults['oii_flux']/fitresults['oii_error'], fitresults['hg_flux']/fitresults['hg_error'], fitresults['hb_flux']/fitresults['hb_error'], 
                fitresults['oiii_flux']/fitresults['oiii_error'], fitresults['hanii_flux']/fitresults['hanii_error'], fitresults['sii_flux']/fitresults['sii_error'], 
                fitresults['siii_9069_flux']/fitresults['siii_9069_error'], fitresults['siii_9532_flux']/fitresults['siii_9532_error'], fitresults['he1_flux']/fitresults['he1_error']])

            plt.ion()
            #plt.figure(1,figsize=(11,8))
            fig = plt.figure(1,figsize=(11,8))
            plt.clf()
            gs = gridspec.GridSpec(3,4)
            ax1 = fig.add_subplot(gs[0:2,:])
            ax2 = fig.add_subplot(gs[2:,:])

            xmin=np.ma.min(spec_lam)-200.0
            xmax=np.ma.max(spec_lam)+200.0
            ymin=np.ma.min(spec_val)
            ymax=1.5*np.ma.max(spec_val)

            ax1.plot(spec_lam, spec_val, 'k',spec_lam, spec_con,'r', ls='steps')
            ax1.axvline(x=config_pars['transition_wave'], c='c',linestyle=':', lw=3)

            ax1trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
            ax2trans = mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
            ### plot observed wavelengths of all the possible lines. 
            for li,lstring, sn_meas in zip(lamobs, suplines_str, snr_meas_array): 
                if (li > xmin+100) & (li < xmax - 100) : 
                    for ax in [ax1,ax2]:
                        ax.axvline(x=li, color ='b')
                    stringplot = lstring + '   (' + str(round(sn_meas, 2)) + ')'
                    # use data coordinates for x-axis and axes coords for y-axis
                    ax1.text(li, 0.85, stringplot, rotation='vertical', ha='right', fontsize='16', transform=ax1trans)
            
            #plt.axvline(x=lamline, color = 'r', lw=2)
            ax1.plot(spec_lam, full_fitmodel, color ='r', lw=1.5)
            ax1.plot(spec_lam, full_contmodel, color = 'b', linestyle = '--', lw=1.5)
            ### plot 0th orders
            for ax in [ax1,ax2]:
                # use data coordinates for x-axis and axes coords for y-axis
                trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
                ax.fill_between(spec_lam, 0, 1, where=spec_zero_bad==1, color = 'red', alpha=0.3, transform=trans, label='Major 0th order contam')
                ax.fill_between(spec_lam, 0, 1, where=spec_zero_mild==1, color = 'orange', alpha=0.3, transform=trans, label='Minor 0th order contam')

            ### plot any masked regions
            for mr,label in zip(['mask_region1','mask_region2','mask_region3'],['masked regions',None,None]):
                if (config_pars[mr][0] != 0.) & (config_pars[mr][1] != 0.):
                    for ax in [ax1,ax2]:
                        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
                        ax.fill_between(config_pars[mr], 0, 1, color = 'grey', alpha=0.3, transform=trans, label=label)
            ax1.legend(bbox_to_anchor=[1.05,1.15])

            ### find values of spec_lam nearest to the nodes 
            nodelam = config_pars['node_wave']  
            nl_arr = [] 
            cont_node = [] 
            for nl in nodelam: 
                #w=np.where(np.abs(spec_lam - nl)  == np.min(np.abs(spec_lam - nl)))
                #w=w[0][0]
                w = np.argmin(np.abs(spec_lam - nl))
                nl_arr.append(spec_lam[w]) 
                cont_node.append(full_contmodel[w]) 
            ax1.plot(nl_arr, cont_node, 'ko', ms=9)

            #### repeat for line_candidates  
            lf_lam  = [] 
            lf_cont = [] 
            for lf in lamlines_found:  
                #w=np.where(np.abs(spec_lam - lf)  == np.min(np.abs(spec_lam -lf)))
                #w=w[0][0]
                w = np.argmin(np.abs(spec_lam - lf))
                lf_lam.append(spec_lam[w]) 
                lf_cont.append(full_contmodel[w]) 
            ax1.plot(lf_lam, lf_cont, 'bo', ms=9)

            #### indicate "current" line
            current_lam = lamlines_found[index_of_strongest_line]
            current_cont = contmodel[np.argmin(np.abs(np.ma.compressed(masked_spec_lam)-current_lam))]
            
            ax1.plot(current_lam, current_cont, 'ro', ms=10)

            ax1.set_ylabel(r'F$_\lambda$ ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$', size='xx-large')
            ax1.set_xlim([xmin, xmax])
            ax1.set_ylim([ymin, ymax])
            ax1.set_title(plotTitle)
 
            #### second panel for s/n
            s2n=(spec_val-full_contmodel)/spec_unc
            s2n_lam=spec_lam
            mask=np.logical_and(s2n>-10000., s2n<10000.)
            s2n=s2n[mask]
            s2n_lam=s2n_lam[mask]
            ax2.plot(s2n_lam,s2n,'k-',linestyle='steps')
            ymin=s2n.min()
            ymax=1.5*s2n.max()
            ax2.axhline(y=config_pars['n_sigma_above_cont'], c='r')
            for li in lamobs :
                ax2.axvline(x=li, color ='b')
            ax2.axvline(x=config_pars['transition_wave'], c='c',linestyle=':', lw=3)
            ax2.set_xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
            ax2.set_ylabel(r'S/N',size='xx-large')
            ax2.set_xlim([xmin, xmax])
            ax2.set_ylim(ymin, ymax)
            #fig = plt.gcf() 
            fig.savefig(plotfilename)   ### saves the figure for everything; junk objects and all;  will repeat/overwrite while iterating on the interactive fit. 
            plt.draw()
#            plt.draw()   ### why is this here twice??? 
            



            ##### options: 
            print "Enter option (read carefully, options have changed): \n \
             \t b=back \n \
             \t a=accept object fit \n \
             \t ac=accept object fit, noting contamination\n  \
             \t r=reject object \n \
             \t z = enter a different z guess  \n \
             \t w = enter a different emission line wavelength guess  \n \
             \t ha,  or hb, o31, o32, o2, s2, s31, s32 =  change redshift guess \n \
             \t n = skip to next brightest line found in this object \n \
             \t fw = change the fwhm guess in pixels \n \
             \t c=add comment \n \
             \t t=change transition wavelength \n \
             \t m1, m2, or m3 =mask up to three discontinuous wavelength regions \n \
             \t nodes = change the wavelengths for the continuum spline nodes \n \
             \t bluecut = change the blue cutoff of the G102 spec \n \
             \t redcut  = change the red cutoff of the G141 spec \n \
             \t reset = reset interactive options back default for this object \n \
             \t lin=linear z-scale \n \
             \t log=logarithmic  \n \
             \t zs102=z1,z2 comma-separated range for G102 zscale \n  \
             \t zs141=z1,z2 comma-separated range for G141 zscale \n  \
             \t dc=recenter direct images \n \
             \t dr=reload direct image reg files\n \
             \t q=quit\n"
            option = raw_input(">")
            if option=='r':
                next=1
                zset =0 
#            elif option=='rc':
#                next=1
#                flagcont = 3  
#                zset=0
            elif option=='a':
                next=1
                zset=1  ### the redshift is stored in "fitresults" 
                flagcont=1 #MR   

            elif option=='ac':
                next=1
                zset=1
                flagcont = 2 


            #### these options keep you iterating in the while loop    
            elif option=='z':
                print "The current redshift guess is: %f\nEnter Redshift Guess:" % zguess 
                zguess = float(raw_input(">")) 
            elif option=='w':
                print "The current emission line wavelength is: %f\nEnter Wavelength Guess in Angstroms:" % lamline
                zguess = float(raw_input(">")) / 6564. - 1
            elif option =='fw': 
                print "Enter a Guess for FWHM in pixels"
                print "The current fwhm_fit is:  " + str(fitresults['fwhm_g141'] /config_pars['dispersion_red']) + " and 2*A_image is: " + str(2 * a_image_obj) 
                fwhm_guess = config_pars['dispersion_red']  * float(raw_input(">"))
	    elif option=='c': # Changed by NR version 0.2.5
                print "Enter your comment here:"
                comment =raw_input(">")
            elif option == 'm1': 
                print "Enter wavelength window to mask out:  blue, red:" 
                maskstr = raw_input(">") 
                try:
                    maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
                except ValueError:
                    print 'Invalid entry. Enter wavelengths separated by commas'
                else:
                    config_pars['mask_region1'] = maskwave 
            elif option == 'm2': 
                print "Enter wavelength window to mask out:  blue, red:" 
                maskstr = raw_input(">") 
                try:
                    maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
                except ValueError:
                    print 'Invalid entry. Enter wavelengths separated by commas'
                else:
                    config_pars['mask_region2'] = maskwave
            elif option == 'm3': 
                print "Enter wavelength window to mask out:  blue, red (Angstroms:" 
                maskstr = raw_input(">") 
                try:
                    maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
                except ValueError:
                    print 'Invalid entry. Enter wavelengths separated by commas'
                else:
                    config_pars['mask_region3'] = maskwave
            elif option == 't': 
                print "The current transition wavelength is: " + str(config_pars['transition_wave']) + "\nEnter the wavelength for the G102 to G141 transition:" 
                try:
                    config_pars['transition_wave'] = float(raw_input(">"))
                except ValueError:
                    print 'Invalid entry. Enter wavelength of grism transition.'
            elif option == 'nodes': 
                strnw = ','.join(str(nw) for nw in config_pars['node_wave'])
                print "Enter Wavelengths for Continuum Spline: w1, w2, w3, w4, ...." 
                print "(current node wavelengths are: %s)" %strnw #+ str(config_pars['node_wave'] )
                nodestr = raw_input(">") 
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
            elif option == 'reset': 
                print "Reset configuration parameters, fwhm guess, and zguess to default values" 
                config_pars = read_config('default.config') 
                fwhm_guess = 2.35 *  a_image_obj * config_pars['dispersion_red'] ### fwhm is defined for the red side, regardless of where line is.
                zguess = lamline /6564. - 1 
            elif option == 'bluecut' :
                print "The current blue cutoff is: " + str(config_pars['lambda_min']) +"\nChange the blue cutoff of G102:"
                try:
                    config_pars['lambda_min'] = float(raw_input(">")) 
                except ValueError:
                    print 'Invalid entry. Enter wavelength of blue cutoff.'
            elif option == 'redcut': 
                print "The current red cutoff is: " + str(config_pars['lambda_max']) + "\nChage the red cutoff of G141:" 
                try:
                    config_pars['lambda_max'] = float(raw_input(">")) 
                except ValueError:
                    print 'Invalid entry. Enter wavelength of red cutoff.'

            elif option =='n':  
                nlines_found_cwt = np.size(lamlines_found)  
                index_of_strongest_line = index_of_strongest_line +1 
               
                if index_of_strongest_line < (nlines_found_cwt): 
                    lamline = lamlines_found[index_of_strongest_line] 
                    zguess = lamline /6564. - 1

                else:
                    print 'There are no other automatically identified peaks. Select another option.'
                    ### stay at current line
                    index_of_strongest_line -= 1

                
            #### other lines      
            elif option=='ha':
                zguess=(lamline/lam_Halpha)-1
            elif option=='hb':
                zguess=(lamline/lam_Hbeta)-1
            elif option=='o2':
                zguess=(lamline/lam_Oii)-1
            elif option=='o31':
                zguess=(lamline/lam_Oiii_1)-1
            elif option=='o32':
                zguess=(lamline/lam_Oiii_2)-1
            elif option=='s2':
                zguess=(lamline/lam_Sii)-1
            elif option=='s31':
                zguess=(lamline/lam_Siii_1)-1
            elif option=='s32':
                zguess=(lamline/lam_Siii_2)-1
            #elif option=='la':
            #    zguess=(lamline/lam_Lya)-1
            elif option=='b':
                print "Going back."
                ### this is needed so going back doesn't write everything to the file
                zset=0
                i=i-2  ### below, it adds 1 to i before moving on.  so subtract 2 to go back one.
                next=1
            elif option=='q':
                print "Quitting; saved through previously completed object."
                return 0
            


            ### redo display things
            elif option=='lin':
                if len(g102zerox)>0:
                    show2dNEW('G102',parnos[i],int(objid_unique[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'linear')
                if len(g141zerox)>0:
                    show2dNEW('G141',parnos[i],int(objid_unique[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'linear')
            elif option=='log':
                if len(g102zerox)>0:
                    show2dNEW('G102',parnos[i],int(objid_unique[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'log')
                if len(g141zerox)>0:
                    show2dNEW('G141',parnos[i],int(objid_unique[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'log')
            elif len(option)>6 and option[0:6]=='zs102=':
                vals=option[6:]
                zran=vals.split(',')
                if len(zran)!=2 or isFloat(zran[0])==False or isFloat(zran[1])==False:
                    print "Invalid zrange."
                elif len(g102zerox)>0:
                    show2dNEW('G102',parnos[i],int(objid_unique[i]),g102firstx,g102firsty,g102firstlen,g102firstwid,g102firstid,g102zerox,g102zeroy,g102zeroid,'linear',zran1=float(zran[0]),zran2=float(zran[1]))
            elif len(option)>6 and option[0:6]=='zs141=':
                vals=option[6:]
                zran=vals.split(',')
                if len(zran)!=2 or isFloat(zran[0])==False or isFloat(zran[1])==False:
                    print "Invalid zrange."
                elif len(g141zerox)>0:
                    show2dNEW('G141',parnos[i],int(objid_unique[i]),g141firstx,g141firsty,g141firstlen,g141firstwid,g141firstid,g141zerox,g141zeroy,g141zeroid,'linear',zran1=float(zran[0]),zran2=float(zran[1]))
            elif option=='dc':
                showDirectNEW(objid_unique[i],i)
                if show_dispersed:  # MB
                    showDispersed(objid_unique[i],i)
            elif option=='dr':
                reloadReg()
            else:
                print "Invalid entry.  Try again."
            print "OK"
        
        
        if zset == 1:
            # I don't even remember what flagcont is. it is a remnant of the old code. 
            WriteToCatalog(linelistoutfile, parnos[0], objid_unique[i],ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont)
            
            write_fitdata(fitdatafilename, spec_lam, spec_val, spec_unc, spec_con, spec_zer, full_fitmodel, full_contmodel, mask_flg)

            fitspec_pickle = open(fitdatafilename + '.pickle', 'wb') 
            output_meta_data = [parnos[0], objid_unique[i], ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont, config_pars]
            pickle.dump(output_meta_data, fitspec_pickle) 
            fitspec_pickle.close() 


        ### write comments always 
        WriteComments(commentsfile, parnos[0], objid_unique[i], comment)    ### if we go back to the previous objects, duplicate comments will still be written 
       
        i=i+1
        if i>=len(objid_unique):
            raw_input("You have reached the last object.\nThis is the end of the field.\nPress any key to quit...")

   # printLLout(linelistoutfile,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment)
    #printLCout(linelistoutfile_cont,parnos,grism,objid,wavelen,npix,ston,flag,flagcont,setzs,comment)
   
   
   # Clean up temp files
   # if save_temp:
   #     if os.path.exists(lltemp)==1:
   #         os.unlink(lltemp)
   #     if os.path.exists(lctemp)==1:
   #         os.unlink(lctemp)
    if os.path.exists('./tempcoo.dat')==1:
        os.unlink('./tempcoo.dat')
    if os.path.exists('./temp_zero_coords.coo')==1:
        os.unlink('./temp_zero_coords.coo')
    if os.path.exists('./temp110.fits')==1:
        os.unlink('./temp110.fits')
    if os.path.exists('./temp140.fits')==1:
        os.unlink('./temp140.fits')
    if os.path.exists('./temp160.fits')==1:
        os.unlink('./temp160.fits')
    if os.path.exists('./temp_zero_coords.reg')==1:
        os.unlink('./temp_zero_coords.reg')
    if os.path.exists('G102_trace.reg') ==True:
        os.unlink('G102_trace.reg') 
    if os.path.exists('G141_trace.reg') == True: 
        os.unlink('G141_trace.reg')
    
def WriteToCatalog(catalogname, parnos, objid, ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont): #parnos, objid are scalar not array. 
    if os.path.exists(catalogname) == False: 
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
        cat.write('#19 hg_flux \n') 
        cat.write('#20 hg_err \n')  
        cat.write('#21 hg_EW_obs \n') 
        cat.write('#22 hb_flux \n') 
        cat.write('#23 hb_err \n')  
        cat.write('#24 hb_EW_obs \n') 
        cat.write('#25 oiii_flux [both lines] \n') 
        cat.write('#26 oiii_err [both lines] \n')  
        cat.write('#27 oiii_EW_obs [both lines] \n') 
        cat.write('#28 hanii_flux \n') 
        cat.write('#29 hanii_err \n') 
        cat.write('#30 hanii_EW_obs \n')  
        cat.write('#31 sii_flux \n') 
        cat.write('#32 sii_err \n') 
        cat.write('#33 sii_EW_obs \n')  
        cat.write('#34 siii_9069_flux \n') 
        cat.write('#35 siii_9069_err \n') 
        cat.write('#36 siii_9069_EW_obs \n') 
        cat.write('#37 siii_9532_flux \n') 
        cat.write('#38 siii_9532_err \n')  
        cat.write('#39 siii_9532_EW_obs \n') 
        cat.write('#40 he1_10830_flux \n') 
        cat.write('#41 he1_10830_err \n') 
        cat.write('#42 he1_10830_EW_obs \n') 
        cat.write('#43 ContamFlag \n') 
       # cat.write('#45 Comment \n') 


    else: 
        cat = open(catalogname , 'a') 
    
    outstr= '{:<8d}'.format(parnos) + \
            '{:<6d}'.format(objid) +\
            '{:<12.6f}'.format(ra_obj) + \
            '{:<12.6f}'.format(dec_obj) + \
            '{:<8.2f}'.format(jmag_obj) + \
            '{:<8.2f}'.format(hmag_obj) + \
            '{:<8.3f}'.format(a_image_obj) + \
            '{:<8.3f}'.format(b_image_obj) + \
            '{:>8.4f}'.format(fitresults['redshift']) + '{:>10.4f}'.format(fitresults['redshift_err']) +\
            '{:>10.4f}'.format(fitresults['dz_oiii'])  + \
            '{:>10.4f}'.format(fitresults['dz_oii'])   + \
            '{:>10.4f}'.format(fitresults['dz_siii_he1']) +\
            '{:>10.2f}'.format(fitresults['fwhm_g141']) + '{:>10.2f}'.format(fitresults['fwhm_g141_err'])  +  \
            '{:>13.2e}'.format(fitresults['oii_flux'])  + '{:>13.2e}'.format(fitresults['oii_error']) +  '{:>13.2e}'.format(fitresults['oii_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['hg_flux']) + '{:>13.2e}'.format(fitresults['hg_error']) + '{:>13.2e}'.format(fitresults['hg_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['hb_flux']) + '{:>13.2e}'.format(fitresults['hb_error']) + '{:>13.2e}'.format(fitresults['hb_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['oiii_flux']) + '{:>13.2e}'.format(fitresults['oiii_error']) + '{:>13.2e}'.format(fitresults['oiii_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['hanii_flux']) + '{:>13.2e}'.format(fitresults['hanii_error']) + '{:>13.2e}'.format(fitresults['hanii_ew_obs']) + \
            '{:>13.2e}'.format(fitresults['sii_flux']) + '{:>13.2e}'.format(fitresults['sii_error']) + '{:>13.2e}'.format(fitresults['sii_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['siii_9069_flux']) + '{:>13.2e}'.format(fitresults['siii_9069_error']) + '{:>13.2e}'.format(fitresults['siii_9069_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['siii_9532_flux']) + '{:>13.2e}'.format(fitresults['siii_9532_error']) + '{:>13.2e}'.format(fitresults['siii_9532_ew_obs']) +\
            '{:>13.2e}'.format(fitresults['he1_flux']) + '{:>13.2e}'.format(fitresults['he1_error']) +  '{:>12.2e}'.format(fitresults['he1_ew_obs']) +\
             '   ' + '{:<6d}'.format(flagcont) + ' \n' 


    cat.write(outstr) 
    cat.close()



def WriteComments(filename, parnos, objid, comment):
    if os.path.exists(filename) == False: 
        cat = open(filename, 'w') 
    else: 
        cat  = open(filename, 'a') 

    outstr = '{:<8d}'.format(parnos) + \
            '{:<6d}'.format(objid) +\
            comment +  '\n' 

    cat.write(outstr) 
    cat.close() 





def UpdateCatalog(linelistoutfile):

    allDirectoryFiles=os.listdir('./fitdata/') 
    objid_list = [] 
    for obj in allDirectoryFiles: 
        x = obj.split('_')[2]
        objid_list.append(int(x)) 
        Parno = obj.split('_')[0]   # this is inefficient, but I don't care.
    objid_list = np.sort(np.unique(np.array(objid_list)))
    for obj in objid_list:
        print obj
        inpickle = './fitdata/' + Parno + '_BEAM_'  + str(obj) +'_fitspec.pickle' 
        fileObject = open(inpickle,'r') 
        alldata = pickle.load(fileObject)  
        ## definition from above 
       #                      0          1                 2      3        4           5            6         7         8           9         10
       # output_meta_data = [parnos[0], objid_unique[i], ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont, config_pars]
        parnos   = alldata[0] 
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

        WriteToCatalog(linelistoutfile, parnos, objid_unique,ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont)
 






