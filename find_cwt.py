from glob import glob
from wisp_analysis import *
import pdb
#import os

def find_cwt(lam, flux, err, zeros, fwhm_est_pix, beam_name, config_pars, plotflag = False):
    ###### inputs for fiddling:
    
    cont_medfilt = int(config_pars['cont_medfilt']) ### window for median filter
    max_width = config_pars['maxwidth'] * fwhm_est_pix    ### not sure if I want disp_blue here.
    min_width = config_pars['minwidth']
    dw  = (max_width - min_width) / config_pars['nwidths'] 

    #dw_log = (np.log10(max_width) - np.log10(min_width)) / config_pars['nwidths'] 
    #wlog = np.log10(min_width)  + np.arange(config_pars['nwidths']) * dw_log 
    widths = min_width + np.arange(config_pars['nwidths']) * dw
    #print widths 
    #widths = ( np.arange(config_pars['nwidths'])  - config_pars['nwidths']/2.) * config_pars['dw']  +  fwhm_est/config_pars['dispersion']  ### kernel widths for use with the ricker filter in the cwt. 
    max_distance_ridge = widths * config_pars['max_dist_ridge_scl'] + config_pars['max_dist_ridge_const']    ### if peak in a row of the cwt matrix is off by more than this many pixels, new/separate ridge
    gap_allowed_between_ridges = config_pars['gap_btw_ridges']  #### gap between ridges can be no more than N pixels, otherwise new/separate ridge/peak
    snr_cwt = config_pars['snr_cwt']  ### snr for the cwt line finder
    noise_cut_cwt = config_pars['noise_cut_cwt'] ### noise cut for cwt to estimate noise 
    min_length_cwt = config_pars['min_length_ridge']  ### minimum length of a cwt ridge to be considered real
    edge_reject = config_pars['edge_reject']      #### reject cwt detections within 5 pixcels of edge
    sn_thresh_cont_check = config_pars['n_sigma_above_cont'] ### in step2, require cwt line candidates to have npix_thresh abvove sn_threh_cont_check
    npix_thresh = config_pars['npix_thresh']  
    min_line_contrast = config_pars['min_line_contrast']  ### minimum allowed ; rejects low EW lines. 

   
        

     ###plotting ####
    if plotflag == True: 
        f, axarr = plt.subplots(2, 1, figsize=(8, 8)) 
        #plt.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.13, hspace=0.2, wspace=0.2)
        w=np.where((lam > config_pars['lambda_min']) & (lam < config_pars['lambda_max'])) 
        spec_max  = np.max(flux[w]) 
        axarr[1].plot(lam,flux, ls='steps-mid', color = 'k') 
        axarr[1].axis([np.min(lam), np.max(lam), -0.5e-18, 1.3 * spec_max]) 




    ### run the order filter
    
    cont_filter = si.medfilt(flux, cont_medfilt) 
    
    if plotflag ==True :
        ###  continuum model
        axarr[1].plot(lam, cont_filter+err*sn_thresh_cont_check, color= 'orange') 
        axarr[1].plot(lam, cont_filter) 

   
    ### calculate and show contiuous wavelet transform array: 
    #w=np.isfinite(cwarray)
    cwarray = si.cwt(flux, si.ricker, widths) 
    if plotflag ==True:
        axarr[0].imshow(cwarray, vmax = .6 * spec_max, vmin  = -0.2 * spec_max, aspect = 'auto') 
 

    ###### find the peaks and overplot them. 
    peaks= si.find_peaks_cwt(flux, widths, wavelet = si.ricker, max_distances=max_distance_ridge, gap_thresh=gap_allowed_between_ridges, 
               min_snr=snr_cwt, noise_perc = noise_cut_cwt, min_length=min_length_cwt )
    if plotflag == True: 
        axarr[1].plot(lam[peaks], flux[peaks], 'ro', ms=7)
    peaks = np.array(peaks)


    ### reject peaks near the edge 
    w= np.where( (peaks > edge_reject) & (peaks < np.size(lam) - edge_reject)) 
    peaks = peaks[w[0]]
    
    ### reject lines with presumably low EWs 
    try:
        peak_contrast = (flux[peaks] - cont_filter[peaks])/cont_filter[peaks]
    except IndexError:
        ### added this step in case peaks fail the above conditional
        peaks = []
    else:
        w=np.where(peak_contrast > min_line_contrast)
        peaks  = peaks[w[0]]

    ### reject peaks that are actually zero orders 
    if np.size(peaks)  > 0: 
        zero_check = [] 
        mef102 = mef = '../G102_DRIZZLE/aXeWFC3_G102_mef_ID'+beam_name + '.fits'
        mef141 = mef = '../G141_DRIZZLE/aXeWFC3_G141_mef_ID'+beam_name + '.fits' 
        hdu102 = fits.open(mef102)
        hdu141 = fits.open(mef141)

        for pk in peaks: 
            
            ### check 1, get the x-pix where we don't zero order check 
            if lam[pk] > config_pars['transition_wave']: 
                xthresh = 820.  
                band = 'G141'
                header = hdu141[1].header
            else: 
                xthresh = 760. 
                band = 'G102'
                header = hdu102[1].header 
     
            xpix_line_fullframe = (lam[pk] - header['CRVAL1']) / header['CDELT1'] + header['CRPIX1']  + header['BB0X']


            ### check 2: if the zeros are flagged in the extract
            zero_clip = zeros[pk-5:pk+5]
            w=np.where(zero_clip == 3.) 

            if ( (len(w[0]) > 0) | (xpix_line_fullframe > xthresh)):
                #print 'Zero order detected and rejected:', lam[pk], xpix_line_fullframe
                zero_check.append(False)   ## fails zero check, so False will remove 
            else: 
                zero_check.append(True) 
        zero_check = np.array(zero_check) 
        peaks = peaks[zero_check] 
        hdu102.close()
        hdu141.close()



     
    if np.size(peaks) > 0:
    
        ### count contiguous pixels above the noise threshold: 

        snr_thresh = sn_thresh_cont_check 
        npix_peak = [] 
        line_snr_guess = []
        for i in peaks:
            ### first of all, is peak above threshold:
            if flux[i] > cont_filter[i] + err[i]  * snr_thresh : 
                pixel_count = 1
            
                cond = 0 
                j = i + 1

                while ((cond == 0.) & (j < np.size(flux) -1)) :
                    if flux[j] > cont_filter[j] + err[j]  * snr_thresh: 
                        pixel_count = pixel_count + 1
                        #cond = 0 
                        j = j + 1 
                    else:  
                        cond = 1.
      
           
                cond = 0
                j = i-1  
                while ((cond == 0) & (j > 0)):
                    if flux[j] > cont_filter[j] + err[j]  * snr_thresh: 
                         pixel_count = pixel_count + 1
                         #cond = 0 
                         j = j-1
                    else:  
                        cond = 1.

            else: 
                pixel_count =0

            npix_peak.append(pixel_count)
            #### crudely estimate the snr of each line candidate
            
            if lam[i] > config_pars['transition_wave']: 
                disp_est = config_pars['dispersion_red'] 
            else: 
                disp_est = config_pars['dispersion_blue']

            fwhm_est = fwhm_est_pix * disp_est
           
            w=np.where( (lam > lam[i] - 0.5 * fwhm_est) & (lam < lam[i] + 0.5 * fwhm_est)) 
            #w=range(i-int(0.5 *fwhm_est_pix), i+ int(0.5 * fwhm_est_pix)) 

    
            line_signal_guess = np.sum( (flux[w] - cont_filter[w]) * disp_est)
            line_noise_guess =  np.sqrt(np.sum((err[w]* disp_est)**2 ) ) 
            line_snr_guess.append(line_signal_guess/line_noise_guess)


        npix_peak = np.array(npix_peak) 
        line_snr_guess = np.array(line_snr_guess) 
        w=np.where(npix_peak >= npix_thresh)# & (line_snr_guess > 1))  
        real_peaks = peaks[w[0]] 
        npix_real= npix_peak[w[0]] 
        snr_real = line_snr_guess[w[0]]
    else:
         real_peaks = [] 
         npix_real = []
         snr_real = []
         peaks = [] 


    if plotflag==True:
            axarr[1].plot(lam[real_peaks], flux[real_peaks], 'rs', ms=9, markeredgecolor= 'r', markerfacecolor = 'none', markeredgewidth=2)
            plt.title(beam_name) 
            #figname = 'temp_'  + beam_name + '.pdf'
            #print 'saving fig for: '  + figname  
            #plt.savefig(figname)
            #plt.close()
            plt.show(block=True)
    return [lam[real_peaks], flux[real_peaks], npix_real, snr_real, cwarray, cont_filter, lam[peaks], flux[peaks]] 





def loop_field_cwt(): 
    #### no inputs.  run from the  inside the data directory.
    if os.path.exists('linelist') == False: 
        os.mkdir('linelist') 

    g102files = glob('*G102_BEAM_*A.dat')
    g102files.sort()
    g141files = glob('*G141_BEAM_*A.dat')
    g141files.sort()

###    os.system('ls *G102_BEAM_*A.dat > linelist/g102_spec.list')
###    os.system('ls *G141_BEAM_*A.dat > linelist/g141_spec.list')
    
    config_pars = read_config('default.config')

###    g102list = asciitable.read('linelist/g102_spec.list', format = 'no_header') 
###    g102files = g102list['col1']
###    g141list = asciitable.read('linelist/g141_spec.list', format = 'no_header')
###    g141files = g141list['col1']
    ### the sizes of the sources are used as a rough estimate 
    ### these are the same in all filters in pixel coordinates
    ### so remove the necessity for both filters
    
    # get list of available catalogs
    catalogs = glob('../DATA/DIRECT_GRISM/fin_F*cat')
    catalogs.sort()
    # the "blue" catalog will be the first filter
    cat = asciitable.read(catalogs[0])
    if len(catalogs) > 1:
        cat2 = asciitable.read(catalogs[1])
    # a_image is the same in both catalogs
    a_images = cat['col5']
    beam_se = cat['col2']

###    blue_se = asciitable.read('../DATA/DIRECT_GRISM/fin_F110.cat') 
###    red_se  = asciitable.read('../DATA/DIRECT_GRISM/fin_F160.cat') 
###
###    a_image_blue = blue_se['col5'] 
###    a_image_red = red_se['col5']
###    beam_se = blue_se['col2']    #### doesn't matter if it comes from red/blue because these cats are matched. 
    
    outfile = open('linelist/temp', 'w') 
    config_pars['transition_wave'] = 11700.
    for filename in g102files:
        #filename = 'Spectra/Par302_G102_BEAM_1A.dat'
        ### get spectral data        
        spdata = asciitable.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero']) 
        trimmed_spec = trim_spec(spdata, None, config_pars) 
        
        ### look up the object in the se catalog and grab the a_image.
        beam = float(filename.split('_')[-1].split('A')[0]) 
        parno = filename.split('Par')[1].split('_')[0]

        w=np.where(beam_se == beam) 
        w=w[0]    # because of the stupid tuple thing
###        a_image = a_image_blue[w][0]
        a_image = a_images[w][0]
        fwhm_est_pix = a_image * 2

        
        #### unpack spectrum and check that it is long enough to proceed. 
        lam = trimmed_spec[0] 
        flux_corr = trimmed_spec[1] - trimmed_spec[3] 
        err = trimmed_spec[2]
        zeros = trimmed_spec[4]
        if len(lam) < config_pars['min_spec_length']: 
            continue 
        
        ### cwt it, unpack and write results. 
        #print config_pars['npix_thresh'] , config_pars['n_sigma_above_cont'], config_pars['transition_wave']
        g102_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
        lam_cwt = g102_cwt[0] 
        flam_cwt = g102_cwt[1] 
        npix_cwt = g102_cwt[2] 
        snr_cwt = g102_cwt[3]
        for i in np.arange(len(lam_cwt)):
            print beam, 'G102', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i]
            outfile.write(parno + '  G102  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) +  ' '  +  str(npix_cwt[i])   +  '   ' + str(snr_cwt[i]) + '\n') 
        if config_pars['n_sigma_for_2pix_lines'] != False: 

            config_pars['npix_thresh'] = 2 
            config_pars['n_sigma_above_cont'] = config_pars['n_sigma_for_2pix_lines']
           # print config_pars['npix_thresh'] , config_pars['n_sigma_above_cont'], config_pars['transition_wave']
            g102_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
            lam_cwt = g102_cwt[0] 
            flam_cwt = g102_cwt[1] 
            npix_cwt = g102_cwt[2] 
            snr_cwt = g102_cwt[3]
            for i in np.arange(len(lam_cwt)):
                print beam, 'G102', lam_cwt[i], npix_cwt[i] , fwhm_est_pix, snr_cwt[i]
                outfile.write(parno + ' G102   '  + str(int(beam)) + ' ' + str(lam_cwt[i]) +  ' '  +  str(npix_cwt[i])  +  '   ' + str(snr_cwt[i]) + '\n') 
        
        ### go back to the beginning with the old config pars 
        config_pars = read_config('default.config')
        config_pars['transition_wave'] = 11700.


    config_pars['transition_wave'] = 11100.
    for filename in g141files:
        #filename = 'Spectra/Par302_G141_BEAM_1A.dat'
        spdata = asciitable.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero']) 
        trimmed_spec = trim_spec(None, spdata, config_pars) 
        beam = float(filename.split('_')[-1].split('A')[0]) 
        parno = filename.split('Par')[1].split('_')[0]
        w=np.where(beam_se == beam) 
        w=w[0]    # because of the stupid tuple thing
###        a_image = a_image_red[w][0]  
        a_image = a_images[w][0]  
        lam = trimmed_spec[0] 
        flux_corr = trimmed_spec[1] - trimmed_spec[3] 
        err = trimmed_spec[2] 
        zeros = trimmed_spec[4] 
        if len(lam) < config_pars['min_spec_length']: 
            continue
        fwhm_est_pix = a_image * 2 
        #print config_pars['npix_thresh'] , config_pars['n_sigma_above_cont'], config_pars['transition_wave']
        config_pars
        g141_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix,str(int(beam)), config_pars, plotflag=False)
        lam_cwt = g141_cwt[0] 
        flam_cwt = g141_cwt[1] 
        npix_cwt = g141_cwt[2] 
        snr_cwt = g141_cwt[3]
        for i in np.arange(len(lam_cwt)):
            print beam, 'G141', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i]
            outfile.write(parno +  '  G141  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) +  ' '  +  str(npix_cwt[i]) +   '   ' + str(snr_cwt[i]) +'\n') 
        if config_pars['n_sigma_for_2pix_lines'] != False: 

            config_pars['npix_thresh'] = 2 
            config_pars['n_sigma_above_cont'] = config_pars['n_sigma_for_2pix_lines'] 
            #print config_pars['npix_thresh'] , config_pars['n_sigma_above_cont'], config_pars['transition_wave']
            g141_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
            lam_cwt = g141_cwt[0] 
            flam_cwt = g141_cwt[1] 
            npix_cwt = g141_cwt[2] 
            snr_cwt = g141_cwt[3] 
            for i in np.arange(len(lam_cwt)):
                print beam, 'G141', lam_cwt[i], npix_cwt[i], snr_cwt[i]
                outfile.write(parno + '  G141  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) +  ' '  +  str(npix_cwt[i])  + '   ' + str(snr_cwt[i]) + '\n') 
        ### go back to the beginning with the old config pars 
        config_pars = read_config('default.config')
        config_pars['transition_wave'] = 11200.

    outfile.close()
    tab = asciitable.read('linelist/temp', format = 'no_header')
    par = tab['col1']
    grism = tab['col2']
    beam = tab['col3']  
    wave = tab['col4'] 
    npix = tab['col5']
    snr = tab['col6']
    s=np.argsort(beam) 
    beam = beam[s] 
    grism = grism[s] 
    wave =wave[s]
    npix = npix[s]
    snr = snr[s]
    par = par[0] 
    beams_unique = np.unique(beam)  
    outfile = open('linelist/Par'+str(par) + 'lines.dat', 'w') 
    for b in beams_unique:
        #### do the g102 for b
        w = (beam == b) & (grism== 'G102')
        
        waves_obj = wave[w] 
        npix_obj = npix[w] 
        snr_obj = snr[w] 
        
        waves_uniq, ind = np.unique(waves_obj, return_index = True) 
        npix_uniq = npix_obj[ind] 
        snr_uniq = snr_obj[ind]
        s = np.argsort(waves_uniq) 
        waves_final_g102 = waves_uniq[s] 
        npix_final_g102 = npix_uniq[s] 
        snr_final_g102 = snr_uniq[s]
       
        for lam, npx, sn in zip(waves_final_g102, npix_final_g102, snr_final_g102): 
            outfile.write(str(par) +  ' G102  ' + str(b) +  '  '   + str(lam) +  '   ' +str(npx)  + '  ' + str(sn)  +  '\n')
        

        ### do the g141 for b
        w = (beam == b) & (grism== 'G141')
        waves_obj = wave[w] 
        npix_obj = npix[w]
        snr_obj = snr[w] 
        
        waves_uniq, ind = np.unique(waves_obj, return_index = True) 
        npix_uniq = npix_obj[ind] 
        snr_uniq = snr_obj[ind] 
        s = np.argsort(waves_uniq) 
        waves_final_g141 = waves_uniq[s] 
        npix_final_g141 = npix_uniq[s] 
        snr_final_g141 = snr_uniq[s] 

        #wave_obj_g141 =  np.sort(np.unique(wave[w]))
        for lam, npx, sn in zip(waves_final_g141, npix_final_g141, snr_final_g141): 
            outfile.write(str(par) + '  G141  ' + str(b) +  ' ' + str(lam) + '  ' + str(npx)  +'  ' +  str(sn) +  '  \n')
    outfile.close()            






def test_obj_cwt(parno, beamno, configfile): 
    blue_se = asciitable.read('../DATA/DIRECT_GRISM/fin_F110.cat') 
    red_se  = asciitable.read('../DATA/DIRECT_GRISM/fin_F160.cat')
    a_image_blue = blue_se['col5'] 
    a_image_red = red_se['col5']
    beam_se = blue_se['col2'] 
    config_pars = read_config(configfile)
    bluefile = 'Par'+str(parno) + '_G102_BEAM_'+str(beamno)+'A.dat'
    redfile =  'Par'+str(parno) + '_G141_BEAM_'+str(beamno)+'A.dat' 
    
    spdata_blue = asciitable.read(bluefile, names = ['lambda', 'flux', 'ferror', 'contam', 'zero']) 
    trimmed_spec_blue= trim_spec(spdata_blue, None, config_pars)  
    
    ### do the blue side 
    lam = trimmed_spec_blue[0] 
    flux_corr = trimmed_spec_blue[1] - trimmed_spec_blue[3] 
    err = trimmed_spec_blue[2]
    zeros = trimmed_spec_blue[4] 
    config_pars['transition_wave']  = 11700.

    if len(lam) < config_pars['min_spec_length']: 
        print 'Short spec. skip it!'
    else: 
        w=np.where(beam_se == beamno) 
        w=w[0]    # because of the stupid tuple thing
        a_image = a_image_blue[w][0] 
        fwhm_est_pix = a_image * 2
        g102_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(beamno), config_pars, plotflag=True)  
        print g102_cwt[0], g102_cwt[1], g102_cwt[2], fwhm_est_pix


    ### do the red side 
    config_pars['transition_wave'] = 11200.
    spdata_red = asciitable.read(redfile, names = ['lambda', 'flux', 'ferror', 'contam', 'zero']) 
    trimmed_spec_red= trim_spec(None, spdata_red, config_pars)
    lam = trimmed_spec_red[0] 
    flux_corr = trimmed_spec_red[1] - trimmed_spec_red[3] 
    err = trimmed_spec_red[2]
    zeros = trimmed_spec_red[4] 
    if len(lam) < config_pars['min_spec_length']: 
        print 'Short spec. skip it!'
    else: 
        w=np.where(beam_se == beamno) 
        w=w[0]    # because of the stupid tuple thing
        a_image = a_image_red[w][0] 
        fwhm_est_pix = a_image * 2
        g141_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(beamno), config_pars, plotflag=True)  
        print g141_cwt[0], g141_cwt[1], g141_cwt[2], fwhm_est_pix 

    #### I do not understand why it seems necessary to repeat this, but repeat calls to the function seem to need it.
    #config_pars['transition_wave'] = 11700. 










        

 

    
