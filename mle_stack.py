from wisp_analysis import * 
import sys

#### maximum likelihood calculation of dust extinction from multiple lines. 

def k_lambda(wave): 

    ##### CALZETTI+2000 
    x = 10000.0/wave   ### wave in Angsroms 
    w1 = np.where((wave >= 6300) & (wave <= 22000))
    w2 = np.where((wave >= 912) & (wave <=  6300))
    klam = wave * 0 
    klam[w1] = 2.659*(-1.857 + 1.040*x[w1]) + 4.05 
    klam[w2] = 2.659*(-2.156 + 1.509*x[w2] -0.198 *x[w2]**2 + 0.011*x[w2]**3) + 4.05 
    return klam 


def calc_levels(arr):
    sorted_likelihoods = np.sort(arr, axis=None)  
    normed_likelihoods = sorted_likelihoods/np.sum(sorted_likelihoods)
    reversed_likelihoods = sorted_likelihoods[:: -1] 
    reversed_normed_likelihoods = normed_likelihoods[::-1] 
    cumulative_likelihoods  = [] 
    for i in np.arange(np.size(sorted_likelihoods)): 
        cumulative_likelihoods.append(np.sum(reversed_normed_likelihoods[:i+1])) 
    cumulative_likelihoods = np.array(cumulative_likelihoods)  
    w=np.where(cumulative_likelihoods > 0.68) 
    lev1= np.max(reversed_likelihoods[w])
    w=np.where(cumulative_likelihoods > 0.95) 
    lev2 = np.max(reversed_likelihoods[w])
    w=np.where(cumulative_likelihoods > 0.997) 
    lev3 = np.max(reversed_likelihoods[w]) 

    return [lev3, lev2, lev1] 



def calc_likelihood(input_meas, likelihood_func, curti=True, strom=False): 
    
    if (curti == True) & (strom == True) : 
        print "You can't use both Curti and Strom Calibrations. Pick one and try again!" 
        sys.exit(1) 

    if (curti == False) & (strom == False): 
        print "Select a metallicity calibration" 
        sys.exit()

    #### input meas comes from measure_stack 
    ### likelihood func is a fits file

    
    ####   intrinsic balmer decrement for 10,000K 
    hahb_int = 2.86 
    hghb_int = 0.468  
    hdhb_int = 0.259 
     
    
    meas = asciitable.read(input_meas, format = 'fixed_width') 
    meas.colnames
    fluxes = meas['flux_norm'] 
    errors = meas['flux_norm_err'] 
    #cont = meas['median_continuum'] 
    #cont_err=  meas['err_continuum'] 
    ew_rest = meas['ew_rest'] 
    fluxes_cgs = meas['flux_scale'] 
    ew_err =  meas['ew_err'] 




    
    haniiflux = fluxes[0] 
    haniierr =  errors[0] 
    haniiew = ew_rest[0] 
    haew_err = ew_err[0] 
    haniiflux_cgs = fluxes_cgs[0] 
    #hanii_cont_err = cont_err[0]  
    #hanii_cont = cont[0]  

    
    hbflux = fluxes[1] 
    hberr = errors[1]
    hbew = ew_rest[1] 
    hbew_err = ew_err[1] 
    hbflux_cgs = fluxes_cgs[1] 
    #hb_cont_err=  cont_err[1] 
    #hb_cont = cont[1] 

       
    hgoiiiflux = fluxes[2] 
    hgoiiierr = errors[2] 
    hgoiiiew = ew_rest[2] 
    hgoiiiew_err = ew_err[2] 
    hgflux_cgs = fluxes_cgs[2] 
    #hg_cont_err=  cont_err[2] 
    #hg_cont = cont[2] 

    hdflux = fluxes[3] 
    hderr = errors[3]  
    hdew = ew_rest[3] 
    hdew_err = ew_err[3]
    hdflux_cgs = fluxes_cgs[3] 
    #hd_cont_err=  cont_err[3] 
    #hd_cont = cont[3] 


    oiiiflux = fluxes[4]     
    oiiierr =  errors[4] 

    oiiflux = fluxes[5] 
    oiierr = errors[5] 

    siiflux = fluxes[6] 
    siierr = errors[6]

    #he1flux = fluxes[7] 
    #he1err = errors[7] 

    #oiflux = fluxes[8] 
    #oierr = errors[8] 

    #neiiiflux = fluxes[9] 
    #neiiierr = errors[9]  

   
    ### calculate observed line ratios :
    if haniiflux > 0 : 
        haniihb = haniiflux / hbflux  
        haniihb_err = haniihb * np.sqrt( (haniierr/haniiflux)**2 + (hberr/hbflux)**2)  
    else :
        haniihb = 0
        haniihb_err = 1e12
        haew_err = 1e12 


    hgoiiihb = hgoiiiflux / hbflux 
    hgoiiihb_err = hgoiiihb * np.sqrt( (hgoiiierr/hgoiiiflux)**2 + (hberr/hbflux)**2)  
    
    if hdflux > 0: 
        hdhb = hdflux/hbflux 
        hdhb_err  = hdhb * np.sqrt((hderr/ hdflux)**2 + (hberr/hbflux)**2) 

    o32 = oiiiflux / oiiflux ## dust extincted version, we will forward model  in MLE
    o32_err = o32 * np.sqrt( (oiiierr/oiiiflux)**2 + (oiierr/oiiflux)**2) 

    r23_num = oiiiflux + oiiflux  ### this is the dustextincted version-- what we predict in the MLE
    r23_num_err = np.sqrt(oiiierr**2 + oiierr**2)   
    r23 = r23_num/hbflux 
    r23_err = r23 * np.sqrt( (r23_num_err/r23_num)**2 + (hberr/hbflux)**2) 

    r2 = oiiflux/hbflux   ### dust extincted 
    r2_err = r2 * np.sqrt( (oiierr/oiiflux)**2 + (hberr/hbflux)**2) 

    r3 = oiiiflux/hbflux  ### dust extincted
    r3_err = r3 * np.sqrt( (oiiierr/oiiiflux)**2 + (hberr/hbflux)**2) 



    #### calculate the range of allowed OIII 4363 ratios, relative to Hg. 

    C1 = -1 * 1/20.  * oiiiflux/hgoiiiflux   ###  1/40  = 4363/(4959 + 5007) for hot temps, from Osterbrock
    C2 = -1 * 1/500. * oiiiflux/hgoiiiflux    #### 1/500 = 4363/(4959 + 5007) for cool temps, from Osterbrock 

    ### because the above numbers involve scaling hg+oiii 4363 fluxes, the solution for OIII 4363/Hg is quadratic 
    oiii4363hg_max = (-1 + np.sqrt(1 - 4 * C1) )/2. 
    oiii4363hg_min = (-1 + np.sqrt(1 - 4  *C2))/2.   
   
     
    #print 'Observed ratios:  HaNII/Hb, HgOIII/Hb, Hb/Hd' 
    #print haniihb, hgoiiihb, hdhb 


    delta_ebv = 0.02 
    delta_oh = 0.02 
    delta_oiii = 0.02 
    delta_hbabs = 0.3 
    
    
    #if haniiflux > 0: 
    #    niiha_mod = np.arange(0, 0.3, delta_nii)  #### ratio ### try a guess for now.  eventually iterate over niiha. 
    #else : 
       # niiha_mod = np.arange(0, 0.1, 0.05)  #### set this to something small to speed up but not break code.  
    ebv_mod = np.arange(0, 0.5, delta_ebv) 
    oh_mod =  np.arange(7, 10, delta_oh) 
    oiiihg_mod = np.arange(oiii4363hg_min, oiii4363hg_max, delta_oiii)    #### ratio ### try a guess for now, eventually iterate. 
    hb_abs_mod = np.arange(0, 5, delta_hbabs)  #### rest-frame stellar absorption ew, to be compared to the rest frame emission ews. 
    

    k_array = k_lambda(np.array([6564., 4861., 4341., 4102., 5007, 3727]))
    k_ha = k_array[0]
    k_hb=  k_array[1] 
    k_hg = k_array[2]  
    k_hd = k_array[3]
    k_oiii = k_array[4] 
    k_oii = k_array[5] 

    x1= np.size(ebv_mod) 
    x2 = np.size(oh_mod) 
    x3= np.size(oiiihg_mod) 
    x4 = np.size(hb_abs_mod) 

    likelihood = np.zeros( (x1, x2, x3, x4))

    for i in np.arange(x1):
        for j in np.arange(x2): 
            for k in np.arange(x3):
                for m in np.arange(x4): 

                    #### calculate model line ratios according to curti 

                    if curti == True : 
                        xoh = oh_mod[j] - 8.69


                        r23_mod = 10** (0.527 - 1.569 * xoh - 1.652*xoh**2 -0.421 * xoh**3 )    #using r2 negates this, I think.
                        r2_mod =  10**(0.418 - 0.961 * xoh - 3.505 *xoh**2 - 1.949 * xoh**3)
                        r3_mod = 1.3 *  10**(-0.277 - 3.549 * xoh - 3.593 * xoh**2 - 0.981 * xoh**3)
                        n2_mod = 1.3 *  10**(-0.489  + 1.513 * xoh - 2.554 * xoh**2 - 5.293 * xoh**3 - 2.867 * xoh**4)
                        o32_mod = 1.3  * 10**(-0.691 + - 2.944 * xoh - 1.308 * xoh**2) 

                    if strom == True :
                        xoh = oh_mod[j] - 8.24 
                        r23_mod = (0.85 - xoh**2)/0.87 
                        
                        xoh_n2  = oh_mod[j]  - 8.77
                        n2_mod = 1.3 *  10**(xoh_n2/0.34)


                        ### this is unsavory, but in order to extinct r23 from strom,  I need to know what the OII/Hb ratio is 
                        ### independent or r23.  so I'm taking the value from Curti.  might want to shift metallicty by 0.2 dex to be consistent 
                        ### with difference from strom. need to think about which way. 
                        xoh_r2 = xoh - 8.69
                        r2_mod =  10**(0.418 - 0.961 * xoh_r2 - 3.505 *xoh_r2**2 - 1.949 * xoh_r2**3)



                    #### apply dust extinction to model line ratios. 
                    if haniiflux > 0: 
                        hahb_mod = hahb_int * 10**(ebv_mod[i] / 2.5 * (k_hb - k_ha)) 
                        
                        
                    hghb_mod = hghb_int * 10**(ebv_mod[i] / 2.5 * (k_hb - k_hg))
                    if hdflux > 0: 
                        hdhb_mod = hdhb_int * 10**(ebv_mod[i] / 2.5 * (k_hb - k_hd)) 


                    ## adding extinction to r23 is more complex. 
                    ### assume extinction at Hb and OIII is the same.
                    ### r2 from curti is used for the strom calibration here. 
                    atten_oiihb = 10**(ebv_mod[i] / 2.5 * (k_hb - k_oii)) 
                    r23_mod = r23_mod - r2_mod * (1 - atten_oiihb) 
   
                    if curti == True : 

                        o32_mod= o32_mod * 10**(ebv_mod[i] / 2.5 * (k_oii - k_oiii)) 
                        r2_mod = r2_mod * 10**(ebv_mod[i] / 2.5 * (k_hb - k_oii)) 
                        r3_mod = r3_mod * 10**(ebv_mod[i] / 2.5 * (k_hb - k_oiii)) 
                


                    #### apply stellar correction; include uncertainty on the stellar correction due to uncertainty on the EW. 
                    if haniiflux > 0 : 
                         ha_stellar_corr = (haniiew + hb_abs_mod[m]/ 1.5) / haniiew
                         ha_stellar_corr_err = hb_abs_mod[m]/1.5 * haew_err/haniiew**2  
                         
                         




                    #### calculate the actual correction  
                    hb_stellar_corr = (hbew + hb_abs_mod[m]) / hbew 
                    hg_stellar_corr = (hgoiiiew + hb_abs_mod[m]) /hgoiiiew 
                    if hdflux > 0 : 
                        hd_stellar_corr = (hdew  + hb_abs_mod[m]) / hdew

                    ## calculate the error on the correction for each line
                    hb_stellar_corr_err = hb_abs_mod[m]   * hbew_err /hbew**2 
                    hg_stellar_corr_err = hb_abs_mod[m] * hgoiiiew_err/hgoiiiew**2 
                    if hdflux > 0: 
                        hd_stellar_corr_err = hb_abs_mod[m] * hdew_err/ hdew**2


                    ### calculate stellar abs corrected line ratios 
                    if haniiflux > 0: 
                         hahb_mod = hahb_mod * hb_stellar_corr / ha_stellar_corr   #### each line gets divided by the correction, because we are un-correcting the model to match the observations


                    hghb_mod = hghb_mod * hb_stellar_corr / hg_stellar_corr
                    if hdflux > 0 :
                        hdhb_mod = hdhb_mod * hb_stellar_corr / hd_stellar_corr  
                    if curti == True: 
                        r2_mod = r2_mod * hb_stellar_corr 
                        r3_mod= r3_mod * hb_stellar_corr
                         ### o32 doesn't get a correction. 


                    ## calculate the error on the line ratios simply due to the stellar absorption correction uncertainty.
                    if haniiflux > 0: 
                        hahb_stellar_err = hb_stellar_corr / ha_stellar_corr * np.sqrt( (ha_stellar_corr_err/ha_stellar_corr)**2 + (hb_stellar_corr_err/ hb_stellar_corr)**2) 

                    hghb_stellar_err = hghb_mod* np.sqrt( (hg_stellar_corr_err/hg_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2)
                    if hdflux > 0: 
                        hdhb_stellar_err = hdhb_mod* np.sqrt( (hd_stellar_corr_err/hd_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2)  
                    
                    if curti == True: 
                        r2_mod_err = r2_mod * hb_stellar_corr_err / hb_stellar_corr 
                        r3_mod_err = r3_mod * hb_stellar_corr_err / hb_stellar_corr 
                   

                    #### add correction for stellar absorption uncertainty to the observed error.
                    if haniiflux > 0 : 
                         haniihb_err_tot = np.sqrt(haniihb_err**2 + hahb_stellar_err**2) 
                    hgoiiihb_err_tot = np.sqrt(hgoiiihb_err**2 + hghb_stellar_err**2) 
                    if hdflux > 0: 
                        hdhb_err_tot = np.sqrt(hdhb_err**2 + hdhb_stellar_err**2) 
                    r2_err_tot = np.sqrt(r2_err**2 + r2_mod_err**2) 
                    r3_err_tot = np.sqrt(r3_err**2 + r3_mod_err**2) 
                   
            
                    #### correct the hahb and hghb ratios for nii and oiii 4363 
                    if haniiflux > 0:
                        haniihb_mod = hahb_mod  * (1 + n2_mod) 

                    hgoiiihb_mod = hghb_mod * (1 + oiiihg_mod[k]) 



                    #### calculate likelihoods 
                    if haniiflux > 0:
                        like1 = np.exp(-1 * (haniihb - haniihb_mod)**2 / (2 * haniihb_err_tot**2)) 
                    else: 
                        like1 =  1.0

                    like2 = np.exp(-1 * (hgoiiihb - hgoiiihb_mod)**2 / (2 * hgoiiihb_err_tot**2)) 
                    if hdflux > 0 : 
                        like3 = np.exp(-1 * (hdhb  - hdhb_mod)**2 / (2 * hdhb_err_tot**2)) 
                    else : 
                        like3 = 1.0 

                    if curti == True: 
                        like4 = np.exp(-1  * (r2_mod - r2)**2 / (2 * r2_err_tot**2) )
                        like5 = np.exp(-1 * (r3_mod - r3)**2 /(2 * r3_err_tot**2)) 
                        like6 = np.exp(-1 * (o32_mod - o32)**2 / (2 * o32_err**2))  #### there is no o32_err_tot, because no propagation of stellar absorption uncertainty.
                    else : 
                        like4 = 1.0
                        like5 = 1.0
                        like6 = 1.0 

                    likelihood[i, j, k, m] = like1 * like2 * like3  * like4 * like5 * like6 
    


    likelihood = likelihood / np.sum(likelihood)
    
        
    ##### evaluate the best fitting model: 
    w=np.where(likelihood == np.max(likelihood))

    best_dust =  ebv_mod[w[0]][0]
    best_oh = oh_mod[w[1]][0]
    best_oiiihg = oiiihg_mod[w[2]][0]
    best_hbabs = hb_abs_mod[w[3]][0]

    if haniiflux > 0: 
        hahb_mod = hahb_int * 10**(best_dust / 2.5 * (k_hb - k_ha)) 
    hghb_mod = hghb_int * 10**(best_dust / 2.5 * (k_hb - k_hg)) 
    if hdflux > 0 : 
        hdhb_mod = hdhb_int * 10**(best_dust / 2.5 * (k_hb - k_hd))  

    #print 'models with dust only, a, g, d' 
    #print hahb_mod, hghb_mod, hdhb_mod 
   
    if haniiflux > 0 : 
        ha_stellar_corr = (haniiew + best_hbabs/ 1.5) / haniiew  #### divide model by this number to match data 
    hb_stellar_corr = (hbew + best_hbabs) / hbew 
    hg_stellar_corr = (hgoiiiew + best_hbabs) /hgoiiiew 
    
    if hdflux > 0: 
        hd_stellar_corr = (hdew  + best_hbabs) / hdew  

    #print 'stellar corrections, a, b, g, d' 
    #print ha_stellar_corr, hb_stellar_corr, hg_stellar_corr, hd_stellar_corr
 
    if haniiflux > 0: 
        ha_stellar_corr_err = best_hbabs/1.5 * haew_err/haniiew**2  
    hb_stellar_corr_err = best_hbabs * hbew_err /hbew**2 
    hg_stellar_corr_err = best_hbabs * hgoiiiew_err/hgoiiiew**2 
    if hdflux > 0 : 
        hd_stellar_corr_err = best_hbabs * hdew_err/ hdew**2  

    #print  'uncertainties on stellar corr'  
    #print ha_stellar_corr_err, hb_stellar_corr_err, hg_stellar_corr_err, hd_stellar_corr_err 

    if haniiflux > 0: 
         hahb_mod = hahb_mod * hb_stellar_corr / ha_stellar_corr ### correction is (Ha / ha_corr) / (Hb/Hb_corr)  
    hghb_mod = hghb_mod * hb_stellar_corr / hg_stellar_corr
    if hdflux > 0: 
        hdhb_mod = hdhb_mod * hb_stellar_corr / hd_stellar_corr  
 
    if haniiflux > 0 : 
        if curti == True :
            xoh_best = best_oh - 8.69 
            best_niiha =   1.3 * 10**(-0.489  + 1.513 * xoh_best - 2.554 * xoh_best**2 - 5.293 * xoh_best**3 - 2.867 * xoh_best**4) 
        if strom  == True : 
            xoh_n2_best = best_oh - 8.77
            best_niiha=  1.3 * 10**(xoh_n2_best/0.34)

        haniihb_mod = hahb_mod  * (1 + best_niiha)
    hgoiiihb_mod = hghb_mod * (1 + best_oiiihg)
    
    #print 'best models with dust, stellar absorption, nii, and oiii4363, a, g, d' 
    #print haniihb_mod, hgoiiihb_mod, hdhb_mod 

    if haniiflux > 0: 
        hahb_stellar_err = hb_stellar_corr / ha_stellar_corr * np.sqrt( (ha_stellar_corr_err/ha_stellar_corr)**2 + (hb_stellar_corr_err/ hb_stellar_corr)**2) 
    hghb_stellar_err = hb_stellar_corr/hg_stellar_corr * np.sqrt( (hg_stellar_corr_err/hg_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2) 
    if hdflux > 0 : 
       hdhb_stellar_err = hb_stellar_corr / hd_stellar_corr * np.sqrt( (hd_stellar_corr_err/hd_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2)  
       hdhb_err_tot=  np.sqrt(hdhb_err**2 + hdhb_stellar_err**2) 


    
    #print 'errors on the line ratio from noise, error on the line ratio from uncertainty in stellar absorption' 

    #print haniihb_err, hahb_stellar_err 
    #print hgoiiihb_err, hghb_stellar_err 
    #print hdhberr, hdhb_stellar_err  

    
    hgoiiihb_err_tot=  np.sqrt(hgoiiihb_err**2 + hghb_stellar_err**2) 
    if haniiflux > 0: 
        haniihb_err_tot = np.sqrt(haniihb_err**2 + hahb_stellar_err**2) 


    #### from here we want to return -- 
    ##  observed line ratios (3), 
    ### statistical errors on observed line ratios (3),
    #### best fit model line ratios (3), 
    ### errors including stellar abs error. 




    #### note that, for whatever reason, astropy.io.fits associates the NAXIS keywords with the array dimensions, bakwards. 
    #### hence, NAXIS1 -- hbabs and NAXIS4 = ebv 
    ### when reading the fits file, the likelihood function still has the dimensions it has here. 
    ### hence np.shape(likelihood)[0] != NAXIS1.  instead, np.shape(likelihood)[0] == NAXIS4, and np.shape(likelihood)[3] = NAXIS1. 
    #### use np.shape(likelihood) to get the array dimensions, and these dimensions work appropriately wiht the CRVAL/CRPIX/CDELT specified here. 
    #### do not use NAXIS from the fits header with the CRPIX/CRVAL/CDELT quantities here. 


    
    hdu = fits.ImageHDU(likelihood)
   



    hdr=  fits.Header() 
    hdr['CRVAL1'] = ebv_mod[0] 
    hdr['CRVAL2'] = oh_mod[0] 
    hdr['CRVAL3'] = oiiihg_mod[0] 
    hdr['CRVAL4'] = hb_abs_mod[0] 
    hdr['CRPIX1'] = 0
    hdr['CRPIX2'] = 0 
    hdr['CRPIX3'] = 0
    hdr['CRPIX4'] = 0  
    hdr['CDELT1'] = delta_ebv
    hdr['CDELT2'] = delta_oh 
    hdr['CDELT3'] = delta_oiii 
    hdr['CDELT4'] = delta_hbabs

    hdr['BESTDUST'] = best_dust 
   
    hdr['BESTOH'] = best_oh

    hdr['BESTOIII'] = best_oiiihg 
    hdr['BESTABS'] = best_hbabs 
    
    if haniiflux > 0: 
        hdr['hahbobs'] = haniihb 
    else : 
        hdr['hahbobs'] = -99 

    hdr['hghbobs'] = hgoiiihb

    if hdflux > 0: 
        hdr['hdhbobs']  = hdhb 
    else : 
        hdr['hdhbobs'] = -99 

    if haniiflux > 0: 
        hdr['ehahbobs'] = haniihb_err 
    else : 
        hdr['ehahbobs'] = -99 
    
    hdr['ehghbobs'] = hgoiiihb_err 

    if hdflux > 0 : 
         
        hdr['ehdhbobs'] = hdhb_err 
    else : 
        hdr['ehdhbobs'] = -99

    if haniiflux > 0: 
        hdr['hahbmod'] = haniihb_mod 
    else : 
        hdr['hahbmod'] = -99 
    hdr['hghbmod'] = hgoiiihb_mod 

    if hdflux > 0: 

        hdr['hdhbmod' ] = hdhb_mod 
    else : 
        hdr['hdhbmod'] = -99 

    if haniiflux>0 :
        hdr['ehahbtot'] = haniihb_err_tot 
    else : 
        hdr['ehahbtot'] = -99 
    hdr['ehghbtot'] = hgoiiihb_err_tot 

    if hdflux > 0 :

        hdr['ehdhbtot'] = hdhb_err_tot 
    else : 
        hdr['ehdhbtot'] = -99.0 


    header_hdu = fits.PrimaryHDU(header=hdr)
    
    hdu1 = fits.HDUList([header_hdu, hdu])
    hdu1.writeto(likelihood_func, overwrite=True) 



    
def dust_plot_likelihood(likelihood_func, plot_rootname, showfig = True):   

    hdu = fits.open(likelihood_func) 
    header = hdu[0].header 
    likelihood = hdu[1].data  

    #### re-create ebv_mod, x1, x2, x3, x4, etc.  from the above function 
    ebv_mod = header['CRVAL1'] + np.arange(np.shape(likelihood)[0]) * header['CDELT1']
    oh_mod = header['CRVAL2']  + np.arange(np.shape(likelihood)[1]) * header['CDELT2'] 
    oiiihg_mod = header['CRVAL3'] + np.arange(np.shape(likelihood)[2]) * header['CDELT3'] 
    hb_abs_mod = header['CRVAL4'] + np.arange(np.shape(likelihood)[3]) * header['CDELT4'] 

    x1= np.size(ebv_mod) 
    x2 = np.size(oh_mod) 
    x3= np.size(oiiihg_mod) 
    x4 = np.size(hb_abs_mod) 

    
    best_dust =  header['BESTDUST']
    best_oh = header['BESTOH'] 
    best_oiiihg = header['BESTOIII'] 
    best_hbabs = header['BESTABS'] 



    like_2d_ebv_oh = np.zeros( (x2, x1)) 
    like_2d_ebv_oiii = np.zeros( (x3, x1)) 
    like_2d_oh_oiii = np.zeros( (x3, x2)) 
    like_2d_ebv_stellar = np.zeros((x4, x1)) 
    like_2d_oh_stellar = np.zeros((x4, x2)) 
    like_2d_oiii_stellar = np.zeros((x4, x3)) 
   

    like_1d_ebv = np.zeros(x1) 
    like_1d_oh = np.zeros(x2)
    like_1d_oiii = np.zeros(x3) 
    like_1d_stellar = np.zeros(x4)

   
    ### this looks super confusing but it is right.  
    ### likelihood is indexed i, j, k, m  
    ### for contour, i need, e.g. j, i, not i, j 

    
    
    
    for j in np.arange(x2):
        for i in np.arange(x1): 
            like_2d_ebv_oh[j, i]  = np.sum(likelihood[i, j, ::, ::])  

    for k in np.arange(x3) : 
       for i in np.arange(x1) : 
            like_2d_ebv_oiii[k, i] = np.sum(likelihood[i, ::, k, ::])   

    for m in np.arange(x4): 
        for i in np.arange(x1): 
            like_2d_ebv_stellar[m, i] = np.sum(likelihood[i, ::, ::, m])
               
    for k in np.arange(x3) :
        for j in np.arange(x2): 
            like_2d_oh_oiii[k, j] = np.sum(likelihood[::, j, k, ::])  

    for m in np.arange(x4): 
        for j in np.arange(x2) : 
            like_2d_oh_stellar[m, j] = np.sum(likelihood[::, j, ::, m])  


    for m in np.arange(x4): 
        for k in np.arange(x3): 
            like_2d_oiii_stellar[m, k] = np.sum(likelihood[::, ::, k, m]) 

    for i in np.arange(x1): 
        like_1d_ebv[i] = np.sum(likelihood[i, ::, ::, ::])  

    for j in np.arange(x2) : 
        like_1d_oh[j] = np.sum(likelihood[::, j, ::, ::])  

    for k in np.arange(x3) : 
        like_1d_oiii[k] = np.sum(likelihood[::, ::, k, ::])  
    
    for m in np.arange(x4) : 
        like_1d_stellar[m] = np.sum(likelihood[::, ::, ::, m]) 


    f, axarr = plt.subplots(4, 4,  figsize=(10, 10))
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    axarr[0][0].plot(ebv_mod, like_1d_ebv) 
    axarr[0][0].get_xaxis().set_visible(False) 
    axarr[0][0].set_ylabel('Likelihood', fontsize=12) 


    levels = calc_levels(like_2d_ebv_oh) 
    axarr[1][0].contour(ebv_mod, oh_mod,  like_2d_ebv_oh, levels = levels)
    axarr[1][0].set_ylabel(r'12 + log (O/H)', fontsize=12) 
    axarr[1][0].plot(best_dust, best_oh, 'bo') 
    axarr[1][0].get_xaxis().set_visible(False)

    
    levels = calc_levels(like_2d_ebv_oiii)
    axarr[2][0].contour(ebv_mod, oiiihg_mod,  like_2d_ebv_oiii, levels = levels) 
    axarr[2][0].set_ylabel(r'[OIII] 4363/H$\gamma$', fontsize=12) 
    axarr[2][0].plot(best_dust, best_oiiihg, 'bo') 
    axarr[2][0].get_xaxis().set_visible(False)



    levels = calc_levels(like_2d_ebv_stellar) 
    axarr[3][0].contour(ebv_mod, hb_abs_mod,  like_2d_ebv_stellar,  levels = levels)
    axarr[3][0].set_xlabel('E(B-V)$_{gas}$', fontsize=12)
    axarr[3][0].set_ylabel(r'Hb stellar abs', fontsize=12)
    axarr[3][0].plot(best_dust, best_hbabs, 'bo') 
 


    axarr[0][1].get_xaxis().set_visible(False)
    axarr[0][1].get_yaxis().set_visible(False) 

    axarr[1][1].plot(oh_mod, like_1d_oh) 
    axarr[1][1].get_xaxis().set_visible(False)
    axarr[1][1].get_yaxis().set_visible(False)  
    
    
    levels = calc_levels(like_2d_oh_oiii)
    axarr[2][1].contour(oh_mod, oiiihg_mod,  like_2d_oh_oiii, levels = levels) 
    axarr[2][1].plot(best_oh, best_oiiihg, 'bo') 
    axarr[2][1].get_xaxis().set_visible(False) 
    axarr[2][1].get_yaxis().set_visible(False)
     
    levels = calc_levels(like_2d_oh_stellar)
    axarr[3][1].contour(oh_mod, hb_abs_mod,  like_2d_oh_stellar, levels = levels) 
    axarr[3][1].plot(best_oh, best_hbabs, 'bo') 
    axarr[3][1].get_yaxis().set_visible(False) 
    axarr[3][1].set_xlabel('12 + Log (O/H)', fontsize=12) 
 

    axarr[0][2].get_xaxis().set_visible(False)
    axarr[0][2].get_yaxis().set_visible(False) 
    axarr[1][2].get_xaxis().set_visible(False)
    axarr[1][2].get_yaxis().set_visible(False)
    
    axarr[2][2].plot(oiiihg_mod, like_1d_oiii) 
    axarr[2][2].get_xaxis().set_visible(False)
    axarr[2][2].get_yaxis().set_visible(False)
 
    levels = calc_levels(like_2d_oiii_stellar) 
    axarr[3][2].contour(oiiihg_mod, hb_abs_mod, like_2d_oiii_stellar, levels= levels) 
    axarr[3][2].plot(best_oiiihg, best_hbabs, 'bo') 
    axarr[3][2].get_yaxis().set_visible(False) 
    axarr[3][2].set_xlabel('OIII/Hg', fontsize=12)


    axarr[0][3].get_xaxis().set_visible(False)
    axarr[0][3].get_yaxis().set_visible(False) 
    axarr[1][3].get_xaxis().set_visible(False)
    axarr[1][3].get_yaxis().set_visible(False)
    axarr[2][3].get_xaxis().set_visible(False)
    axarr[2][3].get_yaxis().set_visible(False)


    axarr[3][3].plot(hb_abs_mod, like_1d_stellar) 
    axarr[3][3].set_xlabel('Hb stellar abs.', fontsize=12)
    axarr[3][3].get_yaxis().set_visible(False)


    plt.savefig(plot_rootname + '_confidence.pdf') 

    
    f2, ax = plt.subplots(1,  figsize=(8, 6))
    #plt.subplots_adjust(hspace=0.1, wspace=0.1)
 
 
    ##### plot balmer decrement:
    if header['HAHBOBS'] != -99 : 
        if header['HDHBOBS'] != -99: 

            lam = np.array([4102, 4341, 6564]) 
            ratios = np.array([header['HDHBOBS'], header['HGHBOBS'], header['HAHBOBS']])  
            errors = np.array([header['EHDHBTOT'], header['EHGHBTOT'], header['EHAHBTOT']]) 
            model = np.array([header['HDHBMOD'], header['HGHBMOD'], header['HAHBMOD']])
        else: 
            lam = np.array([ 4341, 6564]) 
            ratios = np.array([ header['HGHBOBS'], header['HAHBOBS']])  
            errors = np.array([ header['EHGHBTOT'], header['EHAHBTOT']]) 
            model = np.array([header['HGHBMOD'], header['HAHBMOD']])

    else : 
        if header['HDHBOBS'] != -99: 
            lam = np.array([4102, 4341]) 
            ratios = np.array([header['HDHBOBS'], header['HGHBOBS']])  
            errors = np.array([header['EHDHBTOT'], header['EHGHBTOT']]) 
            model = np.array([header['HDHBMOD'], header['HGHBMOD']])
        else : 
            lam = np.array([4341]) 
            ratios = np.array([ header['HGHBOBS']])  
            errors = np.array([header['EHGHBTOT']]) 
            model = np.array([header['HGHBMOD']])




  

    d = ax.plot(lam, ratios, 'bo', markeredgecolor = 'b') 
    ax.errorbar(lam, ratios, yerr = errors, fmt = ' ', ecolor = 'b') 


    m = ax.plot(lam, model, 'bo', ms = 10, markeredgecolor = 'b', markerfacecolor = 'None')
    #ax.plot(4865, 1, 'ro')
    ax.set_xlabel(r'$\lambda_{rest}$  $(\AA)$', fontsize=16)
    ax.set_ylabel(r'$F_{line}$  /  $F(H\beta)$', fontsize=16) 

    ax.legend([d[0], m[0]], ['Data', 'Model']) 

    plt.savefig(plot_rootname + '_dustmod.pdf') 


    #### print out data 
    ebv_levs = calc_levels(like_1d_ebv)
    w=np.where(like_1d_ebv > ebv_levs[2])
    ebv_1sig_low = ebv_mod[w][0] 
    ebv_1sig_high = ebv_mod[w][-1]
    w=np.where(like_1d_ebv > ebv_levs[1])
    ebv_2sig_low = ebv_mod[w][0] 
    ebv_2sig_high = ebv_mod[w][-1]
    w=np.where(like_1d_ebv > ebv_levs[0])
    ebv_3sig_low = ebv_mod[w][0] 
    ebv_3sig_high = ebv_mod[w][-1] 

    w=np.where(like_1d_ebv == np.max(like_1d_ebv)) 
    best_dust_marginalized = ebv_mod[w][0] 


    print 'best fitting E(B-V) global, best fitting E(B-V) marginalized,  68% 95%, and 99.7% confidence intervals (marginalized)' 
    print np.round(best_dust,3), np.round(best_dust_marginalized, 3) ,  np.round(ebv_1sig_low, 3), np.round(ebv_1sig_high,3),\
        np.round(ebv_2sig_low, 3), np.round(ebv_2sig_high,3), np.round(ebv_3sig_low,3), np.round(ebv_3sig_high,3) 


    oh_levs=  calc_levels(like_1d_oh) 
    w=np.where(like_1d_oh > oh_levs[2]) 
    oh_1sig_low = oh_mod[w][0] 
    oh_1sig_high = oh_mod[w][-1]  
    w=np.where(like_1d_oh > oh_levs[1]) 
    oh_2sig_low = oh_mod[w][0] 
    oh_2sig_high = oh_mod[w][-1] 
    w=np.where(like_1d_oh > oh_levs[0]) 
    oh_3sig_low = oh_mod[w][0] 
    oh_3sig_high = oh_mod[w][-1] 
    w=np.where(like_1d_oh == np.max(like_1d_oh)) 
    best_oh_marginalized = oh_mod[w][0]

    print 'best fitting 12+log O/H global, best fitting 12 + log O/H marginalized, 68% 95%, and 99.7% confidence intervals (marginalized)'  
    print np.round(best_oh,4), np.round(best_oh_marginalized,4), np.round(oh_1sig_low, 4 ), np.round(oh_1sig_high, 4), np.round(oh_2sig_low, 4),\
            np.round(oh_2sig_high, 4),  np.round(oh_3sig_low, 4), np.round(oh_3sig_high, 4) 



    oiiihg_levs = calc_levels(like_1d_oiii) 
    w=np.where(like_1d_oiii > oiiihg_levs[2]) 
    oiii_1sig_low = oiiihg_mod[w][0] 
    oiii_1sig_high = oiiihg_mod[w][-1]
    w=np.where(like_1d_oiii > oiiihg_levs[1]) 
    oiii_2sig_low = oiiihg_mod[w][0] 
    oiii_2sig_high = oiiihg_mod[w][-1]
    w=np.where(like_1d_oiii > oiiihg_levs[0]) 
    oiii_3sig_low = oiiihg_mod[w][0] 
    oiii_3sig_high = oiiihg_mod[w][-1]  
    w= np.where(like_1d_oiii == np.max(like_1d_oiii)) 
    best_oiiihg_marginalized = oiiihg_mod[w][0] 

    print 'best fitting oiii/hg global, best fitting oiii/hg marginalized, 68% 95%, and 99.7% confidence intervals (marginalized)'  
    print np.round(best_oiiihg, 4), np.round(best_oiiihg_marginalized, 4), np.round(oiii_1sig_low, 4), np.round(oiii_1sig_high, 4),\
            np.round(oiii_2sig_low, 4), np.round(oiii_2sig_high, 4), np.round(oiii_3sig_low, 4), np.round(oiii_3sig_high, 4) 


    hbabs_levs = calc_levels(like_1d_stellar) 
    w=np.where(like_1d_stellar > hbabs_levs[2]) 
    hbabs_1sig_low  = hb_abs_mod[w][0] 
    hbabs_1sig_high = hb_abs_mod[w][-1] 
    w=np.where(like_1d_stellar > hbabs_levs[1]) 
    hbabs_2sig_low  = hb_abs_mod[w][0] 
    hbabs_2sig_high = hb_abs_mod[w][-1] 
    w=np.where(like_1d_stellar > hbabs_levs[0]) 
    hbabs_3sig_low  = hb_abs_mod[w][0] 
    hbabs_3sig_high = hb_abs_mod[w][-1]  
    w=np.where(like_1d_stellar == np.max(like_1d_stellar)) 
    best_hbabs_marginalized = hb_abs_mod[w][0] 


    print 'best fitting Hb stellar abs global, best fitting Hb stellar abs marginalized, 68% 95%, and 99.7% confidence intervals (marginalized)'    
    print np.round(best_hbabs, 3), np.round(best_hbabs_marginalized, 3), np.round(hbabs_1sig_low, 3), np.round(hbabs_1sig_high,3),\
            np.round(hbabs_2sig_low, 3), np.round(hbabs_2sig_high, 3), np.round(hbabs_3sig_low,3), np.round(hbabs_3sig_high, 3) 







    
    if showfig == True : 


        plt.show()

    plt.close()



    
    



        

    














