from wisp_analysis import * 

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



def dust_calc_likelihood(input_meas, likelihood_func): 
    
    #### input meas comes from measure_stack 
    ### likelihood func is a fits file

    
    ####   intrinsic balmer decrement for 10,000K 
    hahb_int = 2.86 
    hghb_int = 0.468  
    hdhb_int = 0.259 
     
     
    #fluxes = [] 
    #errors = [] 
    #f = open(input_meas) 
    #for line in f: 
    #     x = line.split() 
    #     errors.append(float(x[-1])) 
    #     fluxes.append(float(x[-2])) 
    #f.close()
      ##### fluxes/errors are indexed:  
      ## 0   hanii 
      ## 1   hb 
      ## 2   hg 
      ## 3   hd
      ## 4  OIII 
      ## 5  OII 
      ## 6  SII 
      ## 7  HeI 
      ## 8  OI 
      ## 9  NeIII 
    
    
    #fluxes=  np.array(fluxes) 
    #errors = np.array(errors) 

    meas = asciitable.read(input_meas, format = 'fixed_width') 
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

    he1flux = fluxes[7] 
    he1err = errors[7] 

    oiflux = fluxes[8] 
    oierr = errors[8] 

    neiiiflux = fluxes[9] 
    neiiierr = errors[9]  

   
    #### stuff needed for ews errors 
    #if haniiflux > 0: 
    #    haniierr_cgs= haniiflux_cgs * np.sqrt( (haniierr/haniiflux)**2 + (oiiierr/oiiiflux)**2) 
       # haew_err = haniiew * np.sqrt( (haniierr_cgs/haniiflux_cgs)**2 + (hanii_cont_err/ hanii_cont)**2) 
    if haniiflux < 0 :  
    #    haniierr_cgs = 1e12 
         haew_err = 1e12 

   
    #hberr_cgs= hbflux_cgs * np.sqrt( (hberr/hbflux)**2 + (oiiierr/oiiiflux)**2) 
    #hbew_err = hbew * np.sqrt( (hberr_cgs/hbflux_cgs)**2 + (hb_cont_err/ hb_cont)**2)  

    #hgerr_cgs= hgflux_cgs * np.sqrt( (hgoiiierr/hgoiiiflux)**2 + (oiiierr/oiiiflux)**2) 
    #hgoiiiew_err = hgoiiiew * np.sqrt((hgerr_cgs/hgflux_cgs)**2 + (hg_cont_err/ hg_cont)**2)  

    #hderr_cgs= hdflux_cgs * np.sqrt( (hderr/hdflux)**2 + (oiiierr/oiiiflux)**2) 
    #hdew_err = hdew * np.sqrt( (hderr_cgs/hdflux_cgs)**2 + (hd_cont_err/ hd_cont)**2) 


    #### note hg ==  hg + oiii 4363, and this component may be significant. 

    ### calculate observed line ratios :

    if haniiflux > 0 : 
        haniihb = haniiflux / hbflux  
        haniihb_err = haniihb * np.sqrt( (haniierr/haniiflux)**2 + (hberr/hbflux)**2)  
    else :
        haniihb = 0
        haniihb_err = 1e12 


    hgoiiihb = hgoiiiflux / hbflux 
    hgoiiihb_err = hgoiiihb * np.sqrt( (hgoiiierr/hgoiiiflux)**2 + (hberr/hbflux)**2)  

    hdhb = hdflux/hbflux 
    hdhb_err  = hdhb * np.sqrt((hderr/ hdflux)**2 + (hberr/hbflux)**2) 
     
    #print 'Observed ratios:  HaNII/Hb, HgOIII/Hb, Hb/Hd' 
    #print haniihb, hgoiiihb, hdhb 


    delta_ebv = 0.01 
    delta_nii = 0.01 
    delta_oiii = 0.01 
    delta_hbabs = 0.3 
    
    ebv_mod = np.arange(0, 0.5, delta_ebv) 
    if haniiflux > 0: 
        niiha_mod = np.arange(0, 0.3, delta_nii)  #### ratio ### try a guess for now.  eventually iterate over niiha. 
    else : 
        niiha_mod = np.arange(0, 0.1, 0.05)  #### set this to something small to speed up but not break code.  
    oiiihg_mod = np.arange(0, 0.5, delta_oiii)    #### ratio ### try a guess for now, eventually iterate. 
    hb_abs_mod = np.arange(0, 7, delta_hbabs)  #### rest-frame stellar absorption ew, to be compared to the rest frame emission ews. 
    
    ### also need to include stellar absorption and iterate over that. 




    k_array = k_lambda(np.array([6564., 4861., 4341., 4102.]))
    k_ha = k_array[0]
    k_hb=  k_array[1] 
    k_hg = k_array[2]  
    k_hd = k_array[3]

    x1= np.size(ebv_mod) 
    x2 = np.size(niiha_mod) 
    x3= np.size(oiiihg_mod) 
    x4 = np.size(hb_abs_mod) 

    likelihood = np.zeros( (x1, x2, x3, x4))

    for i in np.arange(x1):
        for j in np.arange(x2): 
            for k in np.arange(x3):
                for m in np.arange(x4): 

                    if haniiflux > 0: 
                        hahb_mod = hahb_int * 10**(ebv_mod[i] / 2.5 * (k_hb - k_ha)) 
                        
                    hghb_mod = hghb_int * 10**(ebv_mod[i] / 2.5 * (k_hb - k_hg)) 
                    hdhb_mod = hdhb_int * 10**(ebv_mod[i] / 2.5 * (k_hb - k_hd))


                    if haniiflux > 0 : 
                         ha_stellar_corr = (haniiew + hb_abs_mod[m]/ 1.5) / haniiew  #### divide model by this number to match data 
                    hb_stellar_corr = (hbew + hb_abs_mod[m]) / hbew 
                    hg_stellar_corr = (hgoiiiew + hb_abs_mod[m]) /hgoiiiew 
                    hd_stellar_corr = (hdew  + hb_abs_mod[m]) / hdew


                    if haniiflux > 0  : 
                        ha_stellar_corr_err = hb_abs_mod[m]/1.5 * haew_err/haniiew**2  
                    hb_stellar_corr_err = hb_abs_mod[m]   * hbew_err /hbew**2 
                    hg_stellar_corr_err = hb_abs_mod[m] * hgoiiiew_err/hgoiiiew**2 
                    hd_stellar_corr_err = hb_abs_mod[m] * hdew_err/ hdew**2


                    if haniiflux > 0 : 
                        hahb_mod = hahb_mod * hb_stellar_corr / ha_stellar_corr ### correction is (Ha / ha_corr) / (Hb/Hb_corr)  
                    hghb_mod = hghb_mod * hb_stellar_corr / hg_stellar_corr 
                    hdhb_mod = hdhb_mod * hb_stellar_corr / hd_stellar_corr  
                    
                    if haniiflux > 0 : 
                        hahb_stellar_err = hb_stellar_corr / ha_stellar_corr * np.sqrt( (ha_stellar_corr_err/ha_stellar_corr)**2 + (hb_stellar_corr_err/ hb_stellar_corr)**2) 
                    hghb_stellar_err = hb_stellar_corr/hg_stellar_corr * np.sqrt( (hg_stellar_corr_err/hg_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2) 
                    hdhb_stellar_err = hb_stellar_corr / hd_stellar_corr * np.sqrt( (hd_stellar_corr_err/hd_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2) 

                    if haniiflux > 0: 
                        haniihb_err_tot = np.sqrt(haniihb_err**2 + hahb_stellar_err**2) 
                    hgoiiihb_err_tot = np.sqrt(hgoiiihb_err**2 + hghb_stellar_err**2) 
                    hdhb_err_tot = np.sqrt(hdhb_err**2 + hdhb_stellar_err**2) 



                    if haniiflux > 0:
                        haniihb_mod = hahb_mod  * (1 + niiha_mod[j]) 
                    hgoiiihb_mod = hghb_mod * (1 + oiiihg_mod[k]) 

                    if haniiflux > 0 : 
                        like1 = np.exp(-1 * (haniihb - haniihb_mod)**2 / (2 * haniihb_err_tot**2)) 
                    else: 
                        like1 = 1.0 

                    like2 = np.exp(-1 * (hgoiiihb - hgoiiihb_mod)**2 / (2 * hgoiiihb_err_tot**2)) 
                    like3 = np.exp(-1 * (hdhb  - hdhb_mod)**2 / (2 * hdhb_err_tot**2)) 

                    likelihood[i, j, k, m] = like1 * like2 * like3  
    


    likelihood = likelihood / np.sum(likelihood)
    
        
    ##### evaluate the best fitting model: 
    w=np.where(likelihood == np.max(likelihood)) 
    best_dust =  ebv_mod[w[0]][0]
    best_niiha = niiha_mod[w[1]][0]
    best_oiiihg = oiiihg_mod[w[2]][0]
    best_hbabs = hb_abs_mod[w[3]][0]

    if haniiflux > 0: 
        hahb_mod = hahb_int * 10**(best_dust / 2.5 * (k_hb - k_ha)) 
    hghb_mod = hghb_int * 10**(best_dust / 2.5 * (k_hb - k_hg)) 
    hdhb_mod = hdhb_int * 10**(best_dust / 2.5 * (k_hb - k_hd))  

    #print 'models with dust only, a, g, d' 
    #print hahb_mod, hghb_mod, hdhb_mod 
   
    if haniiflux > 0 : 
        ha_stellar_corr = (haniiew + best_hbabs/ 1.5) / haniiew  #### divide model by this number to match data 
    hb_stellar_corr = (hbew + best_hbabs) / hbew 
    hg_stellar_corr = (hgoiiiew + best_hbabs) /hgoiiiew 
    hd_stellar_corr = (hdew  + best_hbabs) / hdew 

    #print 'stellar corrections, a, b, g, d' 
    #print ha_stellar_corr, hb_stellar_corr, hg_stellar_corr, hd_stellar_corr
 
    if haniiflux > 0: 
        ha_stellar_corr_err = best_hbabs/1.5 * haew_err/haniiew**2  
    hb_stellar_corr_err = best_hbabs * hbew_err /hbew**2 
    hg_stellar_corr_err = best_hbabs * hgoiiiew_err/hgoiiiew**2 
    hd_stellar_corr_err = best_hbabs * hdew_err/ hdew**2  

    #print  'uncertainties on stellar corr'  
    #print ha_stellar_corr_err, hb_stellar_corr_err, hg_stellar_corr_err, hd_stellar_corr_err 

    if haniiflux > 0: 
         hahb_mod = hahb_mod * hb_stellar_corr / ha_stellar_corr ### correction is (Ha / ha_corr) / (Hb/Hb_corr)  
    hghb_mod = hghb_mod * hb_stellar_corr / hg_stellar_corr 
    hdhb_mod = hdhb_mod * hb_stellar_corr / hd_stellar_corr  
 
    if haniiflux > 0 :
        haniihb_mod = hahb_mod  * (1 + best_niiha) 
    hgoiiihb_mod = hghb_mod * (1 + best_oiiihg)
    
    #print 'best models with dust, stellar absorption, nii, and oiii4363, a, g, d' 
    #print haniihb_mod, hgoiiihb_mod, hdhb_mod 

    if haniiflux > 0: 
        hahb_stellar_err = hb_stellar_corr / ha_stellar_corr * np.sqrt( (ha_stellar_corr_err/ha_stellar_corr)**2 + (hb_stellar_corr_err/ hb_stellar_corr)**2) 
    hghb_stellar_err = hb_stellar_corr/hg_stellar_corr * np.sqrt( (hg_stellar_corr_err/hg_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2) 
    hdhb_stellar_err = hb_stellar_corr / hd_stellar_corr * np.sqrt( (hd_stellar_corr_err/hd_stellar_corr)**2 + (hb_stellar_corr_err/hb_stellar_corr)**2)  

    
    #print 'errors on the line ratio from noise, error on the line ratio from uncertainty in stellar absorption' 

    #print haniihb_err, hahb_stellar_err 
    #print hgoiiihb_err, hghb_stellar_err 
    #print hdhberr, hdhb_stellar_err  

    
    hdhb_err_tot=  np.sqrt(hdhb_err**2 + hdhb_stellar_err**2) 
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
    hdr['CRVAL2'] = niiha_mod[0] 
    hdr['CRVAL3'] = oiiihg_mod[0] 
    hdr['CRVAL4'] = hb_abs_mod[0] 
    hdr['CRPIX1'] = 0
    hdr['CRPIX2'] = 0 
    hdr['CRPIX3'] = 0
    hdr['CRPIX4'] = 0  
    hdr['CDELT1'] = delta_ebv
    hdr['CDELT2'] = delta_nii 
    hdr['CDELT3'] = delta_oiii 
    hdr['CDELT4'] = delta_hbabs

    hdr['BESTDUST'] = best_dust 
    if haniiflux > 0: 
        hdr['BESTNII'] = best_niiha
    else : 
        hdr['BESTNII'] = -99 

    hdr['BESTOIII'] = best_oiiihg 
    hdr['BESTABS'] = best_hbabs 
    
    if haniiflux > 0: 
        hdr['hahbobs'] = haniihb 
    else : 
        hdr['hahbobs'] = -99 

    hdr['hghbobs'] = hgoiiihb 
    hdr['hdhbobs']  = hdhb 

    if haniiflux > 0: 
        hdr['ehahbobs'] = haniihb_err 
    else : 
        hdr['ehahbobs'] = -99 
    
    hdr['ehghbobs'] = hgoiiihb_err 
    hdr['ehdhbobs'] = hdhb_err 

    if haniiflux > 0: 
        hdr['hahbmod'] = haniihb_mod 
    else : 
        hdr['hahbmod'] = -99 
    hdr['hghbmod'] = hgoiiihb_mod 
    hdr['hdhbmod' ] = hdhb_mod 


    if haniiflux>0 :
        hdr['ehahbtot'] = haniihb_err_tot 
    else : 
        hdr['ehahbtot'] = -99 
    hdr['ehghbtot'] = hgoiiihb_err_tot 
    hdr['ehdhbtot'] = hdhb_err_tot 


    header_hdu = fits.PrimaryHDU(header=hdr)
    
    hdu1 = fits.HDUList([header_hdu, hdu])
    hdu1.writeto(likelihood_func, overwrite=True) 



    
def dust_plot_likelihood(likelihood_func, plot_rootname):   

    hdu = fits.open(likelihood_func) 
    header = hdu[0].header 
    likelihood = hdu[1].data  

    #### re-create ebv_mod, x1, x2, x3, x4, etc.  from the above function 
    ebv_mod = header['CRVAL1'] + np.arange(np.shape(likelihood)[0]) * header['CDELT1']
    niiha_mod = header['CRVAL2']  + np.arange(np.shape(likelihood)[1]) * header['CDELT2'] 
    oiiihg_mod = header['CRVAL3'] + np.arange(np.shape(likelihood)[2]) * header['CDELT3'] 
    hb_abs_mod = header['CRVAL4'] + np.arange(np.shape(likelihood)[3]) * header['CDELT4'] 

    x1= np.size(ebv_mod) 
    x2 = np.size(niiha_mod) 
    x3= np.size(oiiihg_mod) 
    x4 = np.size(hb_abs_mod) 

    
    best_dust =  header['BESTDUST']
    best_niiha = header['BESTNII'] 
    best_oiiihg = header['BESTOIII'] 
    best_hbabs = header['BESTABS'] 



    like_2d_ebv_nii = np.zeros( (x2, x1)) 
    like_2d_ebv_oiii = np.zeros( (x3, x1)) 
    like_2d_nii_oiii = np.zeros( (x3, x2)) 
    like_2d_ebv_stellar = np.zeros((x4, x1)) 
    like_2d_nii_stellar = np.zeros((x4, x2)) 
    like_2d_oiii_stellar = np.zeros((x4, x3)) 
   

    like_1d_ebv = np.zeros(x1) 
    like_1d_nii = np.zeros(x2)
    like_1d_oiii = np.zeros(x3) 
    like_1d_stellar = np.zeros(x4)

   
    ### this looks super confusing but it is right.  
    ### likelihood is indexed i, j, k, m  
    ### for contour, i need, e.g. j, i, not i, j 

    
    
    
    for j in np.arange(x2):
        for i in np.arange(x1): 
            like_2d_ebv_nii[j, i]  = np.sum(likelihood[i, j, ::, ::])  

    for k in np.arange(x3) : 
       for i in np.arange(x1) : 
            like_2d_ebv_oiii[k, i] = np.sum(likelihood[i, ::, k, ::])   

    for m in np.arange(x4): 
        for i in np.arange(x1): 
            like_2d_ebv_stellar[m, i] = np.sum(likelihood[i, ::, ::, m])
               
    for k in np.arange(x3) :
        for j in np.arange(x2): 
            like_2d_nii_oiii[k, j] = np.sum(likelihood[::, j, k, ::])  

    for m in np.arange(x4): 
        for j in np.arange(x2) : 
            like_2d_nii_stellar[m, j] = np.sum(likelihood[::, j, ::, m])  


    for m in np.arange(x4): 
        for k in np.arange(x3): 
            like_2d_oiii_stellar[m, k] = np.sum(likelihood[::, ::, k, m]) 

    for i in np.arange(x1): 
        like_1d_ebv[i] = np.sum(likelihood[i, ::, ::, ::])  

    for j in np.arange(x2) : 
        like_1d_nii[j] = np.sum(likelihood[::, j, ::, ::])  

    for k in np.arange(x3) : 
        like_1d_oiii[k] = np.sum(likelihood[::, ::, k, ::])  
    
    for m in np.arange(x4) : 
        like_1d_stellar[m] = np.sum(likelihood[::, ::, ::, m]) 


    f, axarr = plt.subplots(4, 4,  figsize=(10, 10))
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    axarr[0][0].plot(ebv_mod, like_1d_ebv) 
    axarr[0][0].get_xaxis().set_visible(False) 
    axarr[0][0].set_ylabel('Likelihood', fontsize=12) 


    levels = calc_levels(like_2d_ebv_nii) 
    axarr[1][0].contour(ebv_mod, niiha_mod,  like_2d_ebv_nii, levels = levels)
    axarr[1][0].set_ylabel(r'[NII]/H$\alpha$', fontsize=12) 
    axarr[1][0].plot(best_dust, best_niiha, 'bo') 
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

    axarr[1][1].plot(niiha_mod, like_1d_nii) 
    axarr[1][1].get_xaxis().set_visible(False)
    axarr[1][1].get_yaxis().set_visible(False)  
    
    
    levels = calc_levels(like_2d_nii_oiii)
    axarr[2][1].contour(niiha_mod, oiiihg_mod,  like_2d_nii_oiii, levels = levels) 
    axarr[2][1].plot(best_niiha, best_oiiihg, 'bo') 
    axarr[2][1].get_xaxis().set_visible(False) 
    axarr[2][1].get_yaxis().set_visible(False)
     
    levels = calc_levels(like_2d_nii_stellar)
    axarr[3][1].contour(niiha_mod, hb_abs_mod,  like_2d_nii_stellar, levels = levels) 
    axarr[3][1].plot(best_niiha, best_hbabs, 'bo') 
    axarr[3][1].get_yaxis().set_visible(False) 
    axarr[3][1].set_xlabel('NII/Ha', fontsize=12) 
 

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
    if best_niiha != -99 : 
        lam = np.array([4102, 4341, 6564]) 
        ratios = np.array([header['HDHBOBS'], header['HGHBOBS'], header['HAHBOBS']])  
        errors = np.array([header['EHDHBTOT'], header['EHGHBTOT'], header['EHAHBTOT']]) 
        model = np.array([header['HDHBMOD'], header['HGHBMOD'], header['HAHBMOD']])
    else : 
        lam = np.array([4102, 4341]) 
        ratios = np.array([header['HDHBOBS'], header['HGHBOBS']])  
        errors = np.array([header['EHDHBTOT'], header['EHGHBTOT']]) 
        model = np.array([header['HDHBMOD'], header['HGHBMOD']])


  

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


    if best_niiha != -99 : 
        niiha_levs=  calc_levels(like_1d_nii) 
        w=np.where(like_1d_nii > niiha_levs[2]) 
        nii_1sig_low = niiha_mod[w][0] 
        nii_1sig_high = niiha_mod[w][-1]  
        w=np.where(like_1d_nii > niiha_levs[1]) 
        nii_2sig_low = niiha_mod[w][0] 
        nii_2sig_high = niiha_mod[w][-1] 
        w=np.where(like_1d_nii > niiha_levs[0]) 
        nii_3sig_low = niiha_mod[w][0] 
        nii_3sig_high = niiha_mod[w][-1] 
        w=np.where(like_1d_nii == np.max(like_1d_nii)) 
        best_nii_marginalized = niiha_mod[w][0]

        print 'best fitting nii/ha global, best fitting nii/ha marginalized, 68% 95%, and 99.7% confidence intervals (marginalized)'  
        print np.round(best_niiha,4), np.round(best_nii_marginalized,4), np.round(nii_1sig_low, 4 ), np.round(nii_1sig_high, 4), np.round(nii_2sig_low, 4),\
            np.round(nii_2sig_high, 4),  np.round(nii_3sig_low, 4), np.round(nii_3sig_high, 4) 



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







    


    plt.show()   


    
    



        

    














