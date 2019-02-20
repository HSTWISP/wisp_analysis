from wisp_analysis import *
from mpfit import *

import astropy.units as u
from astropy.cosmology import Planck13 as cosmo


def stack_emline_model(pars, x):
    sigma_oii = pars[0]
    sigma_neiii = pars[1] 
    sigma_hd = pars[2] 
    sigma_hg = pars[3]
    sigma_oiiihb = pars[4] 
    sigma_he1_5876 = pars[5] 
    sigma_oi = pars[6] 
    sigma_hasii = pars[7] 



    ha_amp = pars[8] 
    hb_amp = pars[9] 
    hg_amp = pars[10] 
    hd_amp = pars[11]
    oiii_5007_amp = pars[12] 
    oii_amp = pars[13] 
    sii_amp = pars[14] 
    he1_5876_amp = pars[15] 
    oi_amp = pars[16] 
    neiii_amp = pars[17]
    nii_6583_amp = pars[18] 
    he1_6678_amp = pars[19]

    hasii_shift = pars[20] 
    oi_shift = pars[21] 
    he1_5876_shift = pars[22] 
    oiii_shift = pars[23] 
    hb_shift = pars[24]
    hg_shift = pars[25] 
    hd_shift = pars[26] 
    neiii_shift = pars[27] 
    oii_shift = pars[28]

    c_3800_4200 = pars[29] 
    c_5200_6400 = pars[30]
    s_5200_6500  = pars[31]
    cnorm = pars[32] 


    cont = np.zeros(np.size(x)) 
    w=np.where( (x > 3800) & (x < 4600)) 
    cont[w] = c_3800_4200 
    w=np.where( (x> 5400) & (x < 6450)) 
    cont[w] = c_5200_6400 + x[w] * s_5200_6500

    cont = cont + cnorm

  


    model = ha_amp * gaussian(x, 6564.6 + hasii_shift, sigma_hasii) + \
            hb_amp  * gaussian(x, 4862.7 + hb_shift, sigma_oiiihb) + \
            hg_amp * gaussian(x, 4341.7 + hg_shift, sigma_hg) + \
            hd_amp * gaussian(x, 4102.9 + hd_shift, sigma_hd) + \
            oiii_5007_amp * gaussian(x, 5008. + oiii_shift, sigma_oiiihb) + \
            oiii_5007_amp/2.98 * gaussian(x, 4960. + oiii_shift, sigma_oiiihb) +\
            oii_amp * gaussian(x, 3728. + oii_shift, sigma_oii) +\
            sii_amp * gaussian(x, 6725.  + hasii_shift, sigma_hasii) +\
            he1_5876_amp  * gaussian(x, 5877.2 + he1_5876_shift, sigma_he1_5876) +\
            oi_amp *  gaussian(x, 6302. + oi_shift, sigma_oi) + \
            oi_amp /3 * gaussian(x, 6365.5 + oi_shift,sigma_oi) + \
            neiii_amp * gaussian(x, 3870 + neiii_shift, sigma_neiii) + \
            nii_6583_amp * gaussian(x, 6585.23 + hasii_shift, sigma_hasii) + \
            nii_6583_amp / 3 * gaussian(x, 6550. + hasii_shift, sigma_hasii)+\
            he1_6678_amp  * gaussian(x, 6680. + hasii_shift, sigma_hasii) +\
            cont
            
        

    return model


    

def model_resid_stack(pars, fjac=None, lam = None, flux = None, err = None):
    model = stack_emline_model(pars, lam) 
    resid = (flux- model) / err 
   
    status = 0 
    return [status, resid]

    



def measure_stack(input_stack, input_masterlist, output_meas, output_fig, zmax = 2.3, showfig = False): 


    lam_max = 17200/(1+zmax) 

    tab= asciitable.read(input_stack) 
    lam = tab['lam']
    flux = tab['flux_norm']
    err = tab['err']
    median_cont = tab['median_cont'] 
    err_cont = tab['err_cont'] 

    w=np.where(lam < lam_max) 
    w=w[0] 
    lam = lam[w]
    flux = flux[w] 
    err = err[w] 
    median_cont = median_cont[w] 
    err_cont = err_cont[w] 


    w=np.where(err > 0) 
    w=w[0] 
    lam = lam[w]
    flux = flux[w] 
    err = err[w]
    median_cont = median_cont[w] 
    err_cont = err_cont[w]


    w=np.where(err_cont > 0) 
    w=w[0] 
    lam = lam[w]
    flux = flux[w] 
    err = err[w]
    median_cont = median_cont[w] 
    err_cont = err_cont[w]


  
    pguess= np.zeros(33) 
    pguess[0] = 25 
    pguess[1] = 25 
    pguess[2] = 25
    pguess[3] = 25 
    pguess[4] = 25
    pguess[5] = 25 
    pguess[6] = 25 
    pguess[7] = 25 

    pguess[8] = 0.003
    pguess[9] = 0.001
    pguess[10] = 0.0008
    pguess[11] = 0.0003
    pguess[12] = 0.0055
    pguess[13] = 0.015
    pguess[14] = 0.0004
    pguess[15] = 0.001
    pguess[16] = 0.0005
    pguess[17] = 0.0003 
    pguess[18] = 0.0000   #### set nii and he1 to zero, even though they are in the model and can be added. 
    pguess[19] = 0.0000
   
    #### wavelength shifts, in rest wavelengths
    pguess[20]  = 0.
    pguess[21]  = 0.
    pguess[22] = 0.
    pguess[23] = 0.
    pguess[24] = 0. 
    pguess[25] = 0.0 
    pguess[26] = 0.0
    pguess[27] = 0.0
    pguess[28]  = 0.0 

    pguess[29] = 0.00001
    pguess[30]  = -0.01
    pguess[31] = 0.0001 
    pguess[32] = 0.0 



    


    npars = len(pguess) 
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} 
              for i in range(npars)]
    for i in range(npars): parinfo[i]['value'] = pguess[i]  
   
 
    
    parinfo[18]['fixed'] = 1
    parinfo[19]['fixed'] = 1

    #### fwhm positive 
    parinfo[0]['limited'][0] = 1
    parinfo[0]['limits'][0] = 3
    parinfo[0]['limits'][1] = 100
    parinfo[1]['limited'][0] = 1
    parinfo[1]['limits'][0] = 3
    parinfo[1]['limits'][1] = 100 
    parinfo[2]['limited'][0] = 1
    parinfo[2]['limits'][0] = 3
    parinfo[2]['limits'][1] = 100 
    parinfo[3]['limited'][0] = 1
    parinfo[3]['limits'][0] = 3
    parinfo[3]['limits'][1] = 100 
    parinfo[4]['limited'][0] = 1
    parinfo[4]['limits'][0] = 3
    parinfo[4]['limits'][1] = 100 
    parinfo[5]['limited'][0] = 1
    parinfo[5]['limits'][0] = 3
    parinfo[5]['limits'][1] = 100 
    parinfo[6]['limited'][0] = 1
    parinfo[6]['limits'][0] = 3
    parinfo[6]['limits'][1] = 100 
    parinfo[7]['limited'][0] = 1
    parinfo[7]['limits'][0] = 3
    parinfo[7]['limits'][1] = 100 


         
    #### amplitiudes positive     
    parinfo[8]['limited'][0] = 1
    parinfo[8]['limits'][0] = 0
    parinfo[9]['limited'][0] = 1
    parinfo[9]['limits'][0] = 0
    parinfo[10]['limited'][0] = 1
    parinfo[10]['limits'][0] = 0
    parinfo[11]['limited'][0] = 1
    parinfo[11]['limits'][0] = 0
    parinfo[12]['limited'][0] = 1
    parinfo[12]['limits'][0] = 0
    parinfo[13]['limited'][0] = 1
    parinfo[13]['limits'][0] = 0
    parinfo[14]['limited'][0] = 1
    parinfo[14]['limits'][0] = 0
    parinfo[15]['limited'][0] = 1
    parinfo[15]['limits'][0] = 0
    parinfo[16]['limited'][0] = 1
    parinfo[16]['limits'][0] = 0
    parinfo[17]['limited'][0] = 1
    parinfo[17]['limits'][0] = 0




    ### wave shifts
    parinfo[20]['limited']  = [1, 1] 
    parinfo[20]['limits'] = [-10, 10]
    parinfo[21]['limited']  = [1, 1] 
    parinfo[21]['limits'] = [-10, 10]
    parinfo[22]['limited']  = [1, 1] 
    parinfo[22]['limits'] = [-10, 10]
    parinfo[23]['limited']  = [1, 1] 
    parinfo[23]['limits'] = [-10, 10]
    parinfo[24]['limited']  = [1, 1] 
    parinfo[24]['limits'] = [-10, 10]
    parinfo[25]['limited']  = [1, 1] 
    parinfo[25]['limits'] = [-10, 10]
    parinfo[26]['limited']  = [1, 1] 
    parinfo[26]['limits'] = [-10, 10]
    parinfo[27]['limited']  = [1, 1] 
    parinfo[27]['limits'] = [-10, 10]
    parinfo[28]['limited']  = [1, 1] 
    parinfo[28]['limits'] = [-10, 10]


    parinfo[29]['limited']  = [1, 1] 
    parinfo[29]['limits'] = [-0.01, 0.01]
    parinfo[30]['limited']  = [1, 1] 
    parinfo[30]['limits'] = [-0.01, 0.01 ]
    parinfo[30]['limited']  = [1, 1] 
    parinfo[30]['limits'] = [-0.01, 0.01 ]

    parinfo[32]['fixed'] = 1  ### fix continuum normalization to zero for continuum subtracted spectra. 



    #### do fit and evaluate model 
    fa = {'lam':lam, 'flux':flux, 'err':err} 
    out = mpfit(model_resid_stack, pguess, functkw=fa, parinfo = parinfo, quiet=True) 
    model_fit = stack_emline_model(out.params, lam)

    

    #### redo model fit on continuum normalized spectrum
    ### update the pguess based on what the previous out.params gave. 
    pguess[0:8] = out.params[0:8] 
    pguess[8:20]  = out.params[8:20] * 300.   ### scale since the spectra are not line flux normalized. 
    pguess[20:29] = out.params[20:29]
    pguess[10]  = 0.25
    #pguess[8] = 0.003
    #pguess[9] = 0.001
    #pguess[10] = 0.0008
    #pguess[11] = 0.0003
    #pguess[12] = 0.0055
    #pguess[13] = 0.015
    #pguess[14] = 0.0004
    #pguess[15] = 0.001
    #pguess[16] = 0.0005
    #pguess[17] = 0.0003 
    #pguess[18] = 0.0000   #### set nii and he1 to zero, even though they are in the model and can be added. 
    #pguess[19] = 0.0000
    ### set continuum to near 1.
     
    #pguess[29] = 0.00000
    #pguess[30]  = -0.00
    #pguess[31] = 0.0000 
    #parinfo[29]['fixed'] = 1
    #parinfo[30]['fixed']  = 1
    #parinfo[31]['fixed'] = 1
 

    pguess[32] = 1.0 
    parinfo[32]['limited'] = [1, 1] 
    parinfo[32]['limits'] = [0.9, 1.1] 
     
    fa2 = {'lam':lam, 'flux':median_cont, 'err':err_cont} 
    out2 = mpfit(model_resid_stack, pguess, functkw=fa2, parinfo = parinfo, quiet=True) 
    model_fit2 = stack_emline_model(out2.params, lam)
 

     

    




    ### gather line fluxes 
    ### ha  
    ### note that nii  and he1 are set to zero.  This means that this is Ha + NII
    ### and I'm not sure if He 1 is included in the flux or resolved, or something in between 
    ha_flux = np.sqrt(2  * math.pi) *  out.params[7] * out.params[8]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[7][8] / (out.params[7] * out.params[8])
    ha_err =  ha_flux * np.sqrt( (out.perror[7]/out.params[7])**2 + (out.perror[8]/out.params[8])**2 + covar_term) 
     
    ### the "flux" in the continuum normalized model is really the EW, since Flam_cont = 1. 
    ha_ew = np.sqrt(2  * math.pi) *  out2.params[7] * out2.params[8]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[7][8] / (out2.params[7] * out2.params[8])
    ha_ew_err =  ha_ew * np.sqrt( (out2.perror[7]/out2.params[7])**2 + (out2.perror[8]/out2.params[8])**2 + covar_term) 
     

    ##hb 
    hb_flux = np.sqrt(2  * math.pi) *  out.params[4] * out.params[9] ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[4][9] / (out.params[4] * out.params[9])
    hb_err =  hb_flux * np.sqrt( (out.perror[4]/out.params[4])**2 + (out.perror[9]/out.params[9])**2 + covar_term)  

    hb_ew = np.sqrt(2  * math.pi) *  out2.params[4] * out2.params[9] ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[4][9] / (out2.params[4] * out2.params[9])
    hb_ew_err =  hb_ew * np.sqrt( (out2.perror[4]/out2.params[4])**2 + (out2.perror[9]/out2.params[9])**2 + covar_term)

    ##hg 
    hg_flux = np.sqrt(2  * math.pi) *  out.params[3] * out.params[10]     ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[3][10] / (out.params[3] * out.params[10])
    hg_err =  hg_flux * np.sqrt( (out.perror[3]/out.params[3])**2 + (out.perror[10]/out.params[10])**2 + covar_term)  

    hg_ew = np.sqrt(2  * math.pi) *  out2.params[3] * out2.params[10]     ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[3][10] / (out2.params[3] * out2.params[10])
    hg_ew_err =  hg_ew * np.sqrt( (out2.perror[3]/out2.params[3])**2 + (out2.perror[10]/out2.params[10])**2 + covar_term)  

    ### hd 
    hd_flux = np.sqrt(2  * math.pi) *  out.params[2] * out.params[11]   ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[2][11] / (out.params[2] * out.params[11])
    hd_err =  hd_flux * np.sqrt( (out.perror[2]/out.params[2])**2 + (out.perror[11]/out.params[11])**2 + covar_term)

    hd_ew = np.sqrt(2  * math.pi) *  out2.params[2] * out2.params[11]   ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[2][11] / (out2.params[2] * out2.params[11])
    hd_ew_err =  hd_ew * np.sqrt( (out2.perror[2]/out2.params[2])**2 + (out2.perror[11]/out2.params[11])**2 + covar_term)

   
   ### oiii BOTH LINES
    oiii_flux = 1.336 * np.sqrt(2  * math.pi) *  out.params[4] * out.params[12] 
    covar_term = 2  * out.covar[4][12] / (out.params[4] * out.params[12]) 
    oiii_err =  oiii_flux * np.sqrt( (out.perror[4]/out.params[4])**2 + (out.perror[12]/out.params[12])**2 + covar_term) 

    oiii_ew = 1.336 * np.sqrt(2  * math.pi) *  out2.params[4] * out2.params[12] 
    covar_term = 2  * out2.covar[4][12] / (out2.params[4] * out2.params[12]) 
    oiii_ew_err =  oiii_ew * np.sqrt( (out2.perror[4]/out2.params[4])**2 + (out2.perror[12]/out2.params[12])**2 + covar_term) 


    ### oii
    oii_flux = np.sqrt(2  * math.pi) *  out.params[13] * out.params[0]   ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[13][0] / (out.params[13] * out.params[0])
    oii_err =  oii_flux * np.sqrt( (out.perror[13]/out.params[13])**2 + (out.perror[0]/out.params[0])**2 + covar_term) 

    
    oii_ew = np.sqrt(2  * math.pi) *  out2.params[13] * out2.params[0]   ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[13][0] / (out2.params[13] * out2.params[0])
    oii_ew_err =  oii_ew* np.sqrt( (out2.perror[13]/out2.params[13])**2 + (out2.perror[0]/out2.params[0])**2 + covar_term) 


    ### sii
    sii_flux = np.sqrt(2  * math.pi) *  out.params[7] * out.params[14]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[7][14] / (out.params[7] * out.params[14])
    sii_err =  sii_flux * np.sqrt((out.perror[7]/out.params[7])**2 + (out.perror[14]/out.params[14])**2 + covar_term) 

    sii_ew = np.sqrt(2  * math.pi) *  out2.params[7] * out2.params[14]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[7][14] / (out2.params[7] * out2.params[14])
    sii_ew_err =  sii_ew * np.sqrt((out2.perror[7]/out2.params[7])**2 + (out2.perror[14]/out2.params[14])**2 + covar_term) 


    #### he1 5876
    he1_5876_flux = np.sqrt(2  * math.pi) *  out.params[15] * out.params[5]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[15][5] / (out.params[15] * out.params[5])
    he1_5876_err =  he1_5876_flux * np.sqrt((out.perror[15]/out.params[15])**2 + (out.perror[5]/out.params[5])**2 + covar_term)  



    he1_5876_ew = np.sqrt(2  * math.pi) *  out2.params[15] * out2.params[5]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[15][5] / (out2.params[15] * out2.params[5])
    he1_5876_ew_err =  he1_5876_ew * np.sqrt((out2.perror[15]/out2.params[15])**2 + (out2.perror[5]/out2.params[5])**2 + covar_term)  


    ### oi BOTH lines 
    oi_flux = 1.3 * np.sqrt(2  * math.pi) *  out.params[16] * out.params[6]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[16][6] / (out.params[16] * out.params[6])
    oi_err =  oi_flux * np.sqrt((out.perror[16]/out.params[16])**2 + (out.perror[6]/out.params[6])**2 + covar_term)


    oi_ew = 1.3 * np.sqrt(2  * math.pi) *  out2.params[16] * out2.params[6]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[16][6] / (out2.params[16] * out2.params[6])
    oi_ew_err =  oi_ew * np.sqrt((out2.perror[16]/out2.params[16])**2 + (out2.perror[6]/out2.params[6])**2 + covar_term)


    #### Ne III 3870 BLEND
    neiii_flux = np.sqrt(2  * math.pi) *  out.params[1] * out.params[17]     ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[1][17] / (out.params[1] * out.params[17])
    neiii_err =  neiii_flux * np.sqrt((out.perror[1]/out.params[1])**2 + (out.perror[17]/out.params[17])**2 + covar_term)


    neiii_ew = np.sqrt(2  * math.pi) *  out2.params[1] * out2.params[17]     ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out2.covar[1][17] / (out2.params[1] * out2.params[17])
    neiii_ew_err =  neiii_ew * np.sqrt((out2.perror[1]/out2.params[1])**2 + (out2.perror[17]/out2.params[17])**2 + covar_term)
    


    #### gather some metadata on the OIII lines 
    dat = asciitable.read(input_masterlist)  
    #dat2 = asciitable.read(input_masterlist + '.EWlist')

    z = dat['z']  
    foiii = dat['foiii']
    zstack = np.median(z)  

    lum_oiii = foiii * 4 * 3.14159 * cosmo.luminosity_distance(z).cgs.value**2


    median_foiii = np.median(foiii) 
    median_loiii = np.median(lum_oiii) 


    median_fha = median_foiii * ha_flux / oiii_flux 
    median_fhb = median_foiii * hb_flux / oiii_flux 
    median_fhg = median_foiii * hg_flux / oiii_flux 
    median_fhd = median_foiii * hd_flux / oiii_flux 
    median_foii = median_foiii * oii_flux/oiii_flux 
    median_fsii = median_foiii * sii_flux/oiii_flux 
    median_fhe1 = median_foiii * he1_5876_flux / oiii_flux 
    median_foi = median_foiii * oi_flux / oiii_flux 
    median_fneiii = median_foiii * neiii_flux / oiii_flux 

   
    median_lha = median_loiii * ha_flux / oiii_flux 
    median_lhb = median_loiii * hb_flux / oiii_flux 
    median_lhg = median_loiii * hg_flux / oiii_flux 
    median_lhd = median_loiii * hd_flux / oiii_flux 
    median_loii = median_loiii * oii_flux/oiii_flux 
    median_lsii = median_loiii * sii_flux/oiii_flux 
    median_lhe1 = median_loiii * he1_5876_flux / oiii_flux 
    median_loi = median_loiii * oi_flux / oiii_flux 
    median_lneiii = median_loiii * neiii_flux / oiii_flux 



    #w = np.where((lam > 4800) & (lam < 5200)) 
    #median_cont_oiii = np.median(median_cont[w]) 
    #err_cont_oiii = np.median(err_cont[w]) 
    #median_ew_oiii = median_foiii / median_cont_oiii / (1+zstack) 

   # w = np.where((lam > 6400) & (lam < 7000)) 
   # if np.size(w) > 0: 
   #     median_cont_ha = np.median(median_cont[w]) 
   #     err_cont_ha = np.median(err_cont[w]) 
   #     median_ew_ha = median_fha / median_cont_ha  / (1+zstack)

   # w = np.where((lam > 4600) & (lam < 5000)) 
   # median_cont_hb = np.median(median_cont[w]) 
   # err_cont_hb = np.median(err_cont[w]) 
   # median_ew_hb = median_fhb / median_cont_hb / (1+zstack) 


   # w = np.where((lam > 4100) & (lam < 4500)) 
   # median_cont_hg = np.median(median_cont[w]) 
   # err_cont_hg = np.median(err_cont[w]) 
   # median_ew_hg = median_fhg / median_cont_hg  / (1+zstack) 

   # w = np.where((lam > 3800) & (lam < 4200)) 
   # median_cont_hd = np.median(median_cont[w]) 
   # err_cont_hd = np.median(err_cont[w]) 
   # median_ew_hd = median_fhd / median_cont_hd  / (1+ zstack) 

   # w = np.where((lam > 3500) & (lam < 3900)) 
   # median_cont_oii = np.median(median_cont[w]) 
   # err_cont_oii = np.median(err_cont[w]) 
   # median_ew_oii = median_foii / median_cont_oii / (1+zstack) 

   # w = np.where((lam > 6400) & (lam < 7000)) 
   # if np.size(w) > 0: 
   #     median_cont_sii = np.median(median_cont[w]) 
   #     err_cont_sii = np.median(err_cont[w]) 
   #     median_ew_sii = median_fsii / median_cont_sii / (1+zstack) 

    #w = np.where((lam > 5700) & (lam < 6100)) 
    #if np.size(w) > 0 : 
    #    median_cont_he1 = np.median(median_cont[w]) 
    #    err_cont_he1 = np.median(err_cont[w]) 
    #    median_ew_he1 = median_fhe1 / median_cont_he1 / (1+zstack) 

    #w = np.where((lam > 6100) & (lam < 6500))
    #if np.size(w) > 0: 
    #    median_cont_oi = np.median(median_cont[w]) 
    #    err_cont_oi = np.median(err_cont[w]) 
    #    median_ew_oi = median_foi / median_cont_oi / (1+zstack) 

    #w = np.where((lam > 3600) & (lam < 4000)) 
    #median_cont_neiii = np.median(median_cont[w]) 
    #err_cont_neiii = np.median(err_cont[w]) 
    #median_ew_neiii = median_fneiii / median_cont_neiii / (1+zstack) 




    if lam_max < 6564.: 
        ha_flux = -1
        ha_err = -1 
        median_fha = -1 
        median_lha = -1
        #median_cont_ha = -1 
        #err_cont_ha = -1 
        ha_ew = -1 
        ha_ew_err = -1 
    if lam_max < 6725 : 
        sii_flux = -1
        sii_err = -1 
        median_fsii = -1 
        median_lsii = -1 
        #median_cont_sii = -1 
        #err_cont_sii = -1 
        sii_ew = -1 
        sii_ew_err = -1 
    if lam_max < 5876 : 
        he1_5876_flux = -1
        he1_5876_err = -1
        median_fhe1 = -1 
        median_lhe1 = -1 
        #median_cont_he1 = -1 
        #err_cont_he1 = -1 
        he1_5876_ew = -1
        he1_5876_ew_err = -1 
    if lam_max <  6300: 
        oi_flux = -1 
        oi_err = -1 
        median_foi = -1 
        median_loi = -1 
        #median_cont_oi = -1 
        #err_cont_oi = -1 
        oi_ew = -1 
        oi_ew_err = -1






    label  = np.array(['HaNII', 'Hb', 'HgOIII', 'Hd', 'OIIIboth', 'OII', 'SII', 'HeI5876', 'OIboth', 'NeIIIblend']) 
    norm_fluxes = np.array([ha_flux, hb_flux, hg_flux, hd_flux, oiii_flux, oii_flux, sii_flux, he1_5876_flux, oi_flux, neiii_flux]) 
    norm_errs =  np.array([ha_err,  hb_err,  hg_err, hd_err, oiii_err, oii_err, sii_err, he1_5876_err, oi_err, neiii_err]) 
    scaled_fluxes = np.array([median_fha, median_fhb, median_fhg, median_fhd, median_foiii, median_foii, median_fsii, median_fhe1, median_foi, median_fneiii]) 
    scaled_lums   = np.array([median_lha, median_lhb, median_lhg, median_lhd, median_loiii, median_loii, median_lsii, median_lhe1, median_loi, median_lneiii]) 
    #cont_est = np.array([median_cont_ha, median_cont_hb, median_cont_hg, median_cont_hd, median_cont_oiii, median_cont_oii, median_cont_sii, median_cont_he1, median_cont_oi, median_cont_neiii]) 
    #cont_est_err = np.array([err_cont_ha, err_cont_hb, err_cont_hg, err_cont_hd, err_cont_oiii, err_cont_oii, err_cont_sii, err_cont_he1, err_cont_oi, err_cont_neiii])
    ew_stack = np.array([ha_ew, hb_ew, hg_ew, hd_ew, oiii_ew, oii_ew, sii_ew, he1_5876_ew, oi_ew, neiii_ew]) 
    ew_stack_err = np.array([ha_ew_err,hb_ew_err, hg_ew_err, hd_ew_err, oiii_ew_err, oii_ew_err, sii_ew_err, he1_5876_ew_err, oi_ew_err, neiii_ew_err]) 


    out_names = ['line', 'flux_norm', 'flux_norm_err', 'flux_scale', 'luminosity_scale', 'ew_rest', 'ew_err'] 
    out_data = [label, norm_fluxes, norm_errs, scaled_fluxes, scaled_lums,  ew_stack, ew_stack_err]  

    asciitable.write(out_data, output_meas, names= out_names, overwrite=True)

    tab = asciitable.read(output_meas)
    tab['flux_norm'].format = '7.4f' 
    tab['flux_norm_err'].format = '7.4f' 
    tab['flux_scale'].format = '7.4e'
    tab['luminosity_scale'].format = '7.4e'
    tab['ew_rest'].format = '7.3f' 
    tab['ew_err'].format = '7.3f' 
 

    asciitable.write(tab, output_meas, overwrite=True, format = 'fixed_width') 
    




    #output = open(output_meas, 'w') 
    #output.write('# Line               flux_norm       flux_norm_err       flux_scale    luminosity_scale     median_continuum    err_continuum    ew_est \n') 
    #output.write('Ha+ NII            ' + str(ha_flux) + '  ' + str(ha_err) +  '  ' + str(median_fha) + '  ' + str(median_lha) +  '  '  '\n')   ### the ha_flux corresponds only to Ha + NII because the nii amp in the model is fixed at zero.
    #output.write('Hb                 ' + str(hb_flux) + '  ' + str(hb_err) +  '  ' + str(median_fhb) + '  ' + str(median_lha) + '\n') 
    #output.write('Hg + OIII          ' + str(hg_flux) + '  ' + str(hg_err) +  '  ' + str(median_fhg) + '  ' + str(median_lhg) +  '\n')  
    #output.write('Hd                 ' + str(hd_flux) + '  ' + str(hd_err) +  '  ' + str(median_fhd) + '  ' + str(median_lhd) + '\n') 
    #output.write('OIII (both)        ' + str(oiii_flux) + '   ' +str(oiii_err) + '  ' + str(median_foiii) + '  ' + str(median_loiii) +  '  ' + str(median_ew_oiii) +  '\n') 
    #output.write('OII                ' + str(oii_flux) + '   ' +str(oii_err)  + '  ' + str(median_foii) + '  ' + str(median_loii) +'\n') 
    #output.write('SII (both)         ' + str(sii_flux) + '   '  + str(sii_err) +  '  ' + str(median_fsii) + '  ' + str(median_lsii) + '\n')  
    #output.write('HeI 5876           ' + str(he1_5876_flux) + '   ' + str(he1_5876_err) + '  ' + str(median_fhe1) +   + '  ' + str(median_lhe1) + '\n') 
    #output.write('OI 6300+6363       ' + str(oi_flux) + '  ' + str(oi_err) + '  ' + str(median_foi) + '  ' + str(median_loi) + '\n') 
    #output.write('Ne III 3869 blend  ' + str(neiii_flux) + '  ' + str(neiii_err) +  '  ' + str(median_fneiii) + '  ' + str(median_lneiii) + '\n') 

    #output.close() 


    #plt.axvline(x = 5007 - 100) 
    #plt.axvline(x = 5007 + 100)

    


 
    f, axarr = plt.subplots(2, 1,  figsize=(8, 10))

    axarr[0].plot(lam, flux, 'k-', ls='steps-mid') 
    axarr[0].plot(lam, model_fit, 'r--')
    #plt.axhline(y=0, color = 'b', ls = '--')
     

    axarr[1].plot(lam, median_cont, 'k-', ls='steps-mid') 
    axarr[1].plot(lam, model_fit2, 'r--') 

    plt.savefig(output_fig) 
    if showfig == True :
        plt.show()  






    
