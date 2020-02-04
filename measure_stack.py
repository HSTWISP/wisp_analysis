from wisp_analysis import *
from mpfit import *

import astropy.units as u
from astropy.cosmology import Planck13 as cosmo


def stack_emline_model(pars, x, comps = False):
    sigma_oiiihb = pars[4] 
    
    sigma_oii = pars[0] * sigma_oiiihb
    sigma_neiii = pars[1]  * sigma_oiiihb
    sigma_hd = pars[2]  * sigma_oiiihb
    sigma_hg = pars[3] * sigma_oiiihb
    sigma_he1_5876 = pars[5]   * sigma_oiiihb
    sigma_oi = pars[6]  * sigma_oiiihb
    sigma_hasii = pars[7]  * sigma_oiiihb 
    
    sigma_heii = pars[47] * sigma_oiiihb 

    unity = pars[46]  ### need this parameter for propagating errors. fix it to 1. 



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

    heii_amp = pars[48] 


    hasii_shift = pars[20] 
    oi_shift = pars[21] 
    he1_5876_shift = pars[22] 
    oiii_shift = pars[23] 
    hb_shift = pars[24]
    hg_shift = pars[25] 
    hd_shift = pars[26] 
    neiii_shift = pars[27] 
    oii_shift = pars[28]
    heii_shift = pars[49] 

    c_3800_4200 = pars[29] 
    c_5200_6400 = pars[30]
    s_5200_6500  = pars[31]
    cnorm = pars[32]
    c_4500_5100 = pars[51] 




    ha_amp_broad = pars[33] 
    hb_amp_broad = pars[34] 
    hg_amp_broad = pars[35] 
    hd_amp_broad = pars[36]
    oiii_5007_amp_broad = pars[37] 
    oii_amp_broad = pars[38] 
    sii_amp_broad = pars[39] 
    he1_5876_amp_broad = pars[40] 
    oi_amp_broad = pars[41] 
    neiii_amp_broad = pars[42]
    nii_6583_amp_broad = pars[43] 
    he1_6678_amp_broad = pars[44]  

    heii_amp_broad = pars[50] 


    sigma_oiiihb_broad = pars[45] * sigma_oiiihb
    sigma_oii_broad = pars[0] * sigma_oiiihb_broad
    sigma_neiii_broad = pars[1]  * sigma_oiiihb_broad
    sigma_hd_broad = pars[2]  * sigma_oiiihb_broad
    sigma_hg_broad = pars[3] * sigma_oiiihb_broad 
    sigma_he1_5876_broad = pars[5]   * sigma_oiiihb_broad
    sigma_oi_broad = pars[6]  * sigma_oiiihb_broad
    sigma_hasii_broad = pars[7]  * sigma_oiiihb_broad 
    sigma_heii_broad = pars[47] * sigma_oiiihb_broad 



    cont = np.zeros(np.size(x)) 
    w=np.where( (x > 3800) & (x < 4400)) 
    cont[w] = c_3800_4200 
    w=np.where( (x> 6000) & (x < 6450)) 
    cont[w] = c_5200_6400 + x[w] * s_5200_6500 
    w=np.where( (x > 4500) & (x < 5500)) 
    cont[w] = c_4500_5100 

    cont = cont + cnorm

  


    model_narrow  = ha_amp * gaussian(x, 6564.6 + hasii_shift, sigma_hasii) + \
            hb_amp  * gaussian(x, 4862.7 + hb_shift, sigma_oiiihb) +   \
            hg_amp * gaussian(x, 4341.7 + hg_shift, sigma_hg) +  \
            hd_amp * gaussian(x, 4102.9 + hd_shift, sigma_hd) +   \
            oiii_5007_amp * gaussian(x, 5008. + oiii_shift, sigma_oiiihb) +   \
            oiii_5007_amp/3 * gaussian(x, 4960. + oiii_shift, sigma_oiiihb) +  \
            oii_amp * gaussian(x, 3728. + oii_shift, sigma_oii) + \
            sii_amp * gaussian(x, 6725.  + hasii_shift, sigma_hasii) +  \
            he1_5876_amp  * gaussian(x, 5877.2 + he1_5876_shift, sigma_he1_5876) +   \
            oi_amp *  gaussian(x, 6302. + oi_shift, sigma_oi) +  \
            oi_amp /3 * gaussian(x, 6365.5 + oi_shift,sigma_oi) + \
            neiii_amp * gaussian(x, 3870 + neiii_shift, sigma_neiii) +  \
            nii_6583_amp * gaussian(x, 6585.23 + hasii_shift, sigma_hasii) +  \
            nii_6583_amp / 3 * gaussian(x, 6550. + hasii_shift, sigma_hasii)+\
            he1_6678_amp  * gaussian(x, 6680. + hasii_shift, sigma_hasii) +\
            heii_amp * gaussian(x, 4686. + heii_shift, sigma_heii) 
    
    model_broad =   ha_amp_broad * gaussian(x, 6564.6 + hasii_shift, sigma_hasii_broad) + \
            hb_amp_broad * gaussian(x, 4862.7 + hb_shift, sigma_oiiihb_broad)  +  \
            hg_amp_broad* gaussian(x, 4341.7 + hg_shift, sigma_hg_broad) +  \
            hd_amp_broad * gaussian(x, 4102.9 + hd_shift, sigma_hd_broad) +  \
            oiii_5007_amp_broad * gaussian(x, 5008. + oiii_shift, sigma_oiiihb_broad) + \
            oiii_5007_amp_broad / 3. * gaussian(x, 5008. + oiii_shift, sigma_oiiihb_broad) + \
            oii_amp_broad * gaussian(x, 3728. + oii_shift, sigma_oii_broad) + \
            sii_amp_broad * gaussian(x, 6725 + hasii_shift, sigma_hasii_broad) + \
            he1_5876_amp_broad  * gaussian(x, 5877.2 + he1_5876_shift, sigma_he1_5876_broad) +  \
            oi_amp_broad * gaussian(x, 6302 + oi_shift, sigma_oi_broad)  +\
            oi_amp_broad/3 * gaussian(x, 6302 + oi_shift, sigma_oi_broad) +\
            neiii_amp_broad * gaussian(x, 3870 + neiii_shift, sigma_neiii_broad) + \
            nii_6583_amp_broad * gaussian(x, 6585.23 + hasii_shift, sigma_hasii_broad) +  \
            nii_6583_amp_broad / 3 * gaussian(x, 6550. + hasii_shift, sigma_hasii_broad) + \
            he1_6678_amp_broad  * gaussian(x, 6680. + hasii_shift, sigma_hasii_broad)  +\
            heii_amp_broad * gaussian(x, 4686.  + heii_shift, sigma_heii_broad) 

            
#### 52 parameters, 0-51
            

    model = model_narrow + model_broad + cont        
        

    if comps == True : 
        return [model, model_narrow, model_broad, cont]
    else :
        return model


    

def model_resid_stack(pars, fjac=None, lam = None, flux = None, err = None):
    model = stack_emline_model(pars, lam) 
    resid = (flux- model) / err 
   
    status = 0 
    return [status, resid]



def error_2comp_gaussian(params, perror, covar, indices): 

    A1 = params[indices[0]] 
    A1err  = perror[indices[0]] 
    A2 = params[indices[1]] 
    A2err = perror[indices[1]] 
    sigma = params[indices[2]] 
    sigmaerr = perror[indices[2]] 
    S1 = params[indices[3]] 
    S1err = perror[indices[3]] 
    S2 = params[indices[4]] 
    S2err = perror[indices[4]]

    c_A1A2 =    covar[indices[0]][indices[1]]
    c_A1sigma = covar[indices[0]][indices[2]] 
    c_A1S1 =    covar[indices[0]][indices[3]] 
    c_A1S2 =    covar[indices[0]][indices[4]]

    c_A2sigma = covar[indices[1]][indices[2]] 
    c_A2S1 =    covar[indices[1]][indices[3]] 
    c_A2S2 =    covar[indices[1]][indices[4]] 

    c_sigmaS1 = covar[indices[2]][indices[3]] 
    c_sigmaS2 = covar[indices[2]][indices[4]] 
    
    c_S1S2 = covar[indices[3]][indices[4]] 

    
    C = np.sqrt(2 *  math.pi) 
    F = C  * (A1 * sigma * S1 + A2 * sigma * S1 * S2) 

    ### see hand-written notes, 06/17/19, on how to propagate errors for this, including the covariances.  
    term1 = ( C * sigma * S1 * A1err)**2 
    term2 = ( C * sigma * S1 * S2 * A2err)**2 
    term3 = ( C * (A1 * S1 + A2 * S1 * S2) * sigmaerr)**2 
    term4 = ( C * (A1 * sigma + A2 * sigma * S2) * S1err)**2 
    term5 = ( C * A2 * sigma  * S1 * S2err)**2  

    ### now the covariances. 
    term6 = 2 * (C * sigma * S1) * (C * sigma  * S1 * S2) * c_A1A2 
    term7 = 2 * (C * sigma * S1) * (C * (A1*S1 + A2*S1*S2)) * c_A1sigma 
    term8 = 2 * (C * sigma * S1) * (C * (A1 * sigma + A2 * sigma * S2)) * c_A1S1 
    term9 = 2 * (C * sigma * S1) * (C * A2 * sigma * S1) * c_A1S2 

    term10 = 2 * (C * sigma * S1 * S2) * (C  * (A1 * S1 + A2 * S1 * S2)) * c_A2sigma 
    term11 =  2 * (C * sigma * S1 * S2) * (C * (A1 * sigma + A2 * sigma * S2)) * c_A2S1 
    term12 = 2 * (C * sigma * S1 * S2) * (C * A2 * sigma * S1) * c_A2S2 

    term13 = 2 * (C * (A1 * S1 + A2 * S1 * S2)) * (C  * (A1 * sigma + A2 * sigma * S2)) * c_sigmaS1 
    term14=  2 * (C * (A1 * S1+ A2 * S1 * S2)) * (C * A2 * sigma * S1) * c_sigmaS2 

    term15 = 2 * (C * (A1 * sigma + A2 * sigma * S2)) * (C * A2 * sigma * S1) * c_S1S2 

    variance = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + term11 + term12 + term13 + term14 + term15 

    error = np.sqrt(variance)


    return [F, error] 


    



def measure_stack(input_stack, input_masterlist, output_meas, output_fig, zmax = 2.3, showfig = False, twocomp=True): 


    lam_max = 17200/(1+zmax) 

    tab= asciitable.read(input_stack) 
    lam = tab['lam']
    flux = tab['flux_norm']
    err = tab['err']
    #median_cont = tab['median_cont'] 
    #err_cont = tab['err_cont'] 

    w=np.where(lam < lam_max) 
    w=w[0] 
    lam = lam[w]
    flux = flux[w] 
    err = err[w] 
    #median_cont = median_cont[w] 
    #err_cont = err_cont[w] 


    w=np.where(err > 0) 
    w=w[0] 
    lam = lam[w]
    flux = flux[w] 
    err = err[w]
   # median_cont = median_cont[w] 
   # err_cont = err_cont[w]


    #w=np.where(err_cont > 0) 
    #w=w[0] 
    #lam = lam[w]
    #flux = flux[w] 
    #err = err[w]
   # median_cont = median_cont[w] 
   # err_cont = err_cont[w]


  
    pguess= np.zeros(52) 
    pguess[0] = 1  ### line widths are specified relative to oiii+hb line widths now. 
    pguess[1] = 1 
    pguess[2] = 1
    pguess[3] = 1 
    pguess[4] = 10
    pguess[5] = 1 
    pguess[6] = 1 
    pguess[7] = 1 
    pguess[47] = 1   ### heii added 1/28/20

    ### fwhm of broad component relative to the narrow component. 
    pguess[45]  = 2.0 

    #### amplitudes
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
    pguess[19] = 0.000  ### he1 6678
    pguess[48] = 0.0003   ### he ii added 1/28/20 
   
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
    pguess[49] = 0.0    ### he ii added 1/28/20 

    pguess[29] = 0.00001
    pguess[30]  = -0.01
    pguess[31] = 0.0000 
    pguess[32] = 0.0 
    pguess[51] = 0.0
    


    ### broad amplitude guesses
    if twocomp == True :

        pguess[33] = 0.0002
        pguess[34] = 0.0002
        pguess[35]  = 0.0002
        pguess[36] = 0.0002
        pguess[37] = 0.0002
        pguess[38] = 0.0002
        pguess[39] = 0.0002
        pguess[40] = 0.0002
        pguess[41] = 0.0002
        pguess[42] = 0.0002
        pguess[43] = 0.  ###nii 
        pguess[44] = 0.0000 ### he 1 6678
        pguess[50] = 0.0002   ###  he ii added 1/28/20 


    else : 
        pguess[33] = 0.0000
        pguess[34] = 0.0000
        pguess[35]  = 0.0000
        pguess[36] = 0.0000
        pguess[37] = 0.0000
        pguess[38] = 0.0000
        pguess[39] = 0.0000
        pguess[40] = 0.0000
        pguess[41] = 0.0000
        pguess[42] = 0.0000
        pguess[43] = 0.  ###nii 
        pguess[44] = 0.0000 ### he 1 6678
        pguess[50]= 0.000   #### he ii added 1/28/20 



    pguess[46] = 1.0   ### unity, a dummy parameter


    npars = len(pguess) 
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} 
              for i in range(npars)]
    for i in range(npars): parinfo[i]['value'] = pguess[i]  
   
 
    ### set the nii and he1 6678 =0 and fixed it. always. 
    parinfo[18]['fixed'] = 1
    parinfo[19]['fixed'] = 1
    parinfo[43]['fixed'] = 1
    parinfo[44]['fixed'] = 1
   
    ### necessary dummy parameter  
    parinfo[46]['fixed'] = 1


    #### fwhm are specified relative to oiiihb 
    parinfo[0]['limited'][0] = 1 
    parinfo[0]['limited'][1] = 1 
    parinfo[0]['limits'][0] = 0.5
    parinfo[0]['limits'][1] = 2.0

    parinfo[1]['limited'][0] = 1
    parinfo[1]['limited'][1] = 1
    parinfo[1]['limits'][0] = 0.5
    parinfo[1]['limits'][1] = 2.0

    parinfo[2]['limited'][0] = 1
    parinfo[2]['limited'][1] = 1
    parinfo[2]['limits'][0] = 0.5
    parinfo[2]['limits'][1] = 2.0

    parinfo[3]['limited'][0] = 1
    parinfo[3]['limited'][1] = 1
    parinfo[3]['limits'][0] = 0.5
    parinfo[3]['limits'][1] = 2.0

    #### oiii + hb, sets the fwhm 
    parinfo[4]['limited'][0] = 1
    parinfo[4]['limited'][1] = 1
    parinfo[4]['limits'][0] = 5
    parinfo[4]['limits'][1] = 50 

    parinfo[5]['limited'][0] = 1
    parinfo[5]['limited'][1] = 1
    parinfo[5]['limits'][0] = 0.5
    parinfo[5]['limits'][1] = 2.0

    parinfo[6]['limited'][0] = 1
    parinfo[6]['limited'][1] = 1
    parinfo[6]['limits'][0] = 0.5
    parinfo[6]['limits'][1] = 2.0


    parinfo[7]['limited'][0] = 1
    parinfo[7]['limited'][1] = 1
    parinfo[7]['limits'][0] = 0.5
    parinfo[7]['limits'][1] = 2.0

    parinfo[47]['limited'][0] = 1   ### he ii 
    parinfo[47]['limited'][1] = 1
    parinfo[47]['limits'][0] = 0.5
    parinfo[47]['limits'][1] = 4.0




    ### amplitude of broad, relative to amplitude of narrow: 
    parinfo[45]['limited'][0] = 1
    parinfo[45]['limited'][1] = 1
    parinfo[45]['limits'][0] = 1 
    parinfo[45]['limits'][1] = 4. 
   
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
   # parinfo[19]['limited'][0] = 1
   # parinfo[19]['limits'][0] = 0
    parinfo[48]['limited'][0] = 1    ### he ii  
    parinfo[48]['limits'][0] = 0



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
    parinfo[49]['limited']  = [1, 1]    ### he ii 
    parinfo[49]['limits'] = [-10, 10]

   


    ### residual continuum. 
    parinfo[29]['limited']  = [1, 1] 
    parinfo[29]['limits'] = [-0.01, 0.01]
    parinfo[30]['limited']  = [1, 1] 
    parinfo[30]['limits'] = [-0.01, 0.01 ]
    parinfo[30]['limited']  = [1, 1] 
    parinfo[30]['limits'] = [-0.01, 0.01 ]  
    parinfo[51]['limited'] = [1, 1] 
    parinfo[51]['limits'] = [-0.01, 0.01] 


    parinfo[31]['fixed'] = 1  ## fix the slope of the continuum around he1/oi to zero; it is fitting in a wacky way. 
    parinfo[32]['fixed'] = 1  ### fix continuum normalization to zero for continuum subtracted spectra. 


    if twocomp == True: 
    #### broad component amplitudes must be positive. 
       parinfo[33]['limited'][0] = 1
       parinfo[33]['limits'][0] = 0
       parinfo[34]['limited'][0] = 1
       parinfo[34]['limits'][0] = 0
       parinfo[35]['limited'][0] = 1
       parinfo[35]['limits'][0] = 0
       parinfo[36]['limited'][0] = 1
       parinfo[36]['limits'][0] = 0
       parinfo[37]['limited'][0] = 1
       parinfo[37]['limits'][0] = 0
       parinfo[38]['limited'][0] = 1
       parinfo[38]['limits'][0] = 0
       parinfo[39]['limited'][0] = 1
       parinfo[39]['limits'][0] = 0
       parinfo[40]['limited'][0] = 1
       parinfo[40]['limits'][0] = 0
       parinfo[41]['limited'][0] = 1
       parinfo[41]['limits'][0] = 0
       parinfo[42]['limited'][0] = 1
       parinfo[42]['limits'][0] = 0
       parinfo[50]['limited'][0] = 1   ### he ii 
       parinfo[50]['limits'][0] = 0



      # parinfo[44]['limited'][0] = 1
      # parinfo[44]['limits'][0] = 0


    else :
        parinfo[33]['fixed'] = 1
        parinfo[34]['fixed'] = 1
        parinfo[35]['fixed'] = 1
        parinfo[36]['fixed'] = 1
        parinfo[37]['fixed'] = 1
        parinfo[38]['fixed'] = 1
        parinfo[39]['fixed'] = 1
        parinfo[40]['fixed'] = 1
        parinfo[41]['fixed'] = 1
        parinfo[42]['fixed'] = 1   
        parinfo[50]['fixed'] = 1 

    

    #### do fit and evaluate model 
    fa = {'lam':lam, 'flux':flux, 'err':err} 
    out = mpfit(model_resid_stack, pguess, functkw=fa, parinfo = parinfo, quiet=True) 
    model_fit = stack_emline_model(out.params, lam)


    #### write some other output meta_data 
    outfits = stack_emline_model(out.params, lam, comps = True)   
    outfits_data = [lam, outfits[0], outfits[1], outfits[2], outfits[3]] 
    outfits_names = ['lam', 'total', 'narrow', 'broad', 'cont'] 
    asciitable.write(outfits_data, output_meas+ '.fitprof', names = outfits_names, overwrite=True)

    ### parameters
    asciitable.write(out.params, output_meas+ '.fitpars', overwrite=True) 





    #### redo model fit on continuum normalized spectrum
    ### update the pguess based on what the previous out.params gave. 

   # pguess[0:8] = out.params[0:8] 
   # pguess[8:20]  = out.params[8:20] * 300.   ### scale since the spectra are not line flux normalized. 
   # pguess[20:29] = out.params[20:29]
    #pguess[10]  = pguess[9]
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
 

    #pguess[32] = 1.0 
    #parinfo[32]['limited'] = [1, 1] 
    #parinfo[32]['limits'] = [0.7, 1.5]


    #fa2 = {'lam':lam, 'flux':median_cont, 'err':err_cont} 
    #out2 = mpfit(model_resid_stack, pguess, functkw=fa2, parinfo = parinfo, quiet=True) 
    #model_fit2 = stack_emline_model(out2.params, lam)
     

    




    ### gather line fluxes 
    ### ha  
    ### note that nii  and he1 are set to zero.  This means that this is Ha + NII
    ### and I'm not sure if He 1 is included in the flux or resolved, or something in between

    #### ha 
    indices = [8, 33, 4, 7, 45]
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices) 
    ha_flux = measurements[0]
    ha_err = measurements[1] 

    #ha_flux = np.sqrt(2  * math.pi) *  out.params[7] * out.params[8]  * out.params[4]  ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[7][8] / (out.params[7] * out.params[8]) + 2  * out.covar[4][8] / (out.params[4] * out.params[8]) +  2  * out.covar[7][4] / (out.params[7] * out.params[4]) 
    #ha_err =  ha_flux * np.sqrt( (out.perror[7]/out.params[7])**2 + (out.perror[8]/out.params[8])**2   + (out.perror[4]/out.params[4])**2 +  covar_term) 

    #print ha_flux, ha_err
     
    ### the "flux" in the continuum normalized model is really the EW, since Flam_cont = 1. 
   # ha_ew = np.sqrt(2  * math.pi) *  out2.params[7] * out2.params[8]  * out2.params[4]  ## sqrt(2 * pi) * amplitude * sigma 
   # covar_term = 2  * out2.covar[7][8] / (out2.params[7] * out2.params[8]) + 2  * out2.covar[4][8] / (out2.params[4] * out2.params[8]) +  2  * out2.covar[7][4] / (out2.params[7] * out2.params[4]) 
   # ha_ew_err =  ha_ew * np.sqrt( (out2.perror[7]/out2.params[7])**2 + (out2.perror[8]/out2.params[8])**2 + (out2.perror[4]/out2.params[4])**2   +  covar_term) 

    ##hb
    indices = [9, 34, 4, 46, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    hb_flux = measurements[0] 
    hb_err = measurements[1] 
   

    #hb_flux = np.sqrt(2  * math.pi) *  out.params[4] * out.params[9]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[4][9] / (out.params[4] * out.params[9])
    #hb_err =  hb_flux * np.sqrt( (out.perror[4]/out.params[4])**2 + (out.perror[9]/out.params[9])**2 + covar_term)  

   # hb_ew = np.sqrt(2  * math.pi) *  out2.params[4] * out2.params[9] ## sqrt(2 * pi) * amplitude * sigma 
   # covar_term = 2  * out2.covar[4][9] / (out2.params[4] * out2.params[9])
    #hb_ew_err =  hb_ew * np.sqrt( (out2.perror[4]/out2.params[4])**2 + (out2.perror[9]/out2.params[9])**2 + covar_term)

    ##hg 
    indices = [10, 35, 4, 3, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    hg_flux = measurements[0] 
    hg_err = measurements[1] 

    #hg_flux = np.sqrt(2  * math.pi) *  out.params[3] * out.params[10]  * out.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[3][10] / (out.params[3] * out.params[10])  +  2  * out.covar[3][4] / (out.params[3] * out.params[4]) +  2  * out.covar[4][10] / (out.params[4] * out.params[10])
    #hg_err =  hg_flux * np.sqrt( (out.perror[3]/out.params[3])**2 + (out.perror[10]/out.params[10])**2  +  (out.perror[4]/out.params[4])**2  + covar_term)   

    #hg_ew = np.sqrt(2  * math.pi) *  out2.params[3] * out2.params[10] * out2.params[4]     ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out2.covar[3][10] / (out2.params[3] * out2.params[10]) + 2  * out2.covar[4][10] / (out2.params[4] * out2.params[10]) + 2  * out2.covar[3][4] / (out2.params[3] * out2.params[4])
    #hg_ew_err =  hg_ew * np.sqrt( (out2.perror[3]/out2.params[3])**2 + (out2.perror[10]/out2.params[10])**2 + (out2.perror[4]/out2.params[4])**2 + covar_term)  

    ### hd 
    indices = [11, 36, 4, 2, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    hd_flux = measurements[0] 
    hd_err = measurements[1] 
 
    #hd_flux = np.sqrt(2  * math.pi) *  out.params[2] * out.params[11] * out.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[2][11] / (out.params[2] * out.params[11]) + 2  * out.covar[4][11] / (out.params[4] * out.params[11]) + 2  * out.covar[2][4] / (out.params[2] * out.params[4])
    #hd_err =  hd_flux * np.sqrt( (out.perror[2]/out.params[2])**2 + (out.perror[11]/out.params[11])**2  +  (out.perror[4]/out.params[4])**2  + covar_term)

    #hd_ew = np.sqrt(2  * math.pi) *  out2.params[2] * out2.params[11] * out2.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out2.covar[2][11] / (out2.params[2] * out2.params[11]) +   2  * out2.covar[2][4] / (out2.params[2] * out2.params[4]) +  2  * out2.covar[4][11] / (out2.params[4] * out2.params[11])
    #hd_ew_err =  hd_ew * np.sqrt( (out2.perror[2]/out2.params[2])**2 + (out2.perror[11]/out2.params[11])**2 + (out2.perror[4]/out2.params[4])**2 + covar_term)

   ### oiii BOTH LINES
    
    indices = [12, 37, 4, 46, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    oiii_flux = 1.3 *  measurements[0] 
    oiii_err = 1.3 *  measurements[1] 
     
    #oiii_flux = 1.3 * np.sqrt(2  * math.pi) *  out.params[4] * out.params[12] 
    #covar_term = 2  * out.covar[4][12] / (out.params[4] * out.params[12] ) 
    #oiii_err =  oiii_flux * np.sqrt( (out.perror[4]/out.params[4])**2 + (out.perror[12]/out.params[12])**2 + covar_term) 

    #oiii_ew = 1.3 * np.sqrt(2  * math.pi) *  out2.params[4] * out2.params[12] 
    #covar_term = 2  * out2.covar[4][12] / (out2.params[4] * out2.params[12]) 
    #oiii_ew_err =  oiii_ew * np.sqrt( (out2.perror[4]/out2.params[4])**2 + (out2.perror[12]/out2.params[12])**2 + covar_term) 

    ### oii

    indices=  [13, 38,4, 0, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    oii_flux =    measurements[0] 
    oii_err =   measurements[1] 

    #oii_flux = np.sqrt(2  * math.pi) *  out.params[13] * out.params[0] * out.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[13][0] / (out.params[13] * out.params[0]) +  2  * out.covar[13][4] / (out.params[13] * out.params[4]) +  2  * out.covar[4][0] / (out.params[4] * out.params[0])
    #oii_err =  oii_flux * np.sqrt( (out.perror[13]/out.params[13])**2 + (out.perror[0]/out.params[0])**2 + (out.perror[4]/out.params[4])**2  +  covar_term) 

    
   # oii_ew = np.sqrt(2  * math.pi) *  out2.params[13] * out2.params[0] * out2.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
   # covar_term = 2  * out2.covar[13][0] / (out2.params[13] * out2.params[0]) + 2  * out2.covar[4][0] / (out2.params[4] * out2.params[0])  + 2  * out2.covar[13][4] / (out2.params[13] * out2.params[4])
   # oii_ew_err =  oii_ew* np.sqrt( (out2.perror[13]/out2.params[13])**2 + (out2.perror[0]/out2.params[0])**2 + (out2.perror[4]/out2.params[4])**2 + covar_term) 


    ### sii
    indices = [14, 39, 4, 7, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    sii_flux = measurements[0]
    sii_err = measurements[1] 

    #sii_flux = np.sqrt(2  * math.pi) *  out.params[7] * out.params[14]  * out.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[7][14] / (out.params[7] * out.params[14]) +  2  * out.covar[7][4] / (out.params[7] * out.params[4]) +  2  * out.covar[4][14] / (out.params[4] * out.params[14])
    #sii_err =  sii_flux * np.sqrt((out.perror[7]/out.params[7])**2 + (out.perror[14]/out.params[14])**2 + (out.perror[4]/out.params[4])**2  + covar_term) 

    #sii_ew = np.sqrt(2  * math.pi) *  out2.params[7] * out2.params[14] * out2.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out2.covar[7][14] / (out2.params[7] * out2.params[14]) +  2  * out2.covar[7][4] / (out2.params[7] * out2.params[4]) +  2  * out2.covar[4][14] / (out2.params[4] * out2.params[14])
    #sii_ew_err =  sii_ew * np.sqrt((out2.perror[7]/out2.params[7])**2 + (out2.perror[14]/out2.params[14])**2 +  (out2.perror[4]/out2.params[4])**2  + covar_term) 


    #### he1 5876
    indices = [14, 40, 4, 5, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    he1_5876_flux = measurements[0] 
    he1_5876_err = measurements[1] 

    #he1_5876_flux = np.sqrt(2  * math.pi) *  out.params[15] * out.params[5] * out.params[4]    ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[15][5] / (out.params[15] * out.params[5]) + 2  * out.covar[4][5] / (out.params[4] * out.params[5])  + 2  * out.covar[15][4] / (out.params[15] * out.params[4]) 
    #he1_5876_err =  he1_5876_flux * np.sqrt((out.perror[15]/out.params[15])**2 + (out.perror[5]/out.params[5])**2 +  (out.perror[4]/out.params[4])**2  + covar_term)  

    #he1_5876_ew = np.sqrt(2  * math.pi) *  out2.params[15] * out2.params[5] * out2.params[4]   ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out2.covar[15][5] / (out2.params[15] * out2.params[5]) +  2  * out2.covar[4][5] / (out2.params[4] * out2.params[5]) +  2  * out2.covar[15][4] / (out2.params[15] * out2.params[4])
    #he1_5876_ew_err =  he1_5876_ew * np.sqrt((out2.perror[15]/out2.params[15])**2 + (out2.perror[5]/out2.params[5])**2 + (out2.perror[4]/out2.params[4])**2 +  covar_term)  


    ### oi BOTH lines  
    indices= [16, 41, 4, 6, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    oi_flux = 1.3 * measurements[0] 
    oi_err = 1.3 * measurements[1]

    #oi_flux = 1.3 * np.sqrt(2  * math.pi) *  out.params[16] * out.params[6] * out.params[4]    ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[16][6] / (out.params[16] * out.params[6]) + 2  * out.covar[4][6] / (out.params[4] * out.params[6]) + 2  * out.covar[16][4] / (out.params[16] * out.params[4])
    #oi_err =  oi_flux * np.sqrt((out.perror[16]/out.params[16])**2 + (out.perror[6]/out.params[6])**2 + (out.perror[4]/out.params[4])**2  +  covar_term )


    #oi_ew = 1.3 * np.sqrt(2  * math.pi) *  out2.params[16] * out2.params[6] * out2.params[4]    ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out2.covar[16][6] / (out2.params[16] * out2.params[6])  + 2  * out2.covar[16][4] / (out2.params[16] * out2.params[4])  + 2  * out2.covar[4][6] / (out2.params[4] * out2.params[6]) 
    #oi_ew_err =  oi_ew * np.sqrt((out2.perror[16]/out2.params[16])**2 + (out2.perror[6]/out2.params[6])**2 +  (out2.perror[4]/out2.params[4])**2 + covar_term)


    #### Ne III 3870 BLEND
    indices = [17, 42, 4, 1, 45] 
    measurements = error_2comp_gaussian(out.params, out.perror, out.covar, indices)
    neiii_flux = measurements[0] 
    neiii_err = measurements[1] 



    #neiii_flux = np.sqrt(2  * math.pi) *  out.params[1] * out.params[17] * out.params[4]    ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out.covar[1][17] / (out.params[1] * out.params[17]) + 2  * out.covar[4][17] / (out.params[4] * out.params[17]) + 2  * out.covar[1][4] / (out.params[1] * out.params[4])
    #neiii_err =  neiii_flux * np.sqrt( (out.perror[1]/out.params[1])**2 + (out.perror[17]/out.params[17])**2 +  (out.perror[4]/out.params[4])**2 + covar_term)


    #neiii_ew = np.sqrt(2  * math.pi) *  out2.params[1] * out2.params[17] * out2.params[4]     ## sqrt(2 * pi) * amplitude * sigma 
    #covar_term = 2  * out2.covar[1][17] / (out2.params[1] * out2.params[17]) +  2  * out2.covar[4][17] / (out2.params[4] * out2.params[17]) +   2  * out2.covar[1][4] / (out2.params[1] * out2.params[4])
    #neiii_ew_err =  neiii_ew * np.sqrt((out2.perror[1]/out2.params[1])**2 + (out2.perror[17]/out2.params[17])**2 + (out2.perror[4]/out2.params[4])**2  + covar_term)
    

    #### gather some metadata on the OIII lines 
    dat = asciitable.read(input_masterlist)  
    #dat2 = asciitable.read(input_masterlist + '.EWlist')


    z = dat['z']  
    foiii = dat['foiii']
    zstack = np.median(z)

    ew_oiii_rest = dat['EW_oiii_rest'] 
    flam_oiii = dat['flam_oiii'] 
    flam_ha = dat['flam_ha']
    flam_hb = dat['flam_hb']
    flam_hg = dat['flam_hg']
    flam_hd = dat['flam_hd'] 






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


    w=np.where(flam_oiii > 0) 
    #### since oiii ew is rest, these are all rest
    oiii_ew = np.median(ew_oiii_rest[w]) 
    hb_ew = oiii_ew * hb_flux/oiii_flux * np.median(flam_oiii[w]/flam_hb[w]) 
    ha_ew = oiii_ew * ha_flux / oiii_flux * np.median(flam_oiii[w]/flam_ha[w])  
    hg_ew = oiii_ew * hg_flux/ oiii_flux * np.median(flam_oiii[w] / flam_hg[w]) 
    hd_ew = oiii_ew * hd_flux /oiii_flux * np.median(flam_oiii[w]/flam_hd[w]) 

    oiii_ew_err = np.std(ew_oiii_rest[w]) / np.sqrt(np.size(ew_oiii_rest[w])) 

    #### ha_ew_err 
    ### each term is delta A/A, delta B/B, delta C/C. , where A = oiii_ew, B = ha_flux/oiii_flux, C = np.median(flam_oiii/flam_ha)  
    ### term2 would be (ha_flux/oiii_flux) * sqrt((ha_err/ha_flux)** + (oiii_err/oiii_flux)) / (ha_flux/oiii_flux)  
    term1 = oiii_ew_err / oiii_ew 
    term2 = np.sqrt( (ha_err/ha_flux)**2 + (oiii_err/oiii_flux)**2) 
    term3 = np.std(flam_oiii[w]/flam_ha[w]) /np.sqrt(np.size(ew_oiii_rest[w])) / np.median(flam_oiii[w]/flam_ha[w]) 
    ha_ew_err = ha_ew  * np.sqrt(term1**2 + term2**2 + term3**2) 

    term1 = oiii_ew_err / oiii_ew 
    term2 = np.sqrt( (hb_err/hb_flux)**2 + (oiii_err/oiii_flux)**2) 
    term3 = np.std(flam_oiii[w]/flam_hb[w]) /np.sqrt(np.size(ew_oiii_rest[w])) / np.median(flam_oiii[w]/flam_hb[w]) 
    hb_ew_err = hb_ew  * np.sqrt(term1**2 + term2**2 + term3**2)

    term1 = oiii_ew_err / oiii_ew 
    term2 = np.sqrt( (hg_err/hg_flux)**2 + (oiii_err/oiii_flux)**2) 
    term3 = np.std(flam_oiii[w]/flam_hg[w]) /np.sqrt(np.size(ew_oiii_rest[w])) / np.median(flam_oiii[w]/flam_hg[w]) 
    hg_ew_err = hg_ew  * np.sqrt(term1**2 + term2**2 + term3**2)   

    term1 = oiii_ew_err / oiii_ew 
    term2 = np.sqrt( (hd_err/hd_flux)**2 + (oiii_err/oiii_flux)**2) 
    term3 = np.std(flam_oiii[w]/flam_hd[w]) /np.sqrt(np.size(ew_oiii_rest[w])) / np.median(flam_oiii[w]/flam_hd[w]) 
    hd_ew_err = hd_ew  * np.sqrt(term1**2 + term2**2 + term3**2) 



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
    ew_stack = np.array([ha_ew, hb_ew, hg_ew, hd_ew, oiii_ew, -1, -1, -1, -1, -1]) 
    ew_stack_err = np.array([ha_ew_err,hb_ew_err, hg_ew_err, hd_ew_err, oiii_ew_err, -1, -1, -1, -1, -1]) 

    print out.params

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

   






        

    f, axarr = plt.subplots(1, 1,  figsize=(8, 10))

    axarr.plot(lam, flux, 'k-', ls='steps-mid') 
    axarr.plot(lam, model_fit, 'r--')
    #plt.axhline(y=0, color = 'b', ls = '--')
     

   # axarr[1].plot(lam, median_cont, 'k-', ls='steps-mid') 
   # axarr[1].plot(lam, model_fit2, 'r--') 

    plt.savefig(output_fig) 
    if showfig == True :
        plt.show()  






    
