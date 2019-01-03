from wisp_analysis import *
from mpfit import *

def stack_emline_model(pars, x):
    sigma_blue = pars[0]
    sigma_red = pars[1]
    sigma_green = pars[2]

    ha_amp = pars[3] 
    hb_amp = pars[4] 
    hg_amp = pars[5] 
    hd_amp = pars[6]
    oiii_5007_amp = pars[7] 
    oii_amp = pars[8] 
    sii_amp = pars[9] 
    he1_5876_amp = pars[10] 
    oi_amp = pars[11] 
    neiii_amp = pars[12]
    nii_6583_amp = pars[13] 
    he1_6678_amp = pars[14] 


    model = ha_amp * gaussian(x, 6564.6, sigma_red) + \
            hb_amp  * gaussian(x, 4862.7, sigma_green) + \
            hg_amp * gaussian(x, 4341.7, sigma_green) + \
            hd_amp * gaussian(x, 4102.9, sigma_blue) + \
            oiii_5007_amp * gaussian(x, 5008., sigma_green) + \
            oiii_5007_amp/3. * gaussian(x, 4960., sigma_green) +\
            oii_amp * gaussian(x, 3728., sigma_blue) +\
            sii_amp * gaussian(x, 6725., sigma_red) +\
            he1_5876_amp  * gaussian(x, 5877.2, sigma_red) +\
            oi_amp *  gaussian(x, 6302., sigma_red) + \
            oi_amp /3 * gaussian(x, 6365.5,sigma_red) + \
            neiii_amp * gaussian(x, 3870, sigma_blue) + \
            nii_6583_amp * gaussian(x, 6585.23, sigma_red) + \
            nii_6583_amp / 3 * gaussian(x, 6550., sigma_red)+\
            he1_6678_amp  * gaussian(x, 6680., sigma_red) 
    return model


    

def model_resid_stack(pars, fjac=None, lam = None, flux = None, err = None):
    model = stack_emline_model(pars, lam) 
    resid = (flux- model) / err 
   
    status = 0 
    return [status, resid]

    



def measure_stack(input_stack, output_meas, output_fig): 

    tab= asciitable.read(input_stack) 
    lam = tab['lam']
    flux = tab['flux_norm']
    err = tab['err']

    w=np.where(lam < 7000) 
    w=w[0] 
    lam = lam[w]
    flux = flux[w] 
    err = err[w]

    pguess= np.zeros(15) 
    pguess[0] = 50 
    pguess[1] = 50 
    pguess[2] = 50
    pguess[3] = 0.003
    pguess[4] = 0.001
    pguess[5] = 0.0008
    pguess[6] = 0.0003
    pguess[7] = 0.0055
    pguess[8] = 0.015
    pguess[9] = 0.0004
    pguess[10] = 0.0001
    pguess[11] = 0.0001
    pguess[12] = 0.0003 
    pguess[13] = 0.0000   #### set nii and he1 to zero, even though they are in the model and can be added. 
    pguess[14] = 0.0000

    npars = len(pguess) 
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} 
              for i in range(npars)]
    for i in range(npars): parinfo[i]['value'] = pguess[i]  
   
    parinfo[14]['fixed'] = 1
    parinfo[13]['fixed'] = 1

    fa = {'lam':lam, 'flux':flux, 'err':err} 

    out = mpfit(model_resid_stack, pguess, functkw=fa, parinfo = parinfo, quiet=True) 

    model_fit = stack_emline_model(out.params, lam)


    ### gather line fluxes 
    ### ha  
    ### note that nii  and he1 are set to zero.  This means that this is Ha + NII
    ### and I'm not sure if He 1 is included in the flux or resolved, or something in between 
    ha_flux = np.sqrt(2  * math.pi) *  out.params[3] * out.params[1]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[3][1] / (out.params[3] * out.params[1])
    ha_err =  ha_flux * np.sqrt( (out.perror[3]/out.params[3])**2 + (out.perror[1]/out.params[1])**2 + covar_term) 

    ##hb 
    hb_flux = np.sqrt(2  * math.pi) *  out.params[4] * out.params[2] ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[4][2] / (out.params[4] * out.params[2])
    hb_err =  hb_flux * np.sqrt( (out.perror[4]/out.params[4])**2 + (out.perror[2]/out.params[2])**2 + covar_term)  


    ##hg 
    hg_flux = np.sqrt(2  * math.pi) *  out.params[5] * out.params[2]     ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[5][2] / (out.params[5] * out.params[2])
    hg_err =  hg_flux * np.sqrt( (out.perror[5]/out.params[5])**2 + (out.perror[2]/out.params[2])**2 + covar_term)  

    ### hd 
    hd_flux = np.sqrt(2  * math.pi) *  out.params[6] * out.params[0]   ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[6][0] / (out.params[6] * out.params[0])
    hd_err =  hd_flux * np.sqrt( (out.perror[6]/out.params[6])**2 + (out.perror[0]/out.params[0])**2 + covar_term)

    ### oiii BOTH LINES
    oiii_flux = 1.3 * np.sqrt(2  * math.pi) *  out.params[7] * out.params[2] 
    covar_term = 2  * out.covar[7][2] / (out.params[7] * out.params[2]) 
    oiii_err =  oiii_flux * np.sqrt( (out.perror[7]/out.params[7])**2 + (out.perror[2]/out.params[2])**2 + covar_term) 

    ### oii
    oii_flux = np.sqrt(2  * math.pi) *  out.params[8] * out.params[0]   ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[8][0] / (out.params[8] * out.params[0])
    oii_err =  oii_flux * np.sqrt( (out.perror[8]/out.params[8])**2 + (out.perror[0]/out.params[0])**2 + covar_term)

    ### sii
    sii_flux = np.sqrt(2  * math.pi) *  out.params[9] * out.params[1]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[9][1] / (out.params[9] * out.params[1])
    sii_err =  sii_flux * np.sqrt((out.perror[9]/out.params[9])**2 + (out.perror[1]/out.params[1])**2 + covar_term) 


    #### he1 5876
    he1_5876_flux = np.sqrt(2  * math.pi) *  out.params[10] * out.params[1]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[10][1] / (out.params[10] * out.params[1])
    he1_5876_err =  he1_5876_flux * np.sqrt((out.perror[10]/out.params[10])**2 + (out.perror[1]/out.params[1])**2 + covar_term)  


    ### oi BOTH lines 
    oi_flux = 1.3 * np.sqrt(2  * math.pi) *  out.params[11] * out.params[1]    ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[11][1] / (out.params[11] * out.params[1])
    oi_err =  oi_flux * np.sqrt((out.perror[11]/out.params[11])**2 + (out.perror[1]/out.params[1])**2 + covar_term)


    #### Ne III 3870 BLEND
    neiii_flux = np.sqrt(2  * math.pi) *  out.params[12] * out.params[0]     ## sqrt(2 * pi) * amplitude * sigma 
    covar_term = 2  * out.covar[12][0] / (out.params[12] * out.params[0])
    neiii_err =  neiii_flux * np.sqrt((out.perror[12]/out.params[12])**2 + (out.perror[0]/out.params[0])**2 + covar_term)


    output = open(output_meas, 'w') 
    output.write('Ha+ NII            ' + str(ha_flux) + '  ' + str(ha_err) + '\n')   ### the ha_flux corresponds only to Ha + NII because the nii amp in the model is fixed at zero.
    output.write('Hb                 ' + str(hb_flux) + '  ' + str(hb_err) + '\n') 
    output.write('Hg + OIII          ' + str(hg_flux) + '  ' + str(hg_err) + '\n') 
    output.write('Hd                 ' + str(hd_flux) + '  ' + str(hd_err) + '\n') 
    output.write('OIII (both)        ' + str(oiii_flux) + '   ' +str(oiii_err) + '\n') 
    output.write('OII                ' + str(oii_flux) + '   ' +str(oii_err)  + '\n') 
    output.write('SII (both)         ' + str(sii_flux) + '   '  + str(sii_err) + '\n')  
    output.write('HeI 5876           ' + str(he1_5876_flux) + '   ' + str(he1_5876_err) +  '\n') 
    output.write('OI 6300+6363       ' + str(oi_flux) + '  ' + str(oi_err) + '\n') 
    output.write('Ne III 3869 blend  ' + str(neiii_flux) + '  ' + str(neiii_err) + '\n') 

    output.close() 


    #plt.axvline(x = 5007 - 100) 
    #plt.axvline(x = 5007 + 100)

    




    plt.plot(lam, flux, 'k-', ls='steps-mid') 
    plt.plot(lam, model_fit, 'r--')
    plt.savefig(output_fig) 

    plt.show() 






    
