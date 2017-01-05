from wisp import *
def trim_spec(tbdata_blue, tbdata_red,config_pars):
    #ext = hdu.index_of(beam_name) 
    #tbdata = hdu[ext].data  

    if tbdata_blue != None: 
        lam_spec_blue = tbdata_blue['lambda'] 
        flux_spec_blue = tbdata_blue['flux'] 
        error_spec_blue = tbdata_blue['ferror']  
        contam_spec_blue = tbdata_blue['contam']
        zero_spec_blue = tbdata_blue['zero']

        ##### trim the spectrum  ##### 
        ### step 0:  remove first and last 3 pixels  
        ### this is usually redundant, but when spectra fall off edge of 
        ### detector, the latter steps may not catch the messes. 
        lam_spec_blue = lam_spec_blue[3:-3]
        flux_spec_blue = flux_spec_blue[3:-3] 
        error_spec_blue = error_spec_blue[3:-3] 
        contam_spec_blue = contam_spec_blue[3:-3] 
        zero_spec_blue = zero_spec_blue[3:-3] 

        ### only fit fintite data
        w = np.isfinite(flux_spec_blue)  
        lam_spec_blue = lam_spec_blue[w] 
        flux_spec_blue = flux_spec_blue[w]  
        error_spec_blue = error_spec_blue[w]
        contam_spec_blue = contam_spec_blue[w] 
        zero_spec_blue = zero_spec_blue[w]

        w = np.isfinite(error_spec_blue)  
        lam_spec_blue = lam_spec_blue[w] 
        flux_spec_blue = flux_spec_blue[w]  
        error_spec_blue = error_spec_blue[w]
        contam_spec_blue = contam_spec_blue[w]
        zero_spec_blue = zero_spec_blue[w]

        #### clip the edges in wavelength where stuff gets crazy (i.e. low throughput)  
        w=np.where( (lam_spec_blue > config_pars['lambda_min']) & (lam_spec_blue < config_pars['transition_wave'])) 
        lam_spec_blue = lam_spec_blue[w] 
        flux_spec_blue = flux_spec_blue[w] 
        error_spec_blue =error_spec_blue[w]
        contam_spec_blue = contam_spec_blue[w] 
        zero_spec_blue = zero_spec_blue[w] 

    if tbdata_red != None: 
         ###  repeat for red
        lam_spec_red = tbdata_red['lambda'] 
        flux_spec_red = tbdata_red['flux'] 
        error_spec_red = tbdata_red['ferror']  
        contam_spec_red = tbdata_red['contam']
        zero_spec_red = tbdata_red['zero']
        #del hdu[ext].data 

        ##### trim the spectrum  ##### 
        ### step 0:  remove first and last 3 pixels  
        ### this is usually redundant, but when spectra fall off edge of 
       ### detector, the latter steps may not catch the messes. 
        lam_spec_red = lam_spec_red[3:-3]
        flux_spec_red= flux_spec_red[3:-3] 
        error_spec_red = error_spec_red[3:-3] 
        contam_spec_red = contam_spec_red[3:-3] 
        zero_spec_red = zero_spec_red[3:-3]

       ### only fit fintite data
        w = np.isfinite(flux_spec_red)  
        lam_spec_red = lam_spec_red[w] 
        flux_spec_red = flux_spec_red[w]  
        error_spec_red = error_spec_red[w]
        contam_spec_red= contam_spec_red[w] 
        zero_spec_red = zero_spec_red[w] 

        w = np.isfinite(error_spec_red)  
        lam_spec_red = lam_spec_red[w] 
        flux_spec_red = flux_spec_red[w]  
        error_spec_red = error_spec_red[w]
        contam_spec_red = contam_spec_red[w]
        zero_spec_red = zero_spec_red[w]

        #### clip the edges in wavelength where stuff gets crazy (i.e. low throughput)  
        w=np.where( (lam_spec_red > config_pars['transition_wave']) & (lam_spec_red < config_pars['lambda_max'])) 
        lam_spec_red = lam_spec_red[w] 
        flux_spec_red = flux_spec_red[w] 
        error_spec_red =error_spec_red[w]
        contam_spec_red = contam_spec_red[w] 
        zero_spec_red = zero_spec_red[w]

    ##### concatenate.
    if tbdata_blue == None: 
        lam_spec = lam_spec_red 
        flux_spec = flux_spec_red
        error_spec = error_spec_red
        contam_spec = contam_spec_red 
        zero_spec = zero_spec_red
    if tbdata_red == None: 
        lam_spec = lam_spec_blue
        flux_spec = flux_spec_blue 
        error_spec = error_spec_blue 
        contam_spec = contam_spec_blue 
        zero_spec = zero_spec_blue
    
    if (tbdata_red != None) & (tbdata_blue != None): 
        lam_spec = np.append(lam_spec_blue, lam_spec_red) 
        flux_spec = np.append(flux_spec_blue, flux_spec_red) 
        error_spec = np.append(error_spec_blue, error_spec_red) 
        contam_spec = np.append(contam_spec_blue, contam_spec_red) 
        zero_spec = np.append(zero_spec_blue, zero_spec_red) 


    ##### removed masked regions
    w=np.where( (lam_spec < config_pars['mask_region1'][0]) |  (lam_spec > config_pars['mask_region1'][1]))
    lam_spec = lam_spec[w]
    flux_spec  = flux_spec[w] 
    error_spec = error_spec[w]
    contam_spec = contam_spec[w] 
    zero_spec = zero_spec[w]


    w=np.where( (lam_spec < config_pars['mask_region2'][0]) |  (lam_spec > config_pars['mask_region2'][1]))
    lam_spec = lam_spec[w]
    flux_spec  = flux_spec[w] 
    error_spec = error_spec[w]
    contam_spec = contam_spec[w] 
    zero_spec = zero_spec[w]
    
    w=np.where( (lam_spec < config_pars['mask_region3'][0]) |  (lam_spec > config_pars['mask_region3'][1]))
    lam_spec = lam_spec[w]
    flux_spec  = flux_spec[w] 
    error_spec = error_spec[w]
    contam_spec = contam_spec[w]
    zero_spec = zero_spec[w] 
    
    

    #### this should be unncessary.... 
    #w=np.where(contam_spec < 0) 
    #contam_spec[w] = 0
    

    return [lam_spec, flux_spec, error_spec, contam_spec, zero_spec] 


