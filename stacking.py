from wisp_analysis import * 
from scipy import interpolate
import time
import scipy.integrate as integrate 

import astropy.units as u
from astropy.cosmology import Planck13 as cosmo


def robust_sigma(y, zero = False):
    ### checked against IDL for one case and it was numerically identical 
    #### direct transltion from IDL astronomy users library. 
    eps = 1.0e-20 
    
    if zero == True :
        y0 = 0 
    else : 
        y0 = np.median(y) 


    mad = np.median(np.abs(y - y0))/0.6745 

    if mad < eps: 
        mad = np.mean(abs(y - y0))/0.80 
    
    if mad < eps : 
        out = 0. 

    if mad >= eps :
        u= (y  - y0) / (6. * mad) 
        uu = u * u 
        q = np.where(uu < 1.0) 
        count = np.size(q) 
        if count < 3: 
            out = -1 

        N = np.sum(np.isfinite(y)) 
        numerator = np.sum( (y[q] - y0)**2 * (1 - uu[q])**4 ) 
        den1= np.sum( (1 - uu[q]) * (1 - 5. * uu[q])) 
        siggma = N * numerator  / (den1 * (den1 - 1)) 
        out = [np.sqrt(siggma), q] 

    return out 



def resistant_mean(y, Cut) : 
    
    MADscale = 0.6745 
    MADscale2 = 0.80 
    MADlim = 1.0e-24 

    sigcoeff = [ -0.15405, 0.90723, -0.23584, 0.020142 ]
  

    Npts = np.size(y) 
    ymed = np.median(y) 
    absdev = np.abs(y - ymed) 
    medabsdev = np.median(absdev) / MADscale 
    if medabsdev < MADlim : 
        finite_check = np.isfinite(absdev)
        medabsdev = np.mean(absdev[finite_check])/MADscale2 

    Cutoff = Cut * medabsdev 

    goodvec = np.where(absdev < Cutoff) 
    Num_Good = np.size(goodvec) 
    if Num_Good <= 0 : 
        mean  = [float('NaN')]
        return mean 


    
    goodpts = y[goodvec] 
    mean = np.mean(goodpts) 
    sigma = np.sqrt( np.sum( (goodpts - mean)**2 ) / Num_Good ) 
  
    #### apparently this is how we compensate for sigma truncation, in some algornthm 
    if Cut > 1.0: 
        SC = Cut 
    else : 
        SC = 1.0 

    if SC < 4.50 : 
        sigma = sigma/(sigcoeff[0]  + sigcoeff[1] * SC + sigcoeff[2] * SC**2 + sigcoeff[3]* SC**3) 

    Cutoff = Cut  * sigma 

    goodvec = np.where(absdev < Cutoff)
    Num_Good = np.size(goodvec) 
    goodpts = y[goodvec] 
    mean = np.mean(goodpts) 
    sigma = np.sqrt( np.sum( (goodpts - mean)**2 ) / Num_Good )

    if Cut > 1.0: 
        SC = Cut 
    else : 
        SC = 1.0 

    if SC < 4.50 : 
        sigma = sigma/(sigcoeff[0]  + sigcoeff[1] * SC + sigcoeff[2] * SC**2 + sigcoeff[3] *SC**3) 


    sigma = sigma/np.sqrt(Num_Good-1)    

    return [mean, sigma, goodvec] 
    



def stack_spec(inlist, outstack, path_wisp = './', path_3dhst = './', bootstrap = None, dump_EWs = False): 
    t0 = time.time()
    inlist_data = asciitable.read(inlist) 
    fieldname = inlist_data['Field']
    objid = inlist_data['ID'] 
    z_list = inlist_data['z'] 
    foiii_list = inlist_data['foiii'] 

    z_median = np.median(z_list) 
   


    ### find these in the data 
    #f_oiii = inlist_data['f_oiii'] 
    #z = inlist_data['z'] 


    ngals = np.size(objid) 

    ### define arrays to hold output stack and 2d pre-stack "image"

    dlam = 7 ## 7 A bins in the rest frame are reasonable. 
    lam_blue = 3250. 
    lam_red = 7200. 
    nlam = int((lam_red  -lam_blue) / dlam) 
    lam_stack = np.arange(nlam) * dlam + lam_blue

    stack_frame = np.zeros( (ngals, nlam))  - 1   ## fill these placeholders with -1's to mark empty data 
    stack_spec = np.zeros(nlam) - 1 
    stack_err = np.zeros(nlam) - 1 
    #cont_frame = np.zeros( (ngals, nlam)) - 1 
    #econt_frame = np.zeros( (ngals, nlam)) - 1
    #cont_stack_med = np.zeros(nlam) -1 
    #cont_stack_err = np.zeros(nlam) -1 
    #cont_stack_mean = np.zeros(nlam)-1 





    oiii_ew_obs = []
    oiii_ew_lim = [] 
    ### fill the stack frame with one row per spectrum. 
    for i in  np.arange(ngals): 
        
        #### first locate all data
       
        #### 1 deterimine if the source is in WISP or not. 
        if fieldname[i][0:3] == 'Par' : 
            if os.path.exists(path_wisp +  '/' + fieldname[i] + '_output_alaina-mzr/') : 
                specfile = path_wisp + '/' + fieldname[i] + '_output_alaina-mzr/fitdata/' + fieldname[i] + '_BEAM_' + str(objid[i]) + '_fitspec.dat' 
               # catalog  = path_wisp  + '/' + fieldname[i] + '_output_alaina-mzr/' + fieldname[i] + 'list_catalog_alaina-mzr.dat' 

            elif os.path.exists(path_wisp +'/' + fieldname[i] + '_output_marc-mzr/') : 
                specfile =  path_wisp + '/' + fieldname[i] + '_output_marc-mzr/fitdata/' + fieldname[i] + '_BEAM_' + str(objid[i]) + '_fitspec.dat' 
               # catalog  =  path_wisp + '/' + fieldname[i] + '_output_marc-mzr/' + fieldname[i] + 'list_catalog_marc-mzr.dat'  

            else : 
                specfile = None 
                #catalog = None 
                print 'Could not find fit data directory for ' + fieldname[i] 

        #### if not WISP, look for the 3D HST data 
        else :  
            if os.path.exists(path_3dhst+  '/' + fieldname[i] + '_output_final/'): 
                specfile = path_3dhst + '/' + fieldname[i] + '_output_final/fitdata/' + fieldname[i]+  '_' +  '{:05d}'.format(objid[i]) + '_fitspec.dat' 
                #catalog  = path_3dhst + '/' + fieldname[i] + '_output_alaina-mzr/' + fieldname[i] + '_catalog_alaina-mzr.dat' 
            #elif os.path.exists( path_3dhst + '/' + fieldname[i] + '_output_marc-mzr/'): 
            #    specfile = path_3dhst + '/' + fieldname[i] + '_output_marc-mzr/fitdata/' + fieldname[i] + '_' + '{:05d}'.format(objid[i]) + '_fitspec.dat' 
            #    #catalog  = path_3dhst + '/' + fieldname[i] + '_output_marc-mzr/' + fieldname[i] + '_catalog_marc-mzr.dat' 
            else :
                specfile = None 
                #catalog = None 
                print 'Could not find fit data directory for ' + fieldname[i]   

        if specfile is not None: 
            if os.path.exists(specfile) : 
                specdata = asciitable.read(specfile, fill_values=[('--', '-99')]) 
                #catdata = asciitable.read(catalog, fill_values=[('--', '-99')])  

                ##### get spectrum 
                lam_spec = specdata['Lam'] 
                flam_spec = specdata['Flam']
                cont_spec = specdata['Contmodel'] 
                mask_spec = specdata['Masked'] 
                contam_spec = specdata['Contam'] 
                flam_err_spec = specdata['Flam_err']
                


                if dump_EWs == True: 
                    #### used for estimating EWs. 
                    cont_no_contam = cont_spec - contam_spec 
                    snr_spec = cont_no_contam / flam_err_spec  
                    



                ##### locate redshift and oiii flux 
                #cat_objid = catdata['ObjID'] 
                #cat_z = catdata['redshift'] 
                #cat_foiii = catdata['oiii_flux'] 

                #w=np.where(cat_objid == objid[i]) 
                f_oiii = foiii_list[i]
                z = z_list[i]
                lam_spec_rest = lam_spec / (1+z) 

                if dump_EWs == True:
                    w=np.where( np.abs(lam_spec_rest - 5007) == np.min(np.abs(lam_spec_rest - 5007)))
                    w=w[0][0]

                    cont_oiii = cont_no_contam[w] 
                    contsnr_oiii = snr_spec[w] 
                    flam_err_oiii = flam_err_spec[w]
                    if contsnr_oiii > 0.5 :
                        ew_i = f_oiii / cont_oiii
                        oiii_ew_obs.append(f_oiii/cont_oiii) 
                        oiii_ew_lim.append(contsnr_oiii) 
                    else: 
                        oiii_ew_obs.append(f_oiii/(2 * flam_err_oiii))
                        oiii_ew_lim.append(-1)    
                


                #### de-redshift the spectrum and insert into 2d array 
                #print np.size(flam_spec), np.size(cont_spec), objid[i], f_oiii, z
                flam_norm = (flam_spec - cont_spec) / f_oiii 
                #flam_norm2 = (flam_spec - contam_spec) / (cont_spec - contam_spec) ### clean flam and cont-model from contam and normalize
                #err_norm2 = flam_err_spec/cont_spec 
                f = interpolate.interp1d(lam_spec_rest, flam_norm, fill_value=-1, bounds_error=False)
                f2 = interpolate.interp1d(lam_spec_rest, mask_spec, fill_value = -1, bounds_error = False, kind = 'nearest')
                #f3 = interpolate.interp1d(lam_spec_rest, flam_norm2, fill_value = -1, bounds_error = False) 
                #f4 = interpolate.interp1d(lam_spec_rest, err_norm2, fill_value = -1, bounds_error = False) 


                flam_interp = f(lam_stack)
                flam_interp = flam_interp * (1+z) ### this is required to de-redshift and maintain normalization. 

                ### because the spectrum has been continuum normalized, a 1+z factor is not necessary.  
               # cont_interp  = f3(lam_stack)     #  * (1 + z) * 4 * 3.14159 * cosmo.luminosity_distance(z).cgs.value**2
               # econt_interp = f4(lam_stack)     #  * (1 + z) * 4 * 3.14159 * cosmo.luminosity_distance(z).cgs.value**2

                ### integrating these spectra cont. normalized against the observed wavelegths will give observed EWs.  
                ### integrating the cont normalized spectra against rest wavelengths will give rest EWs. 




                #w=np.where( (lam_stack > 4907)  & (lam_stack < 5107)) 
                #print integrate.trapz(flam_interp[w], lam_stack[w]) 


                mask_interp = f2(lam_stack) 
                w=np.where(mask_interp > 0) 
                flam_interp[w] = -1 
                #cont_interp[w] = -1 
                stack_frame[i, :] = flam_interp
                #cont_frame[i, :] = cont_interp 
                #econt_frame[i, :] = econt_interp

                #plt.plot(lam_stack, flam_interp, ls = 'steps-mid') 
                #plt.show()
                #plt.clf()
            else:  
                if os.path.exists(specfile) == False :
                    print 'Could not find : ' + specfile 
               # if os.path.exists(catalog) == False :
                  #print 'Could not find :' + catalog 

    #hdu = fits.PrimaryHDU(cont_frame)
    #hdu1 = fits.HDUList([hdu]) 
    #hdu1.writeto('test.fits')


    



    #plt.imshow(stack_frame) 
    #plt.show() 

    t1 = time.time() 

    print str(t1 - t0) + ' seconds to read and de-redshift spectra'     
    
    for i in np.arange(nlam): 
        w=np.where(stack_frame[:, i] > -1) 
        w=w[0]
        stack_spec[i] = np.median(stack_frame[w, i])   
        stack_err[i] = np.std(stack_frame[w, i]) / np.sqrt(np.size(w))
   
        #w=np.where(cont_frame[:, i] > -1)   
        ### convert stacked spectra back into the flux units for the median redshift. 
        #cont_slice = cont_frame[w, i] 
        #econt_slice = econt_frame[w, i] 

        #cont_stack_med[i] = np.median(cont_slice)  # / (1+ z_median) / (4 * 3.14159 * cosmo.luminosity_distance(z_median).cgs.value**2)
        
        #out = resistant_mean(cont_slice, 2)
        #if len(out) == 3 : 
        #    cont_stack_mean[i] = out[0]   #/ (1+ z_median) / (4 * 3.14159 * cosmo.luminosity_distance(z_median).cgs.value**2)
        #    cont_stack_err[i]   = out[1]  #/ (1+ z_median) / (4 * 3.14159 * cosmo.luminosity_distance(z_median).cgs.value**2)
        #else : 
        #    cont_stack_mean[i] = 0.
        #    cont_stack_err[i]  = 0.


        #cont_stack_mean[i] = np.sum(cont_slice[goodvec] / econt_slice[goodvec]**2) / np.sum(1/econt_slice[goodvec]**2) / (1+ z_median) / (4 * 3.14159 * cosmo.luminosity_distance(z_median).cgs.value**2)
        #cont_stack_err[i]  = cont_sig / np.sqrt(np.size(w)) / (1+ z_median) / (4 * 3.14159 * cosmo.luminosity_distance(z_median).cgs.value**2)
        





    
    t2 = time.time() 
    print str(t2 - t1)  + ' seconds  to take the median' 

    if bootstrap == True : 
        nstraps = 50  ### later make nstraps = ngals. 
        bootstrap_array_2d = np.zeros( (nstraps, nlam))  - 1
        for i in np.arange(nstraps) : 
            stack_frame_bootstrap = np.zeros( (ngals, nlam))  - 1 
            ### generage a set of indices for the individual spectra that we will draw to sample
            sample_indicies = np.round(np.random.uniform(0, ngals-1, ngals)).astype(int) 
            #### insert these individual into the bootstrap stack 
            for j in np.arange(ngals):
                stack_frame_bootstrap[j, :] = stack_frame[sample_indicies[j], :] 

            ### re-write the bootstrapped spectrum into a 2d array, one spectrum per row
            for j in np.arange(nlam): 
                w=np.where(stack_frame_bootstrap[:, j] > -1) 
                w=w[0]
                bootstrap_array_2d[i, j] = np.median(stack_frame_bootstrap[w, j])   
            
        for i in np.arange(nlam) : 
            stack_err[i] = np.std(bootstrap_array_2d[:, i]) 





    print 'writing output file: ' +  outstack 
    outfile = open(outstack, 'w') 
    outfile.write('lam       flux_norm     err\n') 
    for a, b, c in zip(lam_stack, stack_spec, stack_err):  
        outfile.write(str(a) + '  ' + str(b) + '  '  +  str(c) +  '\n') 
    outfile.close() 


    if dump_EWs == True: 
        oiii_ew_obs = np.array(oiii_ew_obs) 
        oiii_ew_lim = np.array(oiii_ew_lim)
        fieldname = np.array(fieldname) 
        objid = np.array(objid) 

        outdata = [fieldname, objid, oiii_ew_obs, oiii_ew_lim] 
        names = ['Field', 'ID', 'OIII_EW_OBS', 'OIII_CONT_SNR_OR_FLAG'] 
        asciitable.write(outdata, inlist + '.EWlist', names  = names,  overwrite=True) 



         

    #plt.plot(lam_stack, stack_spec, ls='steps-mid')  




    #plt.show()








        







