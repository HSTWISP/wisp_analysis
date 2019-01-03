from wisp_analysis import * 
from scipy import interpolate
import time
import scipy.integrate as integrate 



def stack_spec(inlist, outstack, path_wisp = './', path_3dhst = './', bootstrap = None): 
    t0 = time.time()
    inlist_data = asciitable.read(inlist) 
    fieldname = inlist_data['Field']
    objid = inlist_data['ID'] 
    z_list = inlist_data['z'] 
    foiii_list = inlist_data['foiii'] 
   


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
            if os.path.exists(path_3dhst+  '/' + fieldname[i] + '_output_alaina-mzr/'): 
                specfile = path_3dhst + '/' + fieldname[i] + '_output_alaina-mzr/fitdata/' + fieldname[i]+  '_' +  '{:05d}'.format(objid[i]) + '_fitspec.dat' 
                #catalog  = path_3dhst + '/' + fieldname[i] + '_output_alaina-mzr/' + fieldname[i] + '_catalog_alaina-mzr.dat' 
            elif os.path.exists( path_3dhst + '/' + fieldname[i] + '_output_marc-mzr/'): 
                specfile = path_3dhst + '/' + fieldname[i] + '_output_marc-mzr/fitdata/' + fieldname[i] + '_' + '{:05d}'.format(objid[i]) + '_fitspec.dat' 
                #catalog  = path_3dhst + '/' + fieldname[i] + '_output_marc-mzr/' + fieldname[i] + '_catalog_marc-mzr.dat' 
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


                ##### locate redshift and oiii flux 
                #cat_objid = catdata['ObjID'] 
                #cat_z = catdata['redshift'] 
                #cat_foiii = catdata['oiii_flux'] 

                #w=np.where(cat_objid == objid[i]) 
                f_oiii = foiii_list[i]
                z = z_list[i]
                
                ### why the hell are some objects cataloged more than once? 
                ### I can't remember, but use the last one.
                ###  this should be gone. 
                #if np.size(f_oiii) > 1: 
                #    f_oiii = f_oiii[-1] 
                #if np.size(z) > 1 : 
                #    z = z[-1] 

                #### de-redshift the spectrum and insert into 2d array 
                #print np.size(flam_spec), np.size(cont_spec), objid[i], f_oiii, z
                flam_norm = (flam_spec - cont_spec) / f_oiii  ### need to think about whether I need a 1+z here.
                lam_spec_rest = lam_spec / (1+z) 
                f = interpolate.interp1d(lam_spec_rest, flam_norm, fill_value=-1, bounds_error=False)
                f2 = interpolate.interp1d(lam_spec_rest, mask_spec, fill_value = -1, bounds_error = False, kind = 'nearest') 
                flam_interp = f(lam_stack)
                flam_interp = flam_interp * (1+z) ### this is required to de-redshift and maintain normalization. 
                
                #w=np.where( (lam_stack > 4907)  & (lam_stack < 5107)) 
                #print integrate.trapz(flam_interp[w], lam_stack[w]) 


                mask_interp = f2(lam_stack) 
                w=np.where(mask_interp > 0) 
                flam_interp[w] = -1 
                stack_frame[i, :] = flam_interp
                #plt.plot(lam_stack, flam_interp, ls = 'steps-mid') 
                #plt.show()
                #plt.clf()
            else:  
                if os.path.exists(specfile) == False :
                    print 'Could not find : ' + specfile 
               # if os.path.exists(catalog) == False :
                  #print 'Could not find :' + catalog 

   # hdu = fits.PrimaryHDU(stack_frame)
   # hdu1 = fits.HDUList([hdu]) 
   # hdu1.writeto('test.fits')

    #plt.imshow(stack_frame) 
    #plt.show() 

    t1 = time.time() 

    print str(t1 - t0) + ' seconds to read and de-redshift spectra'     
    for i in np.arange(nlam): 
        w=np.where(stack_frame[:, i] > -1) 
        w=w[0]
        stack_spec[i] = np.median(stack_frame[w, i])   
        stack_err[i] = np.std(stack_frame[w, i]) / np.sqrt(np.size(w))
    
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
    outfile.write('lam       flux_norm     err  \n') 
    for a, b, c in zip(lam_stack, stack_spec, stack_err): 
        outfile.write(str(a) + '  ' + str(b) + '  '  +  str(c) + '\n') 
    outfile.close() 
         

    #plt.plot(lam_stack, stack_spec, ls='steps-mid')  




    plt.show()








        







