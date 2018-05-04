from wisp_analysis import * 
from scipy import interpolate

def stack_spec(inlist, outstack, norm_mode = None, stack_mode = None): 

    inlist_data = asciitable.read(inlist) 
    fieldname = inlist_data['col1']
    objid = inlist_data['col2']  
    
    ### find these in the data 
    #f_oiii = inlist_data['f_oiii'] 
    #z = inlist_data['z'] 




    ngals = np.size(objid) 

    ### define arrays to hold output stack and 2d pre-stack "image"

    dlam = 7 ## 7 A bins in the rest frame are reasonable. 
    lam_blue = 3000. 
    lam_red = 7000. 
    nlam = int((lam_red  -lam_blue) / dlam) 
    lam_stack = np.arange(nlam) * dlam + lam_blue

    stack_frame = np.zeros( (ngals, nlam))  - 1   ## fill these placeholders with -1's to mark empty data 
    stack_spec = np.zeros(nlam) - 1 

    ### fill the stack frame with one row per spectrum. 
    for i in  np.arange(ngals): 
        
        #### first locate all data
       
        #### 1 deterimine if the source is in WISP or not. 
        if fieldname[i][0:3] == 'Par' : 
            if os.path.exists('/Volumes/Thunderbay/wisps/mzr_refit/' + fieldname[i] + '_output_alaina-mzr/') : 
                specfile = '/Volumes/Thunderbay/wisps/mzr_refit/' + fieldname[i] + '_output_alaina-mzr/fitdata/' + fieldname[i] + '_BEAM_' + str(objid[i]) + '_fitspec.dat' 
                catalog  = '/Volumes/Thunderbay/wisps/mzr_refit/' + fieldname[i] + '_output_alaina-mzr/' + fieldname[i] + 'list_catalog_alaina-mzr.dat' 

            elif os.path.exists('/Volumes/Thunderbay/wisps/mzr_refit/' + fieldname[i] + '_output_marc-mzr/') : 
                specfile = '/Volumes/Thunderbay/wisps/mzr_refit/' + fieldname[i] + '_output_marc-mzr/fitdata/' + fieldname[i] + '_BEAM_' + str(objid[i]) + '_fitspec.dat' 
                catalog  = '/Volumes/Thunderbay/wisps/mzr_refit/' + fieldname[i] + '_output_marc-mzr/' + fieldname[i] + 'list_catalog_marc-mzr.dat'  

            else : 
                specfile = None 
                catalog = None 
                print 'Could not find fit data directory for ' + fieldname[i] 

        #### if not WISP, look for the 3D HST data 
        else :  
            if os.path.exists('/Volumes/Thunderbay/3DHST/mzr_refit/' + fieldname[i] + '_output_alaina-mzr/'): 
                specfile = '/Volumes/Thunderbay/3DHST/mzr_refit/' + fieldname[i] + '_output_alaina-mzr/fitdata/' + fieldname[i]+  '_' +  '{:05d}'.format(objid[i]) + '_fitspec.dat' 
                catalog  = '/Volumes/Thunderbay/3DHST/mzr_refit/' + fieldname[i] + '_output_alaina-mzr/' + fieldname[i] + '_catalog_alaina-mzr.dat' 
            elif os.path.exists('/Volumes/Thunderbay/3DHST/mzr_refit/' + fieldname[i] + '_output_marc-mzr/'): 
                specfile = '/Volumes/Thunderbay/3DHST/mzr_refit/' + fieldname[i] + '_output_marc-mzr/fitdata/' + fieldname[i] + '_' + '{:05d}'.format(objid[i]) + '_fitspec.dat' 
                catalog  = '/Volumes/Thunderbay/3DHST/mzr_refit/' + fieldname[i] + '_output_marc-mzr/' + fieldname[i] + '_catalog_marc-mzr.dat' 
            else :
                specfile = None 
                catalog = None 
                print 'Could not find fit data directory for ' + fieldname[i]   

        if specfile is not None: 
            if os.path.exists(specfile)  & os.path.exists(catalog) : 
                specdata = asciitable.read(specfile, fill_values=[('--', '-99')]) 
                catdata = asciitable.read(catalog, fill_values=[('--', '-99')])  

                ##### get spectrum 
                lam_spec = specdata['Lam'] 
                flam_spec = specdata['Flam']
                cont_spec = specdata['Contmodel'] 
                mask_spec = specdata['Masked'] 


                #### locate redshift and oiii flux 
                cat_objid = catdata['ObjID'] 
                cat_z = catdata['redshift'] 
                cat_foiii = catdata['oiii_flux'] 

                w=np.where(cat_objid == objid[i]) 
                f_oiii = cat_foiii[w[0]]
                z = cat_z[w[0]]
                
                ### why the hell are some objects cataloged more than once? 
                if np.size(f_oiii) > 1: 
                    f_oiii = f_oiii[0] 
                if np.size(z) > 1 : 
                    z = z[0] 

                #### de-redshift the spectrum and insert into 2d array 
                #print np.size(flam_spec), np.size(cont_spec), objid[i], f_oiii, z
                flam_norm = (flam_spec - cont_spec) / f_oiii  ### need to think about whether I need a 1+z here.
                lam_spec_rest = lam_spec / (1+z) 
                f = interpolate.interp1d(lam_spec_rest, flam_norm, fill_value=-1, bounds_error=False)
                f2 = interpolate.interp1d(lam_spec_rest, mask_spec, fill_value = -1, bounds_error = False, kind = 'nearest') 
                flam_interp = f(lam_stack) 
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
                if os.path.exists(catalog) == False :
                    print 'Could not find :' + catalog 

    #hdu = fits.PrimaryHDU(stack_frame)
    #hdu1 = fits.HDUList([hdu]) 
    #hdu1.writeto('test.fits')

    #plt.imshow(stack_frame) 
    #plt.show() 
    print 'am I doing something new?' 
    for i in np.arange(nlam): 
        w=np.where(stack_frame[:, i] > -1) 
        stack_spec[i] = np.median(stack_frame[w, i])  

    plt.plot(lam_stack, stack_spec, ls='steps-mid')   
    plt.show()








        







