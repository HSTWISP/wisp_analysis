from wisp_analysis import *


def initialize_arrays(data, bluecut, redcut):
    """Define and trim initial arrays from a table of input data.

    Define the wavelength, flux, flux error, contamination and zeroth order 
    flags for an object's spectrum in one of the grisms. The first and 
    last 3 pixels of the spectrum are removed to avoid the mess caused when 
    spectra fall off the edge of the detector. This step is usually 
    redundant, but the trimming and masking steps may not catch these 
    problems at the edges.

        ### mask edges in wavelength where stuff gets crazy (i.e. low throughput)
    The arrays are then trimmed in wavelength to remove the areas where 
    stuff gets crazy (i.e. low throughput). The elements corresponding to
    wavelengths outside of range defined by (bluecut,redcut) are removed.    

    Args:
        data (astropy.io.ascii table): table of data from the 1-D dat file
        bluecut (float): the minimum wavelength 
        redcut (float): the maximum wavelength

    Returns:
        (tuple): tuple containing:
            lam (float): wavelength array
            flux (float): flux array
            error (float): flux error array
            contam (float): contamination array
            zeros (int): zeroth order flags
    """
    lam = data['lambda'][3:-3]
    flux = data['flux'][3:-3]
    error = data['ferror'][3:-3]
    contam = data['contam'][3:-3]
    zeros = data['zero'][3:-3]

    cut = (lam > bluecut) & (lam < redcut)
    lam = lam[cut]
    flux = flux[cut]
    error = error[cut]
    contam = contam[cut]
    zeros = zeros[cut]

    return lam,flux,error,contam,zeros


def newmask(array, maskedarray):
    newarray = np.ma.masked_where(np.ma.getmask(maskedarray), array)
    return newarray


def trim_spec(tbdata_blue, tbdata_red, config_pars, mask_zeros=False, return_masks=False):
    """Create one array of spectra, etc. from multiple grism 1d files
    
    If masking: 
        don't need to mask wavelenght or zeroth orders, keep these for 
        writing out file

    create an array keeping track of the masked regions

    Args:
        tbdata_blue ():
        tbdata_red ():
        config_pars ():
        return_masks (Optional[bool]):

    Returns:
    """
    ### check each grism separately
    if tbdata_blue is not None:
        bluecut = config_pars['lambda_min']
        redcut = config_pars['transition_wave']
        specdata = initialize_arrays(tbdata_blue, bluecut, redcut)
        lam_spec_blue = specdata[0] 
        flux_spec_blue = specdata[1]
        error_spec_blue = specdata[2]
        contam_spec_blue = specdata[3]
        zero_spec_blue = specdata[4]

        ### only fit finite data
        flux_spec_blue = np.ma.masked_where(np.logical_or(~(np.isfinite(flux_spec_blue)),~(np.isfinite(error_spec_blue))), flux_spec_blue)
        error_spec_blue = newmask(error_spec_blue, flux_spec_blue)
        contam_spec_blue = newmask(contam_spec_blue, flux_spec_blue)

    if tbdata_red is not None: 
        bluecut = config_pars['transition_wave']
        redcut = config_pars['lambda_max']
        specdata = initialize_arrays(tbdata_red, bluecut, redcut)
        lam_spec_red = specdata[0]
        flux_spec_red = specdata[1]
        error_spec_red = specdata[2]
        contam_spec_red = specdata[3]
        zero_spec_red = specdata[4]

        ### only fit finite data
        flux_spec_red = np.ma.masked_where(np.logical_or(~(np.isfinite(flux_spec_red)),~(np.isfinite(error_spec_red))), flux_spec_red)
        error_spec_red = newmask(error_spec_red, flux_spec_red)
        contam_spec_red = newmask(contam_spec_red, flux_spec_red)

        ### mask edges in wavelength where stuff gets crazy (i.e. low throughput)  
     
    ### concatenate.
    if tbdata_blue is None: 
        lam_spec = lam_spec_red 
        flux_spec = flux_spec_red
        error_spec = error_spec_red
        contam_spec = contam_spec_red 
        zero_spec = zero_spec_red
    if tbdata_red is None: 
        lam_spec = lam_spec_blue
        flux_spec = flux_spec_blue 
        error_spec = error_spec_blue 
        contam_spec = contam_spec_blue 
        zero_spec = zero_spec_blue
    
    if (tbdata_red is not None) & (tbdata_blue is not None): 
        lam_spec = np.append(lam_spec_blue, lam_spec_red) 
        flux_spec = np.ma.append(flux_spec_blue, flux_spec_red) 
        error_spec = np.ma.append(error_spec_blue, error_spec_red) 
        contam_spec = np.ma.append(contam_spec_blue, contam_spec_red) 
        zero_spec = np.append(zero_spec_blue, zero_spec_red) 

    if mask_zeros:
        ### remove bad 0th orders 
        flux_spec = np.ma.masked_where(zero_spec == 3, flux_spec)
        error_spec = newmask(error_spec, flux_spec)
        contam_spec = newmask(contam_spec, flux_spec)

    ### create an array of flags for masked regions
    masked_regions = np.zeros(lam_spec.shape, dtype=int)
    ### removed masked regions
    for k,v in config_pars.items():
        if 'mask_region' in k:
            bluecut = v[0]
            redcut = v[1]
            mask = np.logical_and(lam_spec > bluecut, lam_spec < redcut)
            flux_spec = np.ma.masked_where(mask, flux_spec)
            error_spec = newmask(error_spec, flux_spec)
            contam_spec = newmask(contam_spec, flux_spec)
            masked_regions[mask] = 1
    
    if return_masks:
        # wavelength and zeroth orders are not masked, 
        # this is fine for measure_z_interactive
        return [lam_spec, flux_spec, error_spec, contam_spec, zero_spec, masked_regions] 
    else:
        # need to apply mask to wavelength and zeroth orders so all the 
        # arrays have the same shape
        # this is necessary for the cwt code
        lam_spec = newmask(lam_spec, flux_spec)
        zero_spec = newmask(zero_spec, flux_spec)
        # compress the arrays 
        outlam = np.ma.compressed(lam_spec)
        outflux = np.ma.compressed(flux_spec)
        outerror = np.ma.compressed(error_spec)
        outcontam = np.ma.compressed(contam_spec)
        outzero = np.ma.compressed(zero_spec)
        return outlam,outflux,outerror,outcontam,outzero


