from wisp_analysis import * 
import glob
import astropy.io.ascii as asciitable



def gather_secure_sample(): 

    catalog_files = glob.glob('Par*/Par*list_catalog*.dat')
    #comment_files = glob.glob('Par*/Par*list_comments*.dat') 


    out = open('wisp_secure', 'w') 

    for cat in catalog_files: 
        z_secure = None
        table = asciitable.read(cat) 
        oiii_snr = table['oiii_flux'] /table['oiii_err'] 
        oii_snr = table['oii_flux'] /table['oii_err']
        hb_snr = table['hb_flux'] / table['hb_err'] 
        hg_snr = table['hg_flux'] / table['hg_err'] 
        ha_snr = table['hanii_flux']/table['hanii_err']
        sii_snr = table['sii_flux'] / table['sii_err'] 
        obj = table['ObjID'] 
        parno  = table['ParNo'] 
        print parno[0]
        z_secure = []
        for i in np.arange(len(obj)): 
            snr_supporting = np.array( [oii_snr[i], hb_snr[i], hg_snr[i], ha_snr[i], sii_snr[i]] )
            w=np.where(snr_supporting >= 3)
            w2 = np.where(snr_supporting >= 2)
            if np.size(w[0]) >=1 : 
                z_secure.append(1) 
            elif np.size(w2[0]) >= 2: 
                z_secure.append(1)
            else: 
                z_secure.append(0) 
        z_secure = np.array(z_secure)
        sample = np.where((oiii_snr >= 5) & (z_secure ==1)) 
        if np.size(sample[0]) > 0: 
            for i in np.arange(np.size(sample[0])):
                out.write('Par' + str(parno[0]) + '   '  + str(obj[sample[0]][i]) + '\n')


    out.close()



