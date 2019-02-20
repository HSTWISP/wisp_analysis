#### this code makes subsamples for using with the stacking code. a host of various cuts will be included.  
### run in directory where the catalogs live. 


import astropy.io.fits as fits
import astropy.io.ascii as asciitable 
import numpy as np 
import os
import copy


########### STEP ZERO ###### CATALOGS 

####### WISP CATALOGS AND PATHS 
### read files 
masscat_wisp = asciitable.read('mederr_MARC_avgal_newparam.dat')
hdu3 = fits.open('selected_photometry_catalog.fits') 
path_wisp = '/Volumes/Thunderbay/wisps/mzr_refit/'

## extra steps for wisp cats.
wispcat = hdu3[1].data   ### this extra list is necessary because the par numbers are not in masscat_wisp 


####### 3DHST CATALOGS AND PATHS 
hdu1 = fits.open('selected_c_cami.fits')
path_3dhst = '/Volumes/Thunderbay/3DHST/mzr_refit/'


### agn matches 
xray = asciitable.read('../agn/xray_matches', format = 'no_header') 
donley = asciitable.read('../agn/donley_agn', format = 'no_header') 


#### extra steps for 3dhst cats. 
masscat_3dhst = hdu1[1].data
w=np.where(masscat_3dhst['field'] == 'GOODS-S') 
masscat_3dhst['field'][w] = 'GOODSS'
w=np.where(masscat_3dhst['field'] == 'GOODS-N') 
masscat_3dhst['field'][w] = 'GOODSN' 
fieldnames_cat = masscat_3dhst['field'] 
objid = masscat_3dhst['phot_id']
####  we're going to have to do this one field at a time. 
### gather the fields, preserve them in the order in which they appear. 
dummy = np.unique(masscat_3dhst['field'], return_index = True)
dummy2 = np.sort(dummy[1]) 
fields = masscat_3dhst['field'][dummy2]


##### STEP ONE --- remove objects from these catalogs that were removed from the sample at some stage, or rejected in the fitting. 
##############################
##############################
#### 3dhst ####################
##############################
##############################

keepers = []
z_3dhst  = [] 
foiii_3dhst = [] 
eoiii_3dhst = [] 
fhb_3dhst = [] 
ehb_3dhst = [] 
fha_3dhst = [] 
eha_3dhst = [] 
fsii_3dhst = [] 
esii_3dhst = [] 
foii_3dhst = []
eoii_3dhst = [] 
ew_oiii_obs_3dhst  = [] 


for fieldname in fields: 
    #### first define the catalog and comments files: 
    if os.path.exists(path_3dhst+  '/' + fieldname + '_output_alaina-mzr/'): 
        catalog  = path_3dhst + '/' + fieldname + '_output_alaina-mzr/' + fieldname + '_catalog_alaina-mzr.dat' 
        comments = path_3dhst + '/' + fieldname + '_output_alaina-mzr/' + fieldname + '_comments_alaina-mzr.dat' 
    elif os.path.exists(path_3dhst + '/' + fieldname + '_output_marc-mzr/'): 
        catalog  = path_3dhst + '/' + fieldname + '_output_marc-mzr/' + fieldname + '_catalog_marc-mzr.dat'  
        comments = path_3dhst + '/' + fieldname + '_output_marc-mzr/' + fieldname + '_comments_marc-mzr.dat'
    else: 
        comments = None 
        catalog = None 
        print 'Could not find fit data directory for ' + fieldname

    #### now read the catalog and the comments files for this field 
    catalog_data =asciitable.read(catalog, format = 'sextractor') 
    catalog_objid = catalog_data['ObjID'] 

    ### the comments file is not easily readable.
    f = open(comments) 
    obj_comments = [] 
    reject_list = [] 
   
    for line in f: 
       x = line.split()
       obj_comments.append(x[1])

       if len(x) > 2 :
          if x[2] == 'rejected':
                reject_list.append('rejected') 
          else: 
                reject_list.append('keep') 
       else: 
            reject_list.append('keep') 
    f.close()


    ind_field= np.where(fieldnames_cat == fieldname)
    ind_field = ind_field[0] 
    obj_field = objid[ind_field]


    for i in np.arange(len(obj_field)) : 
        if os.path.exists(path_3dhst+  '/' + fieldname + '_output_alaina-mzr/'): 
            specfile = path_3dhst + '/' + fieldname + '_output_alaina-mzr/fitdata/' + fieldname+  '_' +  '{:05d}'.format(int(obj_field[i])) + '_fitspec.dat' 
        elif os.path.exists(path_3dhst + '/' + fieldname + '_output_marc-mzr/'): 
            specfile = path_3dhst + '/' + fieldname + '_output_marc-mzr/fitdata/' + fieldname + '_' + '{:05d}'.format(int(obj_field[i])) + '_fitspec.dat'
        else:
            specfile = '~/path/to/nonsense/' 
    

        if os.path.exists(specfile) == True : 
            check1 = True
        else :
            check1 = False 
        
        if (int(obj_field[i]) in catalog_objid) and ( '{:05d}'.format(int(obj_field[i]))  in obj_comments) : 
            check2 = True 
        else :
            check2 = False

       
        obj_comments = np.array(obj_comments)
        w=np.where(obj_comments == '{:05d}'.format(int(obj_field[i]))   )

        if np.size(w[0]) == 0 :
            check3 = False 
        else: 
            if reject_list[w[0][0]] == 'keep' :
                check3 = True 
            else:
                check3 = False 

        #print fieldname, '{:05d}'.format(int(obj_field[i])), check1, check2, check3 
                
        if (check1 == True) & (check2 == True) & (check3 == True) : 
            keepers.append(1) 
            w=np.where(catalog_objid == int(obj_field[i])) 
            z_3dhst.append(catalog_data['redshift'][w[0][-1]])  ###  just in case the object is in the catlag more than once, take the last one. 
            foiii_3dhst.append(catalog_data['oiii_flux'][w[0][-1]]) 
            eoiii_3dhst.append(catalog_data['oiii_err'][w[0][-1]]) 
            fha_3dhst.append(catalog_data['hanii_flux'][w[0][-1]])
            eha_3dhst.append(catalog_data['hanii_err'][w[0][-1]]) 
            fhb_3dhst.append(catalog_data['hb_flux'][w[0][-1]]) 
            ehb_3dhst.append(catalog_data['hb_err'][w[0][-1]]) 
            fsii_3dhst.append(catalog_data['sii_flux'][w[0][-1]])
            esii_3dhst.append(catalog_data['sii_err'][w[0][-1]])
            foii_3dhst.append(catalog_data['oii_flux'][w[0][-1]]) 
            eoii_3dhst.append(catalog_data['oii_err'][w[0][-1]])
            ew_oiii_obs_3dhst.append(catalog_data['oiii_EW_obs'][w[0][-1]])

    

        else: 
            keepers.append(0) 



#### LASTLY, cull the parent catalog. 

print 'Size of Parent 3dhst catalog: '
print np.size(masscat_3dhst['field']) 

print 'Size of retained 3dhst catalog' 
print np.sum(keepers)


keepers = np.array(keepers) 
w=np.where(keepers == 1) 
#### gather what we need from these catalogs: 
id_3dhst = masscat_3dhst['phot_id'][w] 
field_3dhst = masscat_3dhst['field'][w] 
logm_3dhst = masscat_3dhst['logM_50_cami'][w] 


z_3dhst = np.array(z_3dhst) 
foiii_3dhst = np.array(foiii_3dhst) 
eoiii_3dhst = np.array(eoiii_3dhst) 

print np.size(z_3dhst) 
print np.size(id_3dhst) 

##############################
##############################
#### WISP ####################
##############################
##############################
objid_parent = masscat_wisp['ID']   ####  integer numpy array 
par_parent = wispcat['par']


keepers = []
z_wisp = [] 
foiii_wisp = []
eoiii_wisp = []
fha_wisp = []
eha_wisp = []
fhb_wisp = [] 
ehb_wisp = []
fsii_wisp = [] 
esii_wisp = []
foii_wisp = [] 
eoii_wisp  = []

ew_oiii_obs_wisp = [] 





for i  in np.arange(len(par_parent)): 

    if os.path.exists(path_wisp +  '/Par' + str(par_parent[i]) + '_output_alaina-mzr/') : 
            specfile = path_wisp + '/Par' + str(par_parent[i]) + '_output_alaina-mzr/fitdata/Par' + str(par_parent[i]) + '_BEAM_' + str(objid_parent[i]) + '_fitspec.dat' 
            catalog  = path_wisp  + '/Par' + str(par_parent[i])+ '_output_alaina-mzr/Par' + str(par_parent[i]) + 'list_catalog_alaina-mzr.dat' 
            comments = path_wisp  + '/Par' + str(par_parent[i])+ '_output_alaina-mzr/Par' + str(par_parent[i]) + 'list_comments_alaina-mzr.dat'

    elif  os.path.exists(path_wisp +  '/Par' + str(par_parent[i]) + '_output_marc-mzr/') : 
            specfile = path_wisp + '/Par' + str(par_parent[i]) + '_output_marc-mzr/fitdata/Par' + str(par_parent[i]) + '_BEAM_' + str(objid_parent[i]) + '_fitspec.dat' 
            catalog  = path_wisp  + '/Par' + str(par_parent[i])+ '_output_marc-mzr/Par' + str(par_parent[i]) + 'list_catalog_marc-mzr.dat' 
            comments = path_wisp  + '/Par' + str(par_parent[i])+ '_output_marc-mzr/Par' + str(par_parent[i]) + 'list_comments_marc-mzr.dat'  
           
    else:

        specfile = '~/path/to/nonsense/'


    ### check 1, does the spectrum exist 
    if os.path.exists(specfile) == True : 
            check1 = True
    else :
            check1 = False 


    #### check 2, is the object in the comments file and in the catalog?

    if (os.path.exists(catalog)) & (os.path.exists(comments)): 
    
        ## read catalog data. 
        catalog_data =asciitable.read(catalog, format = 'sextractor') 
        catalog_objid = catalog_data['ObjID'] 
    
        ### read comments file, get list of keeps vs. rejects 
        f = open(comments) 
        obj_comments = [] 
        reject_list = [] 
   
        for line in f: 
           x = line.split()
           obj_comments.append(int(x[1]))

           if len(x) > 2 :
              if x[2] == 'rejected':
                    reject_list.append('rejected') 
              else: 
                    reject_list.append('keep') 
           else: 
                reject_list.append('keep') 
        f.close()

    else: 
         check2 = False 
         check3 = False



    
    ####  is objid_parent in obj_comments and catalog_objid 
    if (objid_parent[i] in  obj_comments) and (objid_parent[i] in catalog_objid): 
        check2 = True 
    else: 
        check2 = False 


    obj_comments = np.array(obj_comments)
    w=np.where(obj_comments == objid_parent[i])
    
    if np.size(w[0]) == 0 :
            check3 = False 
    else: 
        if reject_list[w[0][0]] == 'keep' :
            check3 = True 
        else:
            check3 = False 
     
    #print par_parent[i], objid_parent[i], check1, check2, check3 
                
    if (check1 == True) & (check2 == True) & (check3 == True) : 
        keepers.append(1)  
        #### GATHER INFORMATION FOR KEEPERS
        w=np.where(catalog_objid == objid_parent[i])
        z_wisp.append(catalog_data['redshift'][w[0][-1]])  ###  just in case the object is in the catlag more than once, take the last one. 
        foiii_wisp.append(catalog_data['oiii_flux'][w[0][-1]]) 
        eoiii_wisp.append(catalog_data['oiii_err'][w[0][-1]])
        fha_wisp.append(catalog_data['hanii_flux'][w[0][-1]])
        eha_wisp.append(catalog_data['hanii_err'][w[0][-1]]) 
        fhb_wisp.append(catalog_data['hb_flux'][w[0][-1]]) 
        ehb_wisp.append(catalog_data['hb_err'][w[0][-1]]) 
        fsii_wisp.append(catalog_data['sii_flux'][w[0][-1]])
        esii_wisp.append(catalog_data['sii_err'][w[0][-1]])
        foii_wisp.append(catalog_data['oii_flux'][w[0][-1]]) 
        eoii_wisp.append(catalog_data['oii_err'][w[0][-1]])
        ew_oiii_obs_wisp.append(catalog_data['oiii_EW_obs'][w[0][-1]])



    else: 
        keepers.append(0) 
    

z_wisp  = np.array(z_wisp)
foiii_wisp = np.array(foiii_wisp) 
eoiii_wisp  = np.array(eoiii_wisp) 
keepers = np.array(keepers) 

print 'Size of Parent WISP catalog' 
print np.size(masscat_wisp['ID'])


print 'Size of Retained WISP catalog' 
print np.sum(keepers) 



w=np.where(keepers == 1) 

id_wisp = masscat_wisp['ID'][w]
par_wisp = wispcat['par'][w] 
logm_wisp = masscat_wisp['logM_50'][w] 


print np.size(z_wisp), np.size(id_wisp), np.size(par_wisp) 


#### let's make an output/master catalog of these things, so I can play around and decide on binning. 

field = []

for par in par_wisp:
    field.append('Par' + str(par))
for ff  in field_3dhst : 
    field.append(ff)

field = np.array(field) 
        
objid = np.append(id_wisp, id_3dhst) 
logm = np.append(logm_wisp, logm_3dhst) 
z = np.append(z_wisp, z_3dhst) 
foiii = np.append(foiii_wisp, foiii_3dhst) 
eoiii = np.append(eoiii_wisp, eoiii_3dhst) 
fha = np.append(fha_wisp, fha_3dhst) 
eha = np.append(eha_wisp, eha_3dhst) 
fhb = np.append(fhb_wisp, fhb_3dhst) 
ehb = np.append(ehb_wisp, ehb_3dhst) 
fsii = np.append(fsii_wisp, fsii_3dhst) 
esii = np.append(esii_wisp, esii_3dhst) 
foii = np.append(foii_wisp, foii_3dhst) 
eoii  = np.append(eoii_wisp, eoii_3dhst) 

ew_oiii_obs = np.append(ew_oiii_obs_wisp, ew_oiii_obs_3dhst) 



### remove a few more bad objects with zero error or foiii = -1 
w=np.where( (foiii > 0) & (eoiii > 0) & (z > 1.2))  

field = field[w] 
objid = objid[w] 
logm = logm[w] 
z = z[w] 
foiii = foiii[w] 
eoiii = eoiii[w]
fha = fha[w] 
eha = eha[w] 
fhb =fhb[w]
ehb = ehb[w]
fsii = fsii[w]
esii = esii[w]
foii = foii[w]
eoii = eoii[w]
ew_oiii_obs = ew_oiii_obs[w] 


print 'number of objects after OIII no cov,  OIII err = 0, z<1.2  removed'  
print np.size(field), np.size(objid), np.size(logm), np.size(z), np.size(foiii), np.size(eoiii)

objid = objid.astype(int)

#mex_agn_selection
ehb2 = copy.deepcopy(ehb) 
w=np.where(fhb == 0 ) 
ehb2[w] = eoiii[w] 


mex_flag = np.zeros(len(z))
hbsnr = fhb / ehb2 
oiiihb = np.log10(foiii/fhb) 
w=np.where(hbsnr <=3) 
oiiihb[w] =  np.log10(foiii[w] / (3 * ehb2[w])) 

x =logm 
x2 = x-0.75 
ymex_lower_obj = 0.375 / (x2 - 10.5) + 1.14  + 0.11
w=np.where(x2 > 9.6)
ymex_lower_obj[w] = 352.066 - 93.8249 * x2[w] + 8.32651 * x2[w]**2 - 0.246416 * x2[w]**3   + 0.11 
agn_flag = np.where((oiiihb > ymex_lower_obj) & (logm >0)) 
mex_flag[agn_flag] = 1 


xray_agn_flag = np.zeros(len(z)) 

xray_field = xray['col1'] 
xray_id = xray['col2'] 

for i in np.arange(len(xray_id)): 
    w=np.where( (field == xray_field[i]) & (objid == xray_id[i])) 
    if len(w[0]) == 1: 
        xray_agn_flag[w] = 1. 

ir_agn_flag = np.zeros(len(z)) 
ir_field = donley['col1'] 
ir_id = donley['col2'] 


for i in np.arange(len(ir_id)):
    w=np.where((field == ir_field[i]) & (objid == ir_id[i]))
    if len(w[0]) == 1:
        ir_agn_flag[w] = 1.




data = [field, objid, logm, z, foiii, eoiii, ew_oiii_obs, fhb, ehb, fha, eha, fsii, esii, foii, eoii, mex_flag, xray_agn_flag, ir_agn_flag] 
colnames = ['Field', 'ID', 'logM', 'z', 'foiii', 'eoiii', 'EW_oiii_obs', 'fhb', 'ehb', 'fhanii', 'ehanii', 'fsii', 'esii', 'foii', 'eoii', 'MeX_AGN', 'Xray', 'IR_AGN'] 

asciitable.write(data, 'master_stacklist.v013019', names = colnames, overwrite = True) 







    




    
    

        







