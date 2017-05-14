
import astropy.io.fits as fits
import numpy as np
import os 
from wisp_analysis import *
from distutils.sysconfig import *


def show2dNEW (grism,parno,obid,zeroarr,user,trans,zran1=-0.2,zran2=0.75):
# In version 1.0, will first look for wavelength-calibrated stamps in the G1??_DRIZZLE directories; failing this, will default to old stamps
    # zero and first order positions
#    firstx = firstarr['x']
#    firsty = firstarr['y']
#    firstlen = firstarr['len']
#    firstwid = firstarr['width']
#    firstid = firstarr['objid']
    zerox = zeroarr['x']
    zeroy = zeroarr['y']
    zeroid = zeroarr['objid']

    dims=()
    zrad=10.0
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'
    path2dl=par_root_dir + grism + '_DRIZZLE/aXeWFC3_' +grism + '_mef_ID'+str(obid)+'.fits'
    if os.path.exists(path2dl)==1:
        path2d=path2dl
    else:
        path2d=par_root_dir+'Stamps/Par'+ str(parno)+'_'+grism+'_BEAM_'+str(obid)+'A.fits'

    if grism=='G102':
        frameno='1'
    elif grism=='G141':
        frameno='2'
    if os.path.exists(path2d)==1:
        infits=fits.open(path2d)
        ### changing to read in 1st data extension ###
        #darr=infits[-1].data
        hdr=infits[1].header
        darr=infits[1].data
        dims=darr.shape
        infits.close()
    elif os.path.exists(path2d)==0:
        print "%s stamp not found." % (grism)
        return False

    ### USING THE DRIZZLE TRANSFORMATIONS TO GET ZEROTH ORDERS ###
    _cx = np.array([xcoo for xcoo in zerox]) - hdr['BB0X'] - 1
    _cy = np.array([ycoo for ycoo in zeroy]) - hdr['BB0Y'] - 1
    cx = hdr['D001OUXC'] + (hdr['DRZ00'] + hdr['DRZ01']*(_cx-hdr['D001INXC']) + hdr['DRZ02']*(_cy-hdr['D001INYC']))
    cy = hdr['D001OUYC'] + (hdr['DRZ10'] + hdr['DRZ11']*(_cx-hdr['D001INXC']) + hdr['DRZ12']*(_cy-hdr['D001INYC']))
    # convert to (Angs,arcsec) coords
    cx = (cx - hdr['CRPIX1'])*hdr['CDELT1'] + hdr['CRVAL1']
    cy = (cy - hdr['CRPIX2'])*hdr['CDELT2'] + hdr['CRVAL2']
    rad = 5 * hdr['CDELT1']
    outcoo=par_root_dir+"Spectra/temp_zero_coords_%s.reg"%user
    if os.path.exists(outcoo)==1:
        os.unlink(outcoo)
    f = open(outcoo, 'w')
    f.write('wcs;\n')
    for j in range(len(zerox)):
        f.write('circle(%.2f,%.4f,%.1f) # color=red text={%s}\n' % (cx[j],cy[j],rad,zeroid[j]))
    f.close()


    ### THIS WAS ALL FOR THE OLD TRANSFORMATION ###
#    matchind=0
#    i=0
#    for fid in firstid:
#        if fid==obid:
#            matchind=i
#            break
#        i=i+1
#    xmin=firstx[matchind]-firstlen[matchind]/2.0
#    xmax=firstx[matchind]+firstlen[matchind]/2.0
#    ymin=firsty[matchind]-firstwid[matchind]/2.0
#    ymax=firsty[matchind]+firstwid[matchind]/2.0
#    
#    numzer=0
#    if len(dims)>0:
#        xdim=float(max(dims))
#        ydim=float(min(dims))
#    outcoo=par_root_dir+"Spectra/temp_zero_coords.reg"
#    if os.path.exists(outcoo)==1:
#        os.unlink(outcoo)
#    coordout=open(outcoo,'w')
#    print >>coordout, "image"
#    for j in range(len(zerox)):
#        # MB: use dims of stamp, rather than size of 1st order in region file
#        if zerox[j] >= (firstx[matchind] - xdim/2.) and \
#           zerox[j] <= (firstx[matchind] + xdim/2.) and \
#           zeroy[j] >= (firsty[matchind] - ydim/2.) and \
#           zeroy[j] <= (firsty[matchind] + ydim/2.) and grism=='G102':
##    	if zerox[j]>=(xmin-zrad) and zerox[j]<=(xmax+zrad) and zeroy[j]>=(ymin-zrad) and zeroy[j]<=(ymax+zrad) and grism=='G102':
#    		print >>coordout, "circle(%.2f,%.2f,5.0) # text={%s}" % (zerox[j]/1.7-firstx[matchind]/1.7+212./2.0-13,zeroy[j]/1.6-firsty[matchind]/1.6+ydim/2.0+3.6,zeroid[j])
#
#        elif zerox[j] >= (firstx[matchind] - xdim/2.) and \
#             zerox[j] <= (firstx[matchind] + xdim/2.) and \
#             zeroy[j] >= (firsty[matchind] - ydim/2.) and \
#             zeroy[j] <= (firsty[matchind] + ydim/2.) and grism=='G141':
##    	elif  zerox[j]>=(xmin-zrad) and zerox[j]<=(xmax+zrad) and zeroy[j]>=(ymin-zrad) and zeroy[j]<=(ymax+zrad) and grism=='G141':
#    		print >>coordout, "circle(%.2f,%.2f,5.0) # text={%s}" % (zerox[j]/1.7-firstx[matchind]/1.7+184./2.0,zeroy[j]/1.6-firsty[matchind]/1.6+ydim/2.0+0.6,zeroid[j])
#    numzer=numzer+1
#    coordout.close()
    
    if trans=='log':
        zscale='log'
    else:
        zscale='linear'
    cmd='xpaset -p ds9 frame '+frameno
    os.system(cmd)
    cmd='xpaset -p ds9 file '+path2d
    os.system(cmd)
    cmd='xpaset -p ds9 scale limits '+str(zran1)+' '+str(zran2)
    os.system(cmd)
    cmd='xpaset -p ds9 scale '+zscale
    os.system(cmd)
    #cmd='xpaset -p ds9 zoom to fit'
    #os.system(cmd)
    cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/temp_zero_coords_%s.reg'%user
    os.system(cmd)
    # MR
    if frameno=='1':
        cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/G102_trace.reg'
    if frameno=='2':
        cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/G141_trace.reg'
    os.system(cmd)



def showDirectNEW(obid,load_image=False):
    """
    Removed lineno, which was only used to check whether the images 
    should be reloaded.
    """
    obid=int(obid)
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'

    path2direct=par_root_dir+'DATA/DIRECT_GRISM/'
    path110=path2direct+'F110W_drz.fits'
    path140=path2direct+'F140W_drz.fits'
    path160=path2direct+'F160W_drz.fits'
    path140cat=path2direct+'fin_F140.cat'
    path160cat=path2direct+'fin_F160.cat'

    if os.path.exists(path110)==0 and os.path.exists(path140)==0 and os.path.exists(path160)==0:
        print "No Direct Images Found."
        return 0
    if os.path.exists(path140cat)==1:
        infHcat=open(path140cat,'r')
    elif os.path.exists(path160cat)==1:
        infHcat=open(path160cat,'r')
    else:
        return 0
    xcen,ycen=-1,-1
    for line in infHcat:
        if line[0]!='#':
            entries=line.split()
            if int(entries[1])==obid:
                xcenter,ycenter=float(entries[7]),float(entries[8])
                hexcoo=[entries[7],entries[8]]
    infHcat.close()

    # load the direct images
    if load_image:
        if os.path.exists(path110)==1:
            cmd='xpaset -p ds9 frame 3'
            os.system(cmd)
            cmd='xpaset -p ds9 file '+path110
            os.system(cmd)
            ### using F110_drz.reg with F110W_drz.fits
            cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F110_drz.reg'
            os.system(cmd)
            cmd='xpaset -p ds9 pan to '+hexcoo[0]+' '+hexcoo[1]+' fk5'
            os.system(cmd)
        if os.path.exists(path140)==1:
            cmd='xpaset -p ds9 frame 4'
            os.system(cmd)
            cmd='xpaset -p ds9 file '+path140
            os.system(cmd)
            cmd='xpaset -p ds9 regions file '+par_root_dir+ 'DATA/DIRECT_GRISM/F140_drz.reg'
            os.system(cmd)
            cmd='xpaset -p ds9 pan to '+hexcoo[0]+' '+hexcoo[1]+' fk5'
            os.system(cmd)
        elif os.path.exists(path160)==1:
            cmd='xpaset -p ds9 frame 4'
            os.system(cmd)
            cmd='xpaset -p ds9 file '+path160
            os.system(cmd)
            cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F160_drz.reg'
            os.system(cmd)
            cmd='xpaset -p ds9 pan to '+hexcoo[0]+' '+hexcoo[1]+' fk5'
            os.system(cmd)
    # pan to the coordinates of this object
    if os.path.exists(path110):
        panDirect(hexcoo[0],hexcoo[1],grism='G102')
    panDirect(hexcoo[0],hexcoo[1])


def showDispersed(obid,load_image=False):  # MB
    """
    Removed lineno, which was only used to check whether the images 
    should be reloaded.
    """
    obid=int(obid)
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'

    path2dispersed=par_root_dir+'DATA/DIRECT_GRISM/'
    ### Using G102.fits instead of G102_drz.fits ###
    path102=path2dispersed+'G102.fits'
    path141=path2dispersed+'G141.fits'
    path102_0reg = os.path.join(path2dispersed, 'G102_0th.reg')
    path102_1reg = os.path.join(path2dispersed, 'G102_1st.reg')
    path141_0reg = os.path.join(path2dispersed, 'G141_0th.reg')
    path141_1reg = os.path.join(path2dispersed, 'G141_1st.reg')
    if os.path.exists(path102)==0 and os.path.exists(path141)==0:
        print "No Grism Images Found."
        return 0
    # get center of 1st order
    ### Using same syntax as getzeroorders for consistency ###
    if os.path.exists(path102_1reg)==1:
        reg102=open(path102_1reg,'r')
        x102,y102=-1,-1
        for line in reg102:
            # using same syntax as getzeroorders
            linesplit = line.split()
            textid = int(re.search('\d+',linesplit[-1]).group(0))
            if textid==obid:
                x102 = float(linesplit[1].split(',')[0])
                y102 = float(linesplit[2].split(',')[0])
        reg102.close()
    if os.path.exists(path141_1reg)==1:
        reg141=open(path141_1reg,'r')
        x141,y141=-1,-1
        for line in reg141:
            linesplit = line.split()
            textid = int(re.search('\d+',linesplit[-1]).group(0))
            if textid==obid:
                x141 = float(linesplit[1].split(',')[0])
                y141 = float(linesplit[2].split(',')[0])
        reg141.close()
    
    if load_image:
        if os.path.exists(path102)==1:
            cmd='xpaset -p ds9 frame 5'
            os.system(cmd)
            cmd='xpaset -p ds9 file '+path102
            os.system(cmd)
            cmd='xpaset -p ds9 regions file '+path102_0reg
            os.system(cmd)
            cmd='xpaset -p ds9 regions file '+path102_1reg
            os.system(cmd)
            cmd='xpaset -p ds9 pan to %f %f image' % (x102,y102)
            os.system(cmd)
        if os.path.exists(path141)==1:
            cmd='xpaset -p ds9 frame 6'
            os.system(cmd)
            cmd='xpaset -p ds9 file '+path141
            os.system(cmd)
            cmd='xpaset -p ds9 regions file '+path141_0reg
            os.system(cmd)
            cmd='xpaset -p ds9 regions file '+path141_1reg
            os.system(cmd)
            cmd='xpaset -p ds9 pan to %f %f image' % (x141, y141)
            os.system(cmd)

    # pan to the coordinates of this object
    if os.path.exists(path102):
        panDispersed(x102,y102,grism='G102')
    panDispersed(x141,y141)


def createAltGrismRegion(grism):
    workingdir=os.getcwd()
    par = os.path.dirname(workingdir)
    cat110 = os.path.join(par, 'DATA/DIRECT_GRISM/fin_F110.cat')
    cat140 = os.path.join(par, 'DATA/DIRECT_GRISM/fin_F140.cat')
    cat160 = os.path.join(par, 'DATA/DIRECT_GRISM/fin_F160.cat')

    if grism == 'G102':
        if os.path.exists(cat110) == 1:
            cat = np.genfromtxt(cat110, dtype=[('num',int),('a_img',float),
                               ('mag',float)], usecols=(1,4,12))
        else:
            print cat110
            return 0
    if grism == 'G141':
        if os.path.exists(cat140) == 1:
            cat = np.genfromtxt(cat140, dtype=[('num',int),('a_img',float),
                                 ('mag',float)], usecols=(1,4,12))
        elif os.path.exists(cat160) == 1:
            cat = np.genfromtxt(cat160, dtype=[('num',int),('a_img',float),
                                 ('mag',float)], usecols=(1,4,12))
        else:
            print 'nope2'
            return 0

    f = open(os.path.join(workingdir,grism+'_temp.reg'), 'w')
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('image\n')
    cenx = np.zeros(cat['num'].shape[0], dtype=float)
    ceny = np.zeros(cat['num'].shape[0], dtype=float)
    for i in range(cat['num'].shape[0]):
        objid = cat['num'][i]
        gfile = os.path.join(par,'%s_DRIZZLE/aXeWFC3_%s_mef_ID%i.fits'%(grism,grism,objid))
        if os.path.isfile(gfile):
            hdr = fits.getheader(gfile, 1)
            # coords of bounding box
            boxx = np.array([hdr['bb0x'], hdr['bb1x']])
            boxy = np.array([hdr['bb0y'], hdr['bb1y']])
            cenx[i] = boxx[0] + (boxx[1] - boxx[0]) / 2.
            ceny[i] = boxy[0] + (boxy[1] - boxy[0]) / 2.
            slitwidth = hdr['slitwidt']
            sx = 184
            sy = 8
#            sy = slitwidth * cat['a_img'][i]
            f.write('box(%f,%f,%i,%i,0) # text={%i}\n' % (cenx[i],ceny[i],sx,sy,objid))
    f.close()
    return cenx,ceny


def panDirect(ra,dec,grism='G141'):
    # Pan to coords in frame
    if grism=='G141':
        fno='4'
    else:
        fno='3'
    cmd='xpaset -p ds9 frame ' + fno
    os.system(cmd)
    cmd='xpaset -p ds9 pan to '+ra+' '+dec+' fk5'
    os.system(cmd)


def panDispersed(xx,yy,grism='G141'):  # MB
    # Pan to coords in frame
    if grism=='G141':
        fno='6'
    else:
        fno='5'
    cmd='xpaset -p ds9 frame ' + fno
    os.system(cmd)
    cmd='xpaset -p ds9 pan to %f %f image' % (xx,yy)
    os.system(cmd)


def reloadReg():
    workingdir=os.getcwd()
    dirpts=workingdir.split('/')[1:-1]
    par_root_dir='/'
    for pdir in dirpts:
        par_root_dir= par_root_dir +pdir + '/'
    
    # reload direct image region files
    if os.path.exists(par_root_dir+'DATA/DIRECT_GRISM/F110_drz.reg')==1:
        cmd='xpaset -p ds9 frame 3'
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F110_drz.reg'
        os.system(cmd)
    if os.path.exists(par_root_dir+'DATA/DIRECT_GRISM/F160_drz.reg')==1:
        cmd='xpaset -p ds9 frame 4'
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F160_drz.reg'
        os.system(cmd)
    elif os.path.exists(par_root_dir+'DATA/DIRECT_GRISM/F140_drz.reg')==1:
        cmd='xpaset -p ds9 frame 4'
        os.system(cmd)
        cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F140_drz.reg'
        os.system(cmd)


