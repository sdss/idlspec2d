#!/usr/bin/env python

"""
Reformat spectra into a single fits file per object, combining all necessary
pieces from spPlate, spCFrame, spFrame, spFlat, and spZbest

Stephen Bailey, Summer 2011

Bugs/Features:
  - Assumes RUN1D=RUN2D, and gets them from spAll (not env vars), assuming
    that spAll has one and only one RUN2D.

Julien Guy, Sept. 2015 :
 add option --tpcorr to include spectro-photometric corrections from fibers positioned at LAMBDA_EFF=4000A
"""

import sys
import os
import os.path
from glob import glob           #- File pattern globbing (spFrame*.fits)
import re                       #- Regular expressions
from time import asctime
import numpy as N
import pyfits
from astropy.io import fits
import h5py                     #- Needed to read spectro-photometric calibration correction file

def write_readme(filename, header=None, allexp=True):
    """
    Write a file explaining the directory structure and file format.
    
    Optional inputs:
        header : prepend this header string
        allexp : if True, also describe per-exposure HDUs
    """
    import time
    fx = open(filename, 'w')
    print >> fx, "Spectra repackaging generated on " + time.asctime()
    if header is not None:
        print >> fx
        print >> fx, header

    print >> fx, """
For each object, there is one file per plugging (plate-mjd-fiber), 
containing both the coadded spectrum and optionally the individual exposure
frames.  This groups the information from spFrame, spCFrame, spFlat, spPlate,
spZbest, spZline, and spAll so that for each object you only need to read
one file.

    HDU 0  : Header info
    HDU 1  : Coadded spectrum
    HDU 2  : Summary metadata copied from spAll
    HDU 3  : Line fitting metadata from spZline
    HDU 4+ : [Optional] Individual frame spectra [B, R for each exposure]

These are grouped in subdirectories by plate, e.g.:

    README.txt     : this file
    spAll-XXX.fits : subset of spAll*.fits (e.g. qso, star, gal)
    4080/          : dir for objects on plate 4080
        spec-4080-55368-0487.fits  : spec-PLATE-MJD-FIBER
        spec-4080-55471-0485.fits  : different plugging, same plate
    ...

Note that the same THING_ID may have been observed on more than
one plate.  Use the THING_ID, PLATE, MJD, FIBERID metadata in
spAll to find the files you need.

The format of each spec file is:

HDU 0 :
    Header : from input spPlate with additional keywords from spAll:
        
        Keyword     spAllColumn   Comment
        -------     -----------   -------
    [*] PLUG_RA     PLUG_RA       RA of object [deg]
    [*] PLUG_DEC    PLUG_DEC      dec of object [deg]
        THING_ID    THING_ID      Unique object identifier
        FIBERID     FIBERID       Fiber number (1-1000)
       
    [*] Note that RA and DEC already exist in the spPlate headers
        but they are the telescope boresite RA and DEC, not the
        RA and DEC of the object. 

    Keywords that spPlate inherited from a single exposure have
    been removed from the primary HDU.  Other keywords like NEXP
    have been changed to reflect the actual exposures which went
    into this file (e.g. not including the exposures for fibers
    on a different spectrograph)
        
    Data   : None

HDU 1 : Coadded Spectrum
    Header : Minimum to define table
    Data : Binary table with columns taken from spPlate and spZbest:
        flux      : flux in 10^-17 ergs/s/cm^2/A
        loglam    : log10(lambda[A])
        ivar      : inverse variance
        and_mask  : mask bits which affect every spectrum in coadd
        or_mask   : mask bits which affect at least one spectrum in coadd
        wdisp     : wavelength dispersion in dloglam units
        sky       : subtracted sky flux in 10^-17 ergs/s/cm^2/A
        model     : best fit model for classification & redshift (from spZbest)
        calib_corr : (optional) : spectro-photometric correction for fibers positioned at LAMBDA_EFF=4000

HDU 2 : Copy of row for this object from spAll table
    http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/spAll.html

HDU 3 : Copy of rows for this object from spZline table
    http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/RUN1D/spZline.html"""

    #- Continue with explanation of per-exposure HDUs if included
    if allexp:
        print >> fx, """
HDU 4 .. n+4 : Individual frames.
    For each exposure there is one HDU for the red camera and one for the blue.
    These are in the order of the EXPIDnnn keywords in the HDU 0 header, and
    their EXTNAME keywords match the EXPIDnnn keywords from HDU 0.

    Header: Taken from HDU0 of individual spCFrame files
    Keyword YPIX0 defines the y-pixel location of the first flux bin.

    Data: Binary table with columns taken from spCFrame files:
        flux
        loglam
        ivar
        mask
        wdisp
        sky
        calib : conversion between flux and electrons: flux = electrons*calib
        x     : x-pixel location of trace on CCD

How to convert flux, sky, and ivar back to extracted photons (electrons):

    obj_photons  = flux / calib
    sky_photons  = sky / calib
    ivar_photons = ivar * calib^2
      (includes variance from all sources: object, sky, read noise, etc)

You can approximately estimate the non-photon variance with:

    photons = (flux + sky) / calib
    var_photons = 1.0 / (ivar * calib^2)
    var_extra_photons = var_photons - photons
    var_extra_flux = var_extra_photons * calib

"""
    else:
        print >> fx, "\nIndividual exposures are not included in these files"
        
    fx.close()

class CFrame(object):
    """
    Convenience wrapper class for spCFrame, spFrame, and spFlat data
    
    Derives spFrame dir/name from spCFrame, and assumes that it is gzipped.
    """
    def __init__(self, cframefile):

        #- Load original framefile and find out original dimensions
        ### print cframefile
        framefile = cframefile.replace('spCFrame', 'spFrame') + '.gz'
        eflux = fits.getdata(framefile, 0)
        nfiber, npix = eflux.shape
                
        #- Load spCFrame file; trim arrays back to original size
        fx = fits.open(cframefile)
        self.flux = fx[0].data[:, 0:npix]
        self.ivar = fx[1].data[:, 0:npix]
        self.mask = fx[2].data[:, 0:npix]
        self.loglam = fx[3].data[:, 0:npix]
        self.wdisp  = fx[4].data[:, 0:npix]
        self.sky    = fx[6].data[:, 0:npix]
        self.x      = fx[7].data[:, 0:npix]
        self.header = fx[0].header
        
        #- Load superflat spCFrame[8] and fiberflat spFlat[0]
        filedir, basename = os.path.split(cframefile)
        superflat = fx[8].data[:, 0:npix]
        flatfile = fx[0].header['FLATFILE'].replace('sdR', 'spFlat')
        flatfile = flatfile.replace('.fit', '.fits.gz')
        flatfile = os.path.join(filedir, flatfile)
        fiberflat = fits.getdata(flatfile, 0)
        
        #- Calculate calibration vector: flux = electrons * calib
        electrons = eflux * fiberflat * superflat
        ii = N.where(electrons != 0.0)
        self.calib = N.zeros(self.flux.shape)
        self.calib[ii] = self.flux[ii] / electrons[ii]
                
        fx.close()

def load_spCFrame_files(platedir):
    """
    Load all spCFrame files in a given directory.
    Return a dictionary of CFrame objects, keyed by camera-expid string
    """
    print "loading spCFrame files from " + platedir
    cframes = dict()
    for filename in glob(os.path.join(platedir, 'spCFrame-*.fits')):
        print '   ', os.path.basename(filename), asctime()
        expid = get_expid(filename)
        cframes[expid] = CFrame(filename)

    return cframes

def good_ivar(ivar, fiber=None):
    """
    return indices of for array from first non-zero ivar to the last
    non-zero ivar.  i.e. trim off leading and trailing contiguously zero ivars.
    
    If all ivar==0, return indices for full array since code might get
    grumpy with completely blank arrays.
    """
    if N.all(ivar == 0):
        if fiber is not None:
            print 'WARNING: All ivar==0 for fiber', fiber
        else:
            print 'WARNING: All ivar==0'
        return N.arange(len(ivar))

    ivar_good = N.where(ivar > 0.0)[0]
    return N.arange(ivar_good[0], ivar_good[-1]+1)

def get_expid(filename):
     """parse /path/to/spBlat-b1-00123456.fits.gz into b1-00123456"""
     try:
         return re.search('-([br][12]-\d{8}).fits', filename).group(1)
     except AttributeError:  #- search failed
         return None
        
def plate_to_string(plate):

    if plate<10000:
        return "%04d"%plate
    else:
        return "%d"%plate

def process_plate(datadir, outdir, plate, mjd, fibers, spAll, allexp=True, tpcorr_h5=None, plate_f=False, legacy=False):
    """
    Process a plate's worth of objects
    
    Inputs
        datadir : input base directory, e.g. $BOSS_SPECTRO_REDUX/v5_4_40/
        outdir  : output base directory
        plate   : plate to process
        mjd     : mjd to process
        fibers  : *array* of fibers on this plate-mjd to process
        spAll   : metadata from spAll which includes these spectra
        allexp  : if False, don't write individual exposures (default True)
        tpcorr_h5 : if not None, add a column with spectro-photometric correction (not 1 for fibers positionned at LAMBDA_EFF=4000)

    Outputs:
        writes files to outdir/plate/spec-plate-mjd-fiber.fits
    """
    #- Load all C/Frame files for this plate
    platestr = plate_to_string(plate)
    if plate_f or legacy:
        platedir = '%s/%s/' % (datadir, platestr+'p')
    else:
        platedir = '%s/%s/' % (datadir, platestr)
    if allexp:
        cframes = load_spCFrame_files(platedir)

    #- Open spPlate, spZbest, and spZline files
    #spPlateFile = '%s/spPlate-%s-%d.fits' % (platedir, platestr, mjd)
    spPlateFile = '%s/spField-%s-%d.fits' % (platedir, platestr, mjd)
    print 'Processing', os.path.basename(spPlateFile)
    FXplate = fits.open(spPlateFile, memmap=True)

    #- Remove spurious EXPID** if needed
    if 'EXPID**' in FXplate[0].header:
        FXplate[0].header.remove('EXPID**')

    code_version = FXplate[0].header['RUN2D']

    spZbestFile = '%s/%s/spZbest-%s-%d.fits' % \
        (platedir, code_version, platestr, mjd)
    FXzbest = fits.open(spZbestFile, memmap=True)
    
    spZlineFile = '%s/%s/spZline-%s-%d.fits' % \
        (platedir, code_version, platestr, mjd)
    zline = fits.getdata(spZlineFile, 1)

    #- HDU0 will be a modified copy of the spPlate header
    plate_hdu = fits.PrimaryHDU(header=FXplate[0].header)

    #- if tpcorr exist, get wavelength and compute loglam for interpolation
    if tpcorr_h5 is not None :
        tpcorr_loglam = N.log10(tpcorr_h5['wave'].value)

    #- Loop over fibers on this plate on this MJD
    for fiber in fibers:
        #- HDU1 : binary table of coadd flux, log(lambda), ivar, etc.
        flux = FXplate[0].data[fiber-1]
        c0   = FXplate[0].header['COEFF0']
        c1   = FXplate[0].header['COEFF1']
        loglam = c0 + c1*N.arange(len(flux))
        ivar     = FXplate[1].data[fiber-1]
        and_mask = FXplate[2].data[fiber-1]
        or_mask  = FXplate[3].data[fiber-1]
        wdisp    = FXplate[4].data[fiber-1]
        sky      = FXplate[6].data[fiber-1]
        model    = FXzbest[2].data[fiber-1]
        calibcorr = None
        #- get spectro-photometric calibration correction if required
        if tpcorr_h5 is not None :            
            tpcorr_key = '%s/%s/%s' % (plate, mjd, fiber)
            entry=tpcorr_h5.get(tpcorr_key)
            if entry is None :
                #print "didn't find correction for %s"%tpcorr_key
                calibcorr=N.ones((flux.size))
            else :
                #print "FOUND correction for %s"%tpcorr_key
                calibcorr=N.interp(loglam,tpcorr_loglam,entry.value)
        
        #- trim off leading and trailing ivar=0 bins,
        #- but keep ivar=0 bins in the middle of the spectrum
        igood = good_ivar(ivar, fiber=fiber)
        new_coeff0 = round(float(loglam[igood[0]]),4) #- fix float32 rounding

        #- Create coadded spectrum table for HDU 1
        cols = list()
        cols.append( fits.Column(name='flux',     format='E', array=flux[igood]) )
        cols.append( fits.Column(name='loglam',   format='E', array=loglam[igood]) )
        cols.append( fits.Column(name='ivar',     format='E', array=ivar[igood]) )
        cols.append( fits.Column(name='and_mask', format='J', array=and_mask[igood]) )
        cols.append( fits.Column(name='or_mask',  format='J', array=or_mask[igood]) )
        cols.append( fits.Column(name='wdisp',    format='E', array=wdisp[igood]) )
        cols.append( fits.Column(name='sky',      format='E', array=sky[igood]) )
        cols.append( fits.Column(name='model',    format='E', array=model[igood]) )
        if tpcorr_h5 is not None :
            cols.append( fits.Column(name='calib_corr', format='E', array=calibcorr[igood]) )
        
        cols = fits.ColDefs(cols)
        coadd_hdu = fits.BinTableHDU.from_columns(cols)

        #- HDU 2: copy of spAll row
        hdux = [plate_hdu, coadd_hdu]
        ispec = N.where( (spAll.FIELD == plate) & \
                         (spAll.MJD == mjd) & \
                         (spAll.FIBERID == fiber) )[0][0]
                         
        hdux.append( fits.BinTableHDU(data=spAll[ispec:ispec+1]) )
        
        #- HDU 3: copy of rows from spZline
        ii = N.where(zline.FIBERID == fiber)[0]
        hdux.append( fits.BinTableHDU(data=zline[ii]) )
        
        #- HDU 4 .. 4+n : spectra from individual exposures
        #- Loop over individual exposures.  Do this even if we aren't
        #- writing those HDUs, so that we can update the headers with
        #- which exposures went into the coadd
        nexp = 0
        fullexpids = list()
        for iexp in range(1, 100):
            key = 'EXPID%03d' % iexp
            if key not in FXplate[0].header:
                break

            expid = FXplate[0].header[key][0:11]  #- e.g. b1-00123456

            #- check camera for this fiber
            camera = expid[0:2]
            if fiber <= 500 and camera in ('b2', 'r2'):
                continue
            elif fiber > 500 and camera in ('b1', 'r1'):
                continue
        
            if allexp and expid not in cframes:
                raise IOError, '%s not found' % expid
        
            #- If we got this far, we're going to use this exposure
            fullexpids.append(FXplate[0].header[key])
            nexp += 1
        
            if allexp:
                nfiber = cframes[expid].header['NAXIS2']
                ifib = (fiber-1) % nfiber
                d = cframes[expid]

                #- trim off leading and trailing ivar=0 bins,
                #- but keep ivar=0 bins in the middle of the spectrum
                igood = good_ivar(d.ivar[ifib], fiber='%s %d' % (expid, fiber))
        
                cols = list()
                cols.append( fits.Column(name='flux',   format='E', array=d.flux[ifib][igood]) )
                cols.append( fits.Column(name='loglam', format='E', array=d.loglam[ifib][igood]) )
                cols.append( fits.Column(name='ivar',   format='E', array=d.ivar[ifib][igood]) )
                cols.append( fits.Column(name='mask',   format='J', array=d.mask[ifib][igood]) )
                cols.append( fits.Column(name='wdisp',  format='E', array=d.wdisp[ifib][igood]) )
                cols.append( fits.Column(name='sky',    format='E', array=d.sky[ifib][igood]) )
                cols.append( fits.Column(name='calib',  format='E', array=d.calib[ifib][igood]) )
                cols.append( fits.Column(name='x',      format='E', array=d.x[ifib][igood]) )

                #- Place holder - someday we may want to calculate and include
                #- the "extra" variance which isn't proportional to the signal.
                ### n = len(d.flux[ifib])
                ### cols.append( fits.Column(name='var_extra',  format='E', array=N.zeros(n) ) )
        
                cols = fits.ColDefs(cols)
                hdux.append( fits.BinTableHDU.from_columns(cols, header=d.header) )

        #- Convert to fits HDUList
        hdux = fits.HDUList( hdux )
            
        #- Change some keyword headers which don't make sense when
        #- converting a plate header into a single object header

        #- HDU 0 is now a blank image, so fitsverify doesn't like CRPIX1 etc.
        hdr = hdux[0].header
        for key in ['CRPIX1', 'CRVAL1', 'CTYPE1', 'CD1_1']:
            if key in hdr:
                hdr.remove(key) 

        #- We trimmed leading/trailing ivar=0 pixels, so update COEFF0
        hdr['COEFF0']= new_coeff0

        #- Remove original expid list which has both SP1 and SP2
        nexp_orig = hdr['NEXP']
        del hdr['NEXP']
        for iexp in range(nexp_orig):
            expid = "EXPID%03d" % (iexp+1, )
            if expid in hdr:
                del hdr[expid]
            
        #- Add new NEXP, EXPID list for just the exposures in this file
        #- Update EXTNAME of individual exposure HDUs with this expid
        hdr['NEXP']= (nexp, 'Number of individual exposures')
        for iexp, expid in enumerate(fullexpids):
            key = "EXPID%03d" % (iexp+1, )
            hdr[key] = expid
            if allexp:
                ### print "Setting EXTNAME for %d to %s" % (4+iexp, expid)
                hdux[4+iexp].name = expid

        #- Remove mention of the other spectrograph
        #- sp1
        if legacy:
            if fiber <= 500:            #- sp1
                for key in ['NEXP_B2', 'NEXP_R2', 'EXPT_B2','EXPT_R2']:
                    if key in hdr:
                        hdr.remove(key)
            else:                       #- sp2
                for key in ['NEXP_B1', 'NEXP_R1', 'EXPT_B1','EXPT_R1']:
                    if key in hdr:
                        hdr.remove(key)
        else:
            for key in ['NEXP_B1', 'NEXP_R1', 'EXPT_B1','EXPT_R1']:
                if key in hdr:
                    hdr.remove(key)           

        #- Delete a bunch of per-exposure keywords which came along for
        #- the ride in the spPlate header
        for keyword in """
        NGUIDE 
        SEEING20 SEEING50 SEEING80
        RMSOFF20 RMSOFF50 RMSOFF80 AZ       ALT      AIRMASS
        DAQVER   CAMDAQ   SUBFRAME ERRCNT   SYNCERR  SLINES
        PIXERR   PLINES   PFERR    DIDFLUSH TAI-BEG  TAI-END
        DATE-OBS OBJSYS   ROTTYPE  ROTPOS   BOREOFF  ARCOFFX  ARCOFFY
        OBJOFFX  OBJOFFY  CALOFFX  CALOFFY  CALOFFR  GUIDOFFX GUIDEOFFY
        GUIDOFFR FOCUS
        M2PISTON M2XTILT  M2YTILT  M2XTRAN  M2YTRAN
        M1PISTON M1XTILT  M1YTILT  M1XTRAN  M1YTRAN
        SCALE    POINTING GUIDER1  SLITID1  SLIDID2  GUIDERN
        COLLA    COLLB    COLLC
        HARTMANN MC1HUMHT MC1HUMCO MC1TEMDN MC1THT   MC1TRCB  MC1TRCT
        MC1TBCB  MC1TBCT  AUTHOR   TWOPHASE XSIGMA   XSIGMIN  XSIGMAX
        WSIGMA   WSIGMIN  WSIGMAX  LAMPLIST SKYLIST  UNAME
        """.split():
            try:
                del hdr[keyword]
            except:
                pass

        #- Add some additional header keywords
        hdr['PLUG_RA'] =  (spAll.PLUG_RA[ispec],  'RA of object [deg]')
        hdr['PLUG_DEC'] = (spAll.PLUG_DEC[ispec], 'dec of object [deg]')
        hdr['THING_ID'] = (spAll.THING_ID[ispec], 'Unique object identifier')
        hdr['FIBERID'] =  (spAll.FIBERID[ispec],  'Fiber number (1-1000)')

        #- Update other headers with useful comments
        hdux[1].header.add_comment('Coadded spectrum')
        hdux[1].name = 'COADD'
        hdux[2].header.add_comment('Metadata from spAll row')
        hdux[2].name = 'SPALL'
        hdux[3].header.add_comment('Line fits from spZline')
        hdux[3].name = 'SPZLINE'

        #- BUNIT is invalid for binary table HDUs
        for i in range(1, len(hdux)):
            if 'BUNIT' in hdux[i].header:
                del hdux[i].header['BUNIT']

        #- Write final file
        outfile = '%s/spec-%s-%d-%04d.fits' % (outdir, platestr, mjd, fiber)
        ### print mjd, os.path.basename(outfile)
        try:
            hdux.writeto(outfile, overwrite=True, output_verify='fix')
        except fits.VerifyError, err:
            print "Unable to write %s" % outfile
            raise err
        
    #- Done with this plate-mjd; close input files
    FXplate.close()
    FXzbest.close()

def get_selection_doc(opts):
    """
    Return a documentation string about which data cuts were applied.
    """
    doc = list()
    doc.append("Object selection criteria:")
    if opts.plates is not None:
        if opts.subset == "ALL":
            doc.append("    All plates")
        else:
            doc.append("    Plates: " + ", ".join(map(str, opts.plates)) )
    if opts.fibers is not None:
        doc.append("    Fibers: " + opts.fibers_orig )
    else:
        if opts.subset == "ALL":
            doc.append("    All objects kept")
            doc.append("    Optional spAll-XXX subsets defined by")
            doc.append("      qso:")
            doc.append("        - Targetted as QSOs")
            doc.append("        - Targetted as GALAXY but CLASS=QSO")
            doc.append("        - FPG scan IDed as QSO")
            doc.append("        - QSO ancillary programs")
            doc.append("      gal:  OBJTYPE=GALAXY or CLASS=GALAXY")
            doc.append("      star: OBJTYPE=SPECTROPHOTO_STD or CLASS=STAR")
            doc.append("      std:  OBJTYPE=SPECTROPHOTO_STD")
            doc.append("      sky:  OBJTYPE=SKY")
        elif opts.subset == 'QSO':
            doc.append("    Only quasar targets:")
            doc.append("      - Targetted as QSOs")
            doc.append("      - Targetted as GALAXY but CLASS=QSO")
            doc.append("      - FPG scan IDed as QSO")
            doc.append("      - QSO ancillary programs")
        elif opts.subset == 'GALAXY':
            doc.append("    Only galaxies: OBJTYPE=GALAXY or CLASS=GALAXY")
        elif opts.subset == 'STAR':
            doc.append("    Only stars: OBJTYPE=SPECTROPHOTO_STD or CLASS=STAR")
        elif opts.subset in ('STD', 'SPECTROPHOTO_STD'):
            doc.append("    Only SpecPhoto standard stars: OBJTYPE=SPECTROPHOTO_STD")
        elif opts.subset == 'SKY':
            doc.append("    Only sky fibers: OBJTYPE=SKY")

    return "\n".join(doc)

def write_file_list(filename, spectra):
    FX = open(filename, 'w')
    for plate, mjd, fiber in sorted( zip(spectra.FIELD, spectra.MJD, spectra.FIBERID) ):
        platestr = plate_to_string(plate)
        specfile = "%s/spec-%s-%05d-%04d.fits" % (platestr, platestr, mjd, fiber)
        print >> FX, specfile

    FX.close()

def parse_string_range(s):
    """
    e.g. "1,2,5-8,20" -> [1,2,5,6,7,8,20]

    modified from Sven Marnach,
    http://stackoverflow.com/questions/5704931/parse-string-of-integer-sets-with-intervals-to-list

    Feature/Bug: Only works with positive numbers
    """
    ranges = (x.split("-") for x in s.split(","))
    x = [i for r in ranges for i in range(int(r[0]), int(r[-1]) + 1)]
    return x

def check_options(opts, args):
    """Sanity check options"""
    if opts.spall is None:
        parser.print_help()
        print "You must specify --spall"
        sys.exit(1)

    if opts.outdir is None:
        parser.print_help()
        print "You must specify --outdir"
        sys.exit(1)
        
    if opts.fibers is not None:
        if opts.plates is None or len(opts.plates) != 1:
            print "If you specify fibers, you must specify one and only one plate"
            sys.exit(1)
    
#-----------------------------------------------------------------------------
#- Parse command line options and call subroutines

import optparse

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-s", "--spall",  type="string",  help="input spAll file")
parser.add_option("-i", "--indir",  type="string",  help="input directory [$BOSS_SPECTRO_REDUX/$RUN2D/]")
parser.add_option("-o", "--outdir", type="string",  help="output directory")
parser.add_option("-m", "--meta",   action='store_true', help="only write top level metadata (README, spSome)")
parser.add_option("-u", "--update", action='store_true', help="update missing plates; don't overwrite others")
parser.add_option("-p", "--plates", type="string",  help="plates to process (comma separated, no spaces) default to all plates")
parser.add_option("-M", "--mjd", type="string",  help="mjds to process (comma separated, no spaces)")
parser.add_option("-c", "--coadd",  action='store_true', help="Only write coadded spectrum (no individual exposures)")
parser.add_option("-f", "--fibers", type="string", help="Comma separated list of fibers")
parser.add_option("-S", "--subset", type="string", default='ALL', help="Subset of objects to process [ALL, QSO, GALAXY, STAR, STD, or SKY]")
parser.add_option("-C", "--tpcorr", type="string", default=None, help="add a column with the spectrophotometric calibration correction for targets with LAMBDA_EFF=4000A, argument is the path to the tpcorr.hdf5 file, see http://darkmatter.ps.uci.edu/tpcorr/")
parser.add_option("-P", "--platef", action='store_true',  help="set the plate format input")
parser.add_option("-L", "--legacy", action='store_true',  help="set the sdss legacy format input")

opts, args = parser.parse_args()

#- Expand comma separated lists into arrays of integers
if opts.plates is not None:
    opts.plates = [int(x) for x in opts.plates.split(',')]
if opts.mjd is not None:
    opts.mjd = [int(x) for x in opts.mjd.split(',')]
if opts.fibers is not None:
    opts.fibers_orig = opts.fibers
    opts.fibers = parse_string_range(opts.fibers)

#- Open tpcorr file is option is set
tpcorr_h5=None
if opts.tpcorr is not None:
    print "Reading spectro-photometric calibration correction file",opts.tpcorr
    if not os.path.isfile(opts.tpcorr) :
        print "ERROR: cannot open",opts.tpcorr
        sys.exit(1)
    tpcorr_h5     = h5py.File(opts.tpcorr, 'r')
    

#- Sanity check
check_options(opts, args)

if not os.path.isdir(opts.outdir):
    os.makedirs(opts.outdir)

#- Load spAllFile
print "Reading spAll file", asctime()
try:
    # spectra = fits.getdata(opts.spall).view(N.recarray)
    spectra = fits.open(opts.spall)[1].data
except MemoryError:
    print "ERROR: Not enough memory to read the spAll file."
    print "If you are on riemann, try again from an interactive batch job:"
    print "    qsub -I -q fast -X -V"
    sys.exit(1)

#- Default input directory is BOSS_SPECTRO_REDUX/RUN2D,
#- with RUN2D from spAll (*not* environment variable)
#- Assumes spAll has one and only one RUN2D
if len(set(spectra.RUN2D)) == 1:
    run2d = spectra.RUN2D[0]
else:
    print >> sys.stderr, "ERROR: spAll file has more than one RUN2D."
    print >> sys.stderr, set(spectra.RUN2D)
    sys.exit(2)
    
if opts.indir is None:
    datadir = os.path.join(os.environ['BOSS_SPECTRO_REDUX'], run2d)
else:
    datadir = opts.indir

#- If plates aren't specified, use all of them
if opts.plates is None:
    opts.plates = sorted(set(spectra.FIELD))
    print "Using all %d plates" % len(opts.plates)

QSO_A1 = QSO_A2 = 0
QSO_A1  |= 2**22  # QSO_AAL
QSO_A1  |= 2**23  # QSO_AALS
QSO_A1  |= 2**24  # QSO_IAL
QSO_A1  |= 2**25  # QSO_RADIO
QSO_A1  |= 2**26  # QSO_RADIO_AAL
QSO_A1  |= 2**27  # QSO_RADIO_IAL
QSO_A1  |= 2**28  # QSO_NOAALS
QSO_A1  |= 2**29  # QSO_GRI
QSO_A1  |= 2**30  # QSO_HIZ
QSO_A1  |= 2**31  # QSO_RIZ
QSO_A2  |= 2**3   # QSO_VAR
QSO_A2  |= 2**4   # QSO_VAR_FPG
QSO_A2  |= 2**5   # RADIO_2LOBE_QSO
QSO_A2  |= 2**7   # QSO_SUPPZ
QSO_A2  |= 2**8   # QSO_VAR_SDSS
QSO_A2  |= 2**9   # QSO_WISE_SUPP

#- Keep only target type subset
if opts.fibers is None:
    if opts.subset == 'ALL':
        print "Keeping all objects"
    elif opts.subset == 'QSO':
        print "Trimming to just QSO targets"
        ii  = (spectra.OBJTYPE == 'QSO') 
        ii |= ((spectra.OBJTYPE == 'GALAXY') & (spectra.CLASS == 'QSO'))
        ii |= (spectra.CLASS_PERSON == 3)   # 3 == FPG IDed as QSO # key doesn't exist?
        #- Ancillary QSO programs
        ii |= (spectra.ANCILLARY_TARGET1 & QSO_A1)
        ii |= (spectra.ANCILLARY_TARGET2 & QSO_A2)
        spectra = spectra[ii]
    elif opts.subset == 'GALAXY' or opts.subset == 'GAL':
        opts.subset = 'GAL'
        print "Trimming to just GALAXY targets"
        ii  = (spectra.OBJTYPE == 'GALAXY') 
        ii |= (spectra.CLASS == 'GALAXY')
        spectra = spectra[ii]
    elif opts.subset == 'STAR':
        print "Trimming to just STAR targets"
        ii  = (spectra.OBJTYPE == 'SPECTROPHOTO_STD') 
        ii |= (spectra.CLASS == 'STAR')
        spectra = spectra[ii]
    elif opts.subset in ('STD', 'SPECTROPHOTO_STD'):
        print "Trimming to just spectro-photometric standard stars"
        ii  = (spectra.OBJTYPE == 'SPECTROPHOTO_STD') 
        spectra = spectra[ii]
    elif opts.subset == 'SKY':
        print "Trimming to just SKY fibers (!)"
        ii  = (spectra.OBJTYPE == 'SKY') 
        spectra = spectra[ii]
    elif opts.subset is not None:
        print "FATAL: subclass must be ALL, QSO, GALAXY, STAR, STD, or SKY"
        sys.exit(1)    
else:
    print "Fibers specified; not trimming by target type"
    ii = N.zeros(len(spectra), dtype='bool')
    for fiber in opts.fibers:
        ii |= (spectra.FIBERID == fiber)
    spectra = spectra[ii]

#- Write README and spSome files
if opts.meta:
    print "Writing README.txt"
    header = "Input data from:\n    %s\n" % datadir
    header += get_selection_doc(opts)
    write_readme(opts.outdir + '/README.txt', header=header)
    if opts.subset != "ALL":
        spSomeName = opts.outdir+'/spAll-%s-%s.fits' % \
                                        (opts.subset.lower(), run2d)
        print "Writing", os.path.basename(spSomeName)
        fits.writeto(spSomeName, spectra, overwrite=True)
        filelist = opts.outdir+'/specfiles-%s-%s.txt' % \
                                        (opts.subset.lower(), run2d)
        print "Writing", os.path.basename(filelist)
        write_file_list(filelist, spectra)
    else:
        import shutil
        shutil.copy(opts.spall, opts.outdir+'/spAll-%s.fits' % run2d)
    sys.exit(0)

#- For efficiency, process one plate at a time
print "Starting plate processing", asctime()
for plate in sorted(set(opts.plates)):    
    print 'Plate %d : %s' % (plate, asctime())
    platestr = plate_to_string(plate)
    outdir = '%s/%s/' % (opts.outdir, platestr)
    
    #- find MJDs for this plate
    ii = N.where(spectra.FIELD == plate)[0]
    plate_mjds = sorted(set(spectra.MJD[ii]))

    #- Filter by mjd option if given
    if opts.mjd is not None:
        plate_mjds = sorted(set(plate_mjds) & set(opts.mjd))

    #- Create output directory if needed
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    for mjd in plate_mjds:
        #- Process fibers for just this PLATE-MJD
        ii = N.where((spectra.FIELD == plate) & (spectra.MJD == mjd))
        fibers = spectra.FIBERID[ii]
        #- If --update option is True, select only fibers where there is no spec files
        if opts.update:
            fibers = [f for f in fibers if not os.path.exists(os.path.join(outdir, 'spec-%s-%d-%04d.fits'%(platestr, mjd, f)))]
            fibers =  N.array(fibers)
            print "Updating only missing files. %d fibers found for this plate" % fibers.size

        process_plate(datadir, outdir, plate, mjd, fibers, spectra, allexp=not opts.coadd,tpcorr_h5=tpcorr_h5, plate_f=opts.platef, legacy=opts.legacy)
            
print "Wrote files to " + opts.outdir
print "Done", asctime()
            
        
