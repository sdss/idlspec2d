#!/usr/bin/env python

"""
Reformat spectra into a single fits file per object, combining all necessary
pieces from spPlate, spCFrame, spFrame, spFlat, and spZbest

Stephen Bailey, Summer 2011

Bugs/Features:
  - Assumes RUN1D=RUN2D, and gets them from spAll (not env vars), assuming
    that spAll has one and only one RUN2D.
"""

import sys
import os
import os.path
from glob import glob           #- File pattern globbing (spFrame*.fits)
import re                       #- Regular expressions
from time import asctime
import numpy as N
import pyfits
import fitsio

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

These are grouped in subdirectories by plate, e.g.:

    README.txt    : this file
    spXXX.fits    : subset of spAll.fits (e.g. Qso, Star, Gal)
    4080/         : dir for objects on plate 4080
        spec-4080-55368-0487.fits  : spec-PLATE-MJD-FIBER
        spec-4080-55471-0485.fits  : different plugging, same plate
    ...

Note that the same THING_ID may have been observed on more than
one plate.  Use the THING_ID, PLATE, MJD, FIBERID metadata in
spSome to find the files you need.

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

HDU 2 : Copy of row for this object from spAll table
    http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/spAll.html

HDU 3 : Copy of rows for this object from spZline table
    http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/RUN1D/spZline.html"""

    #- Continue with explanation of per-exposure HDUs if included
    if allexp:
        print >> fx, """
HDU 4 .. n+4 : Individual frames.
    For each exposure there is one HDU for the red camera and one for the blue.
    These are in the order of the EXPIDnn keywords in the HDU0 header.

    Header: Taken from HDU0 of individual spCFrame files
    Data: Binary table with columns taken from spCFrame files:
        flux
        loglam
        ivar
        mask
        wdisp
        sky
        calib : conversion between flux and electrons: flux = electrons*calib

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
        eflux = pyfits.getdata(framefile, 0)
        nfiber, npix = eflux.shape
                
        #- Load spCFrame file; trim arrays back to original size
        fx = pyfits.open(cframefile)
        self.flux = fx[0].data[:, 0:npix]
        self.ivar = fx[1].data[:, 0:npix]
        self.mask = fx[2].data[:, 0:npix]
        self.loglam = fx[3].data[:, 0:npix]
        self.wdisp  = fx[4].data[:, 0:npix]
        self.sky    = fx[6].data[:, 0:npix]
        self.header = fx[0].header
        
        #- Load superflat spCFrame[8] and fiberflat spFlat[0]
        filedir, basename = os.path.split(cframefile)
        superflat = fx[8].data[:, 0:npix]
        flatfile = fx[0].header['FLATFILE'].replace('sdR', 'spFlat')
        flatfile = flatfile.replace('.fit', '.fits.gz')
        flatfile = os.path.join(filedir, flatfile)
        fiberflat = pyfits.getdata(flatfile, 0)
        
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

def good_ivar(ivar):
    """
    return indices of for array from first non-zero ivar to the last
    non-zero ivar.  i.e. trim off leading and trailing contiguously zero ivars.
    
    If all ivar==0, return indices for full array since code might get
    grumpy with completely blank arrays.
    """
    if N.all(ivar == 0):
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
        
def process_plate(datadir, outdir, plate, mjd, fibers, spAll, allexp=True):
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
        
    Outputs:
        writes files to outdir/plate/spec-plate-mjd-fiber.fits
    """
    #- Load all C/Frame files for this plate
    platedir = '%s/%d/' % (datadir, plate)
    if allexp:
        cframes = load_spCFrame_files(platedir)

    #- Open spPlate, spZbest, and spZline files
    spPlateFile = '%s/spPlate-%d-%d.fits' % (platedir, plate, mjd)
    print 'Processing', os.path.basename(spPlateFile)
    FXplate = pyfits.open(spPlateFile, memmap=True)
    code_version = FXplate[0].header['RUN2D']

    spZbestFile = '%s/%s/spZbest-%d-%d.fits' % \
        (platedir, code_version, plate, mjd)
    FXzbest = pyfits.open(spZbestFile, memmap=True)
    
    spZlineFile = '%s/%s/spZline-%d-%d.fits' % \
        (platedir, code_version, plate, mjd)
    zline = pyfits.getdata(spZlineFile, 1)

    #- HDU0 will be a modified copy of the spPlate header
    plate_hdu = pyfits.PrimaryHDU(header=FXplate[0].header)
    
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

        #- trim off leading and trailing ivar=0 bins,
        #- but keep ivar=0 bins in the middle of the spectrum
        igood = good_ivar(ivar)

        #- Create coadded spectrum table for HDU 1
        cols = list()
        cols.append( pyfits.Column(name='flux',     format='E', array=flux[igood]) )
        cols.append( pyfits.Column(name='loglam',   format='E', array=loglam[igood]) )
        cols.append( pyfits.Column(name='ivar',     format='E', array=ivar[igood]) )
        cols.append( pyfits.Column(name='and_mask', format='J', array=and_mask[igood]) )
        cols.append( pyfits.Column(name='or_mask',  format='J', array=or_mask[igood]) )
        cols.append( pyfits.Column(name='wdisp',    format='E', array=wdisp[igood]) )
        cols.append( pyfits.Column(name='sky',      format='E', array=sky[igood]) )
        cols.append( pyfits.Column(name='model',    format='E', array=model[igood]) )
        
        cols = pyfits.ColDefs(cols)
        coadd_hdu = pyfits.new_table(cols)

        #- HDU 2: copy of spAll row
        hdux = [plate_hdu, coadd_hdu]
        ispec = N.where( (spAll.PLATE == plate) & \
                         (spAll.MJD == mjd) & \
                         (spAll.FIBERID == fiber) )[0][0]
                         
        hdux.append( pyfits.BinTableHDU(data=spAll[ispec:ispec+1]) )
        
        #- HDU 3: copy of rows from spZline
        ii = N.where(zline.FIBERID == fiber)[0]
        hdux.append( pyfits.BinTableHDU(data=zline[ii]) )
        
        #- HDU 4 .. 4+n : spectra from individual exposures
        #- Loop over individual exposures.  Do this even if we aren't
        #- writing those HDUs, so that we can update the headers with
        #- which exposures went into the coadd
        nexp = 0
        fullexpids = list()
        for iexp in range(1, 100):
            key = 'EXPID%02d' % iexp
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
                igood = good_ivar(d.ivar[ifib])
        
                cols = list()
                cols.append( pyfits.Column(name='flux',   format='E', array=d.flux[ifib][igood]) )
                cols.append( pyfits.Column(name='loglam', format='E', array=d.loglam[ifib][igood]) )
                cols.append( pyfits.Column(name='ivar',   format='E', array=d.ivar[ifib][igood]) )
                cols.append( pyfits.Column(name='mask',   format='J', array=d.mask[ifib][igood]) )
                cols.append( pyfits.Column(name='wdisp',  format='E', array=d.wdisp[ifib][igood]) )
                cols.append( pyfits.Column(name='sky',    format='E', array=d.sky[ifib][igood]) )
                cols.append( pyfits.Column(name='calib',  format='E', array=d.calib[ifib][igood]) )

                #- Place holder - someday we may want to calculate and include
                #- the "extra" variance which isn't proportional to the signal.
                ### n = len(d.flux[ifib])
                ### cols.append( pyfits.Column(name='var_extra',  format='E', array=N.zeros(n) ) )
        
                cols = pyfits.ColDefs(cols)
                hdux.append( pyfits.new_table(cols, header=d.header) )

        #- Convert to pyfits HDUList
        hdux = pyfits.HDUList( hdux )
            
        #- Change some keyword headers which don't make sense when
        #- converting a plate header into a single object header

        #- HDU 0 is now a blank image, so fitsverify doesn't like CRPIX1 etc.
        hdr = hdux[0].header
        del hdr['CRPIX1']
        del hdr['CRVAL1']
        del hdr['CTYPE1']
        del hdr['CD1_1']
        
        #- Remove original expid list which has both SP1 and SP2
        nexp_orig = hdr['NEXP']
        del hdr['NEXP']
        for iexp in range(nexp_orig):
            expid = "EXPID%02d" % (iexp+1, )
            del hdr[expid]
            
        #- Add new NEXP, EXPID list for just the exposures in this file
        #- Update EXTNAME of individual exposure HDUs with this expid
        hdr.update('NEXP', nexp, 'Number of individual exposures')
        for iexp, expid in enumerate(fullexpids):
            key = "EXPID%02d" % (iexp+1, )
            hdr.update(key, expid)
            if allexp:
                ### print "Setting EXTNAME for %d to %s" % (4+iexp, expid)
                hdux[4+iexp].update_ext_name(expid)

        #- Remove mention of the other spectrograph
        #- sp1
        if fiber <= 500:            #- sp1
            del hdr['NEXP_B2']
            del hdr['NEXP_R2']
            del hdr['EXPT_B2']
            del hdr['EXPT_R2']                
        else:                       #- sp2
            del hdr['NEXP_B1']
            del hdr['NEXP_R1']
            del hdr['EXPT_B1']
            del hdr['EXPT_R1']

        #- Delete a bunch of per-exposure keywords which came along for
        #- the ride in the spPlate header
        for keyword in """
        MJD      NGUIDE 
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
            del hdr[keyword]

        #- Add some additional header keywords
        hdr.update('PLUG_RA',  spAll.PLUG_RA[ispec],  'RA of object [deg]')
        hdr.update('PLUG_DEC', spAll.PLUG_DEC[ispec], 'dec of object [deg]')
        hdr.update('THING_ID', spAll.THING_ID[ispec], 'Unique object identifier')
        hdr.update('FIBERID',  spAll.FIBERID[ispec],  'Fiber number (1-1000)')

        #- Update other headers with useful comments
        hdux[1].header.add_comment('Coadded spectrum')
        hdux[1].update_ext_name('COADD')
        hdux[2].header.add_comment('Metadata from spAll row')
        hdux[2].update_ext_name('SPALL')
        hdux[3].header.add_comment('Line fits from spZline')
        hdux[3].update_ext_name('SPZLINE')

        #- BUNIT is invalid for binary table HDUs
        for i in range(1, len(hdux)):
            if 'BUNIT' in hdux[i].header:
                del hdux[i].header['BUNIT']

        #- Write final file
        outfile = '%s/spec-%d-%d-%04d.fits' % (outdir, plate, mjd, fiber)
        ### print mjd, os.path.basename(outfile)
        try:
            hdux.writeto(outfile, clobber=True)
        except pyfits.core.VerifyError:
            print "Unable to write %s" % outfile
        
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
            doc.append("    Optional spXXX subsets defined by")
            doc.append("      Qso:")
            doc.append("        - Targetted as QSOs")
            doc.append("        - Targetted as GALAXY but CLASS=QSO")
            doc.append("        - FPG scan IDed as QSO")
            doc.append("        - QSO ancillary programs")
            doc.append("      Gal:  OBJTYPE=GALAXY or CLASS=GALAXY")
            doc.append("      Star: OBJTYPE=SPECTROPHOTO_STD or CLASS=STAR")
            doc.append("      Std:  OBJTYPE=SPECTROPHOTO_STD")
            doc.append("      Sky:  OBJTYPE=SKY")
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
parser.add_option("-c", "--coadd",  action='store_true', help="Only write coadded spectrum (no individual exposures)")
parser.add_option("-f", "--fibers", type="string", help="Comma separated list of fibers")
parser.add_option("-S", "--subset", type="string", default='ALL', help="Subset of objects to process [ALL, QSO, GALAXY, STAR, STD, or SKY]")

opts, args = parser.parse_args()

#- Expand comma separated lists into arrays of integers
if opts.plates is not None:
    opts.plates = [int(x) for x in opts.plates.split(',')]
if opts.fibers is not None:
    opts.fibers_orig = opts.fibers
    opts.fibers = parse_string_range(opts.fibers)

#- Sanity check
check_options(opts, args)

if not os.path.isdir(opts.outdir):
    os.makedirs(opts.outdir)

#- Load spAllFile
print "Reading spAll file", asctime()
try:
    # spectra = pyfits.getdata(opts.spall).view(N.recarray)
    spectra = fitsio.read(opts.spall, 1).view(N.recarray)
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
    opts.plates = sorted(set(spectra.PLATE))
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
        ii |= (spectra.CLASS_PERSON == 3)   # 3 == FPG IDed as QSO
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
    print "Writing spSome and README files"
    header = "Input data from:\n    %s\n" % datadir
    header += get_selection_doc(opts)
    write_readme(opts.outdir + '/README.txt', header=header)
    if opts.subset != "ALL":
        spSomeName = opts.outdir+'/sp%s-%s.fits' % (opts.subset.title(), run2d)
        pyfits.writeto(spSomeName, spectra, clobber=True)
    else:
        import shutil
        shutil.copy(opts.spall, opts.outdir+'/spAll-%s.fits' % run2d)
    sys.exit(0)

#- For efficiency, process one plate at a time
print "Starting plate processing", asctime()
for plate in sorted(set(opts.plates)):    
    print 'Plate %d : %s' % (plate, asctime())
    outdir = '%s/%04d/' % (opts.outdir, plate)
    
    #- find MJDs for this plate
    ii = N.where(spectra.PLATE == plate)[0]
    plate_mjds = sorted(set(spectra.MJD[ii]))

    #- Create output directory if needed
    #- If it already exists and --update option is True, move to next plate
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    elif opts.update:
        continue
    
    for mjd in plate_mjds:
        #- Process fibers for just this PLATE-MJD
        ii = N.where((spectra.PLATE == plate) & (spectra.MJD == mjd))
        fibers = spectra.FIBERID[ii]
        process_plate(datadir, outdir, plate, mjd, fibers, spectra, allexp=not opts.coadd)
            
print "Wrote files to " + opts.outdir
print "Done", asctime()
            
        
