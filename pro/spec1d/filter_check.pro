;+
; NAME:
;   filter_check
;
; PURPOSE:
;   
;   Gather spectra based on an input file of the form 
;   created by platemerge (the spAll file). Calculate the 
;   ugriz throughput for each object in the plates, possibly
;   putting limits on target type, MJD, or signal-to-noise
;   (essentially by requiring survey quality). 
;
; CALLING SEQUENCE:
;   res = filter_check( spallfile, outfile, [mjdlimits= , primtarget=,
;   filter_prefix=, mingisn2=])
;
; INPUTS:
;   spallfile  - spAll.fit file as created by platemerge
;   filter_prefix  - Use alternate prefix for filter curves to use
;                    (allowed are sdss or doi) 
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;   mjdlimits  - Only look in a certain range of MJDs
;   primtarget - Require a certain target type
;   mingisn2 - Minimum plate SN^2 in g AND i
;
; OUTPUTS:
;   outfile    - Fits file with all the spAll.fit info, but with
;                synthetic ugriz replaced with the desired filter 
;                curves
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   filter_thru()
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/sdss_u_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_g_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_r_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_i_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_z_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_u_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_g_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_r_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_i_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_z_atm.dat
;
; REVISION HISTORY:
;   05-APr-2000  Written by M. Blanton, Fermilap
;-
;------------------------------------------------------------------------------
pro filter_check, spallfile, outfile, filter_prefix, mjdlimits=mjdlimits, $
                  primtarget=primtarget, mingisn2=mingisn2

;--------
; Read in the platemerge file 
spall=mrdfits(spallfile,1)
nall=(size(spall))[1]

;--------
; Select desired spectra
select=intarr(nall)
select[*]=1
if(keyword_set(mjdlimits)) then $
  select[where(spall.mjd lt mjdlimits[0] or spall.mjd gt mjdlimits[1]]=0
if(keyword_set(primtarget)) then $
  select[where((spall.primtarget and primtarget) eq 0)]=0
if(keyword_set(mingisn2)) then $
  select[where(spall.spec1_g lt mingisn2 or $
               spall.spec2_g lt mingisn2 or $
               spall.spec1_i lt mingisn2 or $
               spall.spec2_i lt mingisn2 )]=0
spselect=spall[where(select gt 0)]
nselect=(size(spselect))[1]

;--------
; Gather all of the spectra 
readspec,spselect.plate,spselect.fiberid,mjd=spselect.mjd, $
  flux=flux,flerr=flerr,invvar=invvar,loglam=loglam,synflux=synflux

;--------
; Send the spectra and fluxes to the synthesizer;
; an attempt to do here what is done in spreduce1d.pro
waveimg=10d^loglam 
npixobj=n_elements(loglam)
flambda2fnu = waveimg^2 / 2.99792e18
fthru=filter_thru(flux*rebin(flambda2fnu, npixobj, nselect), waveimg=waveimg, $
                  mask=(invvar eq 0), /norm, filter_prefix=filter_prefix, $
                  /converttoair)
synfthru=filter_thru(synflux*rebin(flambda2fnu, npixobj, nselect), $
                     waveimg=waveimg, /norm, filter_prefix=filter_prefix, $
                     /converttoair)

spselect.counts_spectro=fthru
spselect.counts_synth=synfthru

mwrfits,outfile,spselect

end
;------------------------------------------------------------------------------
