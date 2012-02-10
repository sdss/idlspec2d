;+
; NAME:
;   filter_select
;
; PURPOSE:
;   Gather spectra based on an input file of the form created by platemerge (the spAll file).
;   Calculate the ugriz throughput for each object in the plates, possibly
;   putting limits on target type, MJD, or signal-to-noise
;   (essentially by requiring survey quality). 
;
; CALLING SEQUENCE:
;   filter_select, spallfile, outfile, [mjdlimits= , primtarget=,
;   filter_prefix=, mingisn2=, rpsfmodel=]
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
; BUGS:
;   Depends on spec_append, internal routine of readspec
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   filter_thru()
;   spec_append
;   readspec
;   mrdfits()
;   mwrfits
;
; REVISION HISTORY:
;   05-APr-2000  Written by M. Blanton, Fermiland
;-
;------------------------------------------------------------------------------
pro filter_select, spallfile, outbase, filter_prefix, mjdlimits=mjdlimits, $
                  primtarget=primtarget, mingisn2=mingisn2, firstn=firstn, $
                  rfluxlimits=rfluxlimits, rpsfmodel=rpsfmodel, $
                  flux=flux,invvar=invvar,loglam=loglam, spselect=spselect, $
                  outspectra=outspectra

  ;--------
  ; Read in the file
  spall=mrdfits(spallfile,1)
  nall=n_elements(spall)
  splog,string(nall)+' objects total'

  ;--------
  ; Select desired spectra
     
  ; set all initially to selected 
  select=intarr(nall)
  select[*]=1

  ; throw out those outside MJD limits (inclusive)
  if(keyword_set(mjdlimits)) then begin
    indx=where(spall.mjd lt mjdlimits[0] or spall.mjd gt mjdlimits[1],count)
    if(count gt 0) then select[indx]=0
  endif

  ; throw out those outside flux limit (inclusive)
  if(keyword_set(rfluxlimits)) then begin
    indx=where(spall.psfCounts[2] lt rfluxlimits[0] or $
               spall.psfCounts[2] gt rfluxlimits[1],count)
    if(count gt 0) then select[indx]=0
  endif

  ; throw out those outside flux limit (inclusive)
  if(keyword_set(rpsfmodel)) then begin
    indx=where(spall.psfcounts[2]-spall.counts_model[2] gt rpsfmodel,count)
    if(count gt 0) then select[indx]=0
  endif

  if(keyword_set(primtarget)) then begin
    indx=where((spall.primtarget and primtarget) eq 0,count) 
    if(count gt 0) then select[indx]=0
  endif

  if(keyword_set(mingisn2)) then begin
    indx=where(spall.spec1_g lt mingisn2 or $
	       spall.spec2_g lt mingisn2 or $
	       spall.spec1_i lt mingisn2 or $
	       spall.spec2_i lt mingisn2,count )
    if(count gt 0) then select[indx]=0
  endif
  indx=where(select gt 0,count)

  if(keyword_set(firstn)) then begin 
    if (count gt firstn) then count=firstn
  endif

  if(count eq 0) then begin
    splog,'no objects '
    return
  endif 

  spselect=spall[indx[0:count-1]]
  nselect=n_elements(spselect)

  ;--------
  ; Gather all of the spectra
  splog,string(nselect)+' objects selected	

  readspec,spselect.plate, $
       spselect.fiberid, $
       mjd=spselect.mjd, $
       flux=compflux,invvar=compinvvar, $
       loglam=comploglam,andmask=compandmask, ormask=compormask
  compinvvar = skymask(compinvvar, compandmask, compormask)
  compandmask = 0 ; Free memory
  compormask = 0 ; Free memory

  npixobj=n_elements(comploglam)/nselect
;  loglammin=max(comploglam[0,*])
;  loglammax=min(comploglam[npixobj-1,*])
  loglammin=3.5855
  loglammax=3.9610
  nloglam=round(10000l*(loglammax-loglammin))
  flux=fltarr(nloglam,nselect)
  invvar=fltarr(nloglam,nselect)
  loglam=fltarr(nloglam,nselect)
  for i = 0, nselect-1 do begin
    ilmin=10000*(loglammin-comploglam[0,i])
    flux[*,i]=compflux[ilmin:ilmin+nloglam-1,i]
    invvar[*,i]=compinvvar[ilmin:ilmin+nloglam-1,i]
    loglam[*,i]=comploglam[ilmin:ilmin+nloglam-1,i]
  end
  
  ;--------
  ; Send the spectra and fluxes to the synthesizer;
  ; an attempt to do here what is done in spreduce1d.pro
  waveimg=10d^loglam 
  flambda2fnu = 1.e+6 * waveimg^2 / (3631.*2.99792e18)
  flux=flux*rebin(flambda2fnu, nloglam, nselect)
  fthru=filter_thru(flux, waveimg=waveimg, mask=(invvar eq 0), $
   filter_prefix=filter_prefix, /toair)
       
  for band=0, 4 do begin 
    spselect[*].counts_spectro[band]=fthru[*,band]
  endfor

  if (keyword_set(outspectra)) then begin
    mwrfits,spselect,outbase+'.fits',/create
    openw,11,outbase+'.flux'
    openw,12,outbase+'.invvar'
    openw,13,outbase+'.loglam'
    writeu,11,float(flux)
    writeu,12,float(invvar)
    writeu,13,float(loglam)
    close,11
    close,12
    close,13
  endif

end
;------------------------------------------------------------------------------
