; Select spectroscopically-confirmed LRGs, and compute photo-z's
; for comparison to spectroscopic redshifts.

pro lrg_photoz_test

   ; Read the spAll file
   spfile = filepath('spAll.fits', root_dir=getenv('SPECTRO_DATA'))
   columns = ['PROGNAME', 'PLATEQUALITY', 'PLATE', 'FIBERID', 'MJD', $
    'CLASS', 'SPECPRIMARY', 'PRIMTARGET', 'Z', 'Z_ERR', 'ZWARNING', $
    'MODELFLUX', 'MODELFLUX_IVAR', 'EXTINCTION']
   spall = hogg_mrdfits(spfile, 1, columns=columns, $
    nrowchunk=10000L) ;, range=[100000,200000])

   ; Trim to LRGs
   itrim = where( strmatch(spall.class,'GALAXY*') $
;    AND strmatch(spall.progname,'main*') $
    AND strmatch(spall.platequality,'good*') $
    AND spall.specprimary EQ 1 $ ; Best spectroscopic observations
    AND total(spall.modelflux_ivar GT 0,1) EQ 5 $ ; Good photom in all bands
    AND (spall.primtarget AND (2L^5+2L^26)) NE 0 $ ; Targetted as LRG
    AND spall.zwarning EQ 0, ntrim )
   spall = spall[itrim]

   ; Do the photo-z fits
   zfit = lrg_photoz(spall.modelflux, spall.modelflux_ivar, z_err=zfit_err, $
    extinction=spall.extinction, /abcorrect, chi2=chi2)
   zdiff = zfit - spall.z

   splot, spall.z, zfit, psym=3, $
    xtitle='Spectroscopic Z', ytitle='Photometric Z'
   soplot, [0,1], [0,1], color='red'
stop

   return
end
