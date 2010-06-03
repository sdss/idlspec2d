; Select spectroscopically-confirmed LRGs, and compute photo-z's
; for comparison to spectroscopic redshifts.

pro lrg_photoz_test

   ; Read the spAll file
   spfile = filepath('spAll.fits', root_dir=getenv('BOSS_SPECTRO_REDUX'))
   columns = ['PROGNAME', 'PLATEQUALITY', 'PLATE', 'FIBERID', 'MJD', $
    'RUN', 'RERUN', 'CAMCOL', 'FIELD', 'ID', $
    'CLASS', 'SPECPRIMARY', 'PRIMTARGET', 'Z', 'Z_ERR', 'ZWARNING', $
    'MODELFLUX', 'MODELFLUX_IVAR', 'EXTINCTION','NMGYPERCOUNT','COLC']
   spall = hogg_mrdfits(spfile, 1, columns=columns, $
    nrowchunk=10000L) ;, range=[100000,200000])

   ; Trim to LRGs
   itrim = where(strmatch(spall.platequality,'good*') $
;    AND strmatch(spall.progname,'main*') $
;    AND strmatch(spall.class,'GALAXY*') $
;    AND total(spall.modelflux_ivar GT 0,1) EQ 5 $ ; Good photom in *all* bands
    AND spall.modelflux_ivar[1] $ ; Good photom in g-band
    AND spall.modelflux_ivar[2] $ ; Good photom in r-band
    AND spall.modelflux_ivar[3] $ ; Good photom in i-band
    AND (spall.primtarget AND (2L^5+2L^26)) NE 0 $ ; Targetted as LRG
;    AND spall.zwarning EQ 0 $
    AND spall.specprimary EQ 1, ntrim) ; Best spectroscopic observations
   spall = spall[itrim]

   ; Do the photo-z fits
   zfit = lrg_photoz(spall.modelflux, spall.modelflux_ivar, z_err=zfit_err, $
    extinction=spall.extinction, /abcorrect, chi2=chi2)
   zdiff = zfit - spall.z

   azdiff = abs(zfit - spall.z)
   print, " zMin zMax   Nobj  <z-zFit>   68%cf   90%cf   95%cf   99%cf"
   for i=0, 13 do begin
      if (i LT 13) then begin
         zmin = 0.05 * i
         zmax = 0.05 * (i+1)
      endif else begin
         zmin = 0.10
         zmax = 0.60
      endelse
      indx = where(spall.z GT zmin AND spall.z LE zmax, ct)
      isort = indx[sort(azdiff[indx])]
      i68 = isort[0.68*ct]
      i95 = isort[0.95*ct]
      i99 = isort[0.99*ct]
      print, zmin, zmax, ct, median(zdiff[indx]), $
       azdiff[i68], azdiff[i95], azdiff[i99], $
       format='(2f5.2,i7,f10.4,3f8.3)'
   endfor

   splot, spall.z, zfit, psym=3, yr=[0,0.6], $
    xtitle='Spectroscopic Z', ytitle='Photometric Z'
   soplot, [0,1], [0,1], color='red'
save

   return
end
