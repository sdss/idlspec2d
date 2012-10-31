; Look at spectro-photo residuals from QSOs
;   TOPDIR - Top-level directory for plate list and spectra
;   RUN2D  - RUN2D version for plate list and spectra
;   SYNTOPDIR - Top-level directory for SYNFLUX
;   SYNRUN2D  - RUN2D version for SYNFLUX
;   SYNRUN1D  - RUN1D version for SYNFLUX
;   MINFLUX   - Use data only where SYNFLUX is greater than this value;
;               default to 0.5
;   RLIMIT    - Trim to only objects where the spectroscopic r-band
;               magnitude is brighter than this; default to 20.5
pro sphoto_resid_qso, topdir=topdir, run2d=run2d, $
 syntopdir=syntopdir, synrun2d=synrun2d, synrun1d=synrun1d, $
 minflux=minflux1, rlimit=rlimit1

   if (keyword_set(minflux1)) then minflux = minflux1 $
    else minflux = 0.5
   if (keyword_set(rlimit1)) then rlimit = rlimit1 $
    else rlimit = 20.5

; Look at Bailey's new test reductions
topdir=getenv('BOSS_SPECTRO_REDUX')+'/test/boss'
run2d='v5_5_11a'
syntopdir=getenv('BOSS_SPECTRO_REDUX')
synrun2d='v5_5_11'
synrun1d='v5_5_11'

; Look at the DR9 reductions
;topdir=getenv('BOSS_SPECTRO_REDUX')
;run2d='v5_4_45'
;syntopdir=getenv('BOSS_SPECTRO_REDUX')
;synrun2d='v5_4_45'
;synrun1d='v5_4_45'

; Look at the DR8 reductions
;topdir='/clusterfs/riemann/raid006/dr8/common/sdss-spectro/redux'
;run2d='26'
;syntopdir=topdir
;synrun2d='26'
;synrun1d=''

   ; The plate must be good in both sets of reductions if there are two
   if (run2d EQ '26') then begin
      plist = mrdfits(topdir+'/plates-dr8.fits',1)
      plist2 = plist
   endif else begin
      platelist, plist=plist, topdir=topdir, run2d=run2d
      platelist, plist=plist2, topdir=syntopdir, run2d=synrun2d
   endelse
   maxplate = (max(plist.plate) > max(plist2.plate)) + 1
   match, plist.plate+plist.mjd*maxplate, $
    plist2.plate+plist2.mjd*maxplate, i1, i2
   igood = i1[ where(strmatch(plist[i1].platequality,'good*') $
    AND plist[i1].platesn2 GT 0 $
    AND strmatch(plist2[i2].platequality,'good*') $
    AND plist2[i2].platesn2 GT 0, ngood) ]
   splog, 'Trim to ', ngood, ' good plates'
   plist = plist[igood]

   readspec, topdir=syntopdir, run2d=synrun2d, run1d=synrun1d, $
    plist.plate, mjd=plist.mjd, $
    zans=zans, plug=plug
   rmag = 22.5 - 2.5*alog10(zans.spectrosynflux[2]>0.1)
   qq = strmatch(zans.class,'QSO*') AND zans.zwarning EQ 0 $
    AND strmatch(plug.objtype,'QSO*') AND rmag LT rlimit

   ii=where(qq,ct)
   splog, 'Reading ', ct, ' objects'
   readspec, zans[ii].plate, mjd=zans[ii].mjd, zans[ii].fiberid, $
    topdir=topdir, run2d=run2d, /align, $
    wave=wave, flux=flux, invvar=invvar
   readspec, zans[ii].plate, mjd=zans[ii].mjd, zans[ii].fiberid, $
    topdir=syntopdir, run2d=synrun2d, run1d=synrun1d, /align, $
    wave=swave, synflux=synflux

   ; Trim to the same wavelengths if necessary
   if (min(wave) NE min(swave) OR max(wave) NE max(swave)) then begin
      minwave = min(wave) > min(swave)
      maxwave = max(wave) < max(swave)
      indx1 = where(wave GE minwave AND wave LE maxwave)
      indx2 = where(swave GE minwave AND swave LE maxwave)
      wave = wave[indx1]
      flux = flux[indx1,*]
      invvar = invvar[indx1,*]
      swave = swave[indx2]
      synflux = synflux[indx2,*]
   endif

   ; Mask all lines within 100 Ang of our line list, and blueward of LyA
   thisz = zans[ii].z
   lines = yanny_readone(djs_filepath('emlines.par', $
    root_dir=getenv('IDLSPEC2D_DIR'),subdir='etc'))
   wstart = lines.lambda-100
   wend = lines.lambda+100
   wstart[0]=0. ; mask everything blueward of Lyman-alpha
   mask=(invvar GT 0) ; more aggressive masking???
   for i=0L, ct-1 $
    do for j=0, n_elements(lines)-1 do $
     mask[*,i] *= (wave/(thisz[i]+1) LT wstart[j] $
      OR wave/(thisz[i]+1) GT wend[j])

   ratio = 0*wave
   for i=0L, n_elements(wave)-1L do begin
      djs_iterstat, flux[i,*]/(synflux[i,*]>minflux), $
       invvar=invvar[i,*]*(synflux[i,*]>minflux)^2*(synflux[i,*] GT minflux), $
       mean=mn1
      ratio[i] = mn1
   endfor

   ; Output the correction vector
   datfile = 'sphoto_resid_qso-'+run2d+'.dat'
   ratio_interp = djs_maskinterp(ratio, $
    finite(ratio) EQ 0 OR ratio LT 0.5 OR ratio GT 1.5, /const)
   splog, filename=datfile
   splog, /noname, '# Multiplicative spectro-photometry errors'
   splog, /noname, '# Number of plates = ', n_elements(plist)
   splog, /noname, '# Number of objects = ', n_elements(ii)
   splog, /noname, '# TOPDIR =', topdir
   splog, /noname, '# RUN2D =', run2d
   splog, /noname, '# SYNTOPDIR =', syntopdir
   splog, /noname, '# SYNRUN2D =', synrun2d
   splog, /noname, '# SYNRUN1D =', synrun1d
   splog, /noname, '# MINFLUX =', minflux
   splog, /noname, '# RLIMIT =', rlimit
   splog, /noname, '# '
   splog, /noname, '# Log(wave)  Resid'
   for j=0L, n_elements(wave)-1L do $
    splog, /noname, alog10(wave[j]), ratio_interp[j]
   splog,/close

   ssfile = 'sphoto_resid_qso-'+run2d+'.ss'
   save, file=ssfile, wave, ratio


;; Read the LRGs on a single plate
;readspec, 4010+lindgen(20), zans=zans2, flux=flux2, invvar=invvar2, wave=wave2, $
; synflux=synflux2, plug=plug2, /align
;k=where(strmatch(zans2.class,'GAL*') AND zans2.zwarning EQ 0 $
; AND strmatch(plug2.objtype,'GAL*'))
;ratio2 = 0*wave2
;for i=0L, n_elements(wave2)-1 do begin $
; djs_iterstat, flux2[i,k]/synflux2[i,k], $
;  invvar=invvar2[i,k]*(synflux2[i,k]>0.1)^2*(synflux2[i,k] GT 0.1), mean=mn1 & ratio2[i]=mn1
;splot,wave,calibmed,yr=[0.9,1.1],xr=[3600,5000]
;soplot,wave2,ratio2,color='red'
;soplot,wave,calibmed,yr=[0.9,1.1],xr=[3600,5000]

stop
   return
end
