pro elodie_best, objflux, objivar, objloglam, objdloglam

; ???
readspec,406,zans=zans
i=where(strtrim(zans.class,2) EQ 'STAR' $
 AND strtrim(zans.subclass,2) EQ 'F' $
 AND zans.zwarning EQ 0 AND zans.sn_median GT 10)
zans=zans[i]
readspec,406,zans.fiberid,flux=objflux,invvar=objivar,loglam=loglam
objloglam0 = loglam[0]
objdloglam = loglam[1] - loglam[0]
hdr = ['COEFF0  = ' + string(objloglam0), $
       'COEFF1  = ' + string(objdloglam) ]

zmin = -0.0015
zmax = 0.0015

   stime0 = systime(1)

   ;----------
   ; Read all the Elodie spectra

   elodie_path = getenv('ELODIE_PATH')
   allfiles = findfile(filepath('00*', root_dir=elodie_path), count=nstar)
; ???
nstar=25
allfiles = allfiles[0:nstar-1]
   t0 = systime(1)
   starhdr = replicate(ptr_new(), nstar)
   for istar=0, nstar-1 do begin
      splog, 'Reading file ', istar+1, ' of ', nstar
      thisflux = read_elodie(allfiles[istar], loglam=starloglam, hdr=hdr)
      if (NOT keyword_set(starflux)) then starflux = thisflux $
       else starflux = [[starflux],[thisflux]]
      starhdr[istar] = ptr_new(hdr)
   endfor
   npix = n_elements(starloglam)
   splog, 'Time to read all files = ', systime(1) - t0
stop

   ;----------
   ; Trim wavelengths to those covered by the majority of the objects

   fracgpix = total(starflux NE 0, 2) / nstar
   igood = where(fracgpix GT 0.95)
   i1 = igood[0]
   i2 = (reverse(igood))[0]
   starloglam = starloglam[i1:i2]
   starflux = starflux[i1:i2,*]

   ;----------
   ; Interpolate over bad data, of which there is very little

   starflux = djs_maskinterp(starflux, starflux EQ 0, /const, iaxis=0)

   ;----------
   ; Compute the best-fit between each object and each star

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]

   for istar=0, nstar-1 do begin
splog, 'Star number ', istar, ' of ', nstar
      res1 = zfind(objflux, objivar, hdr=hdr, starflux=starflux[*,istar], $
       starloglam0=starloglam[0], npoly=3, zmin=zmin, zmax=zmax)
      if (istar EQ 0) then res_all = replicate(res1[0], nobj, nstar)
      res_all[*,istar] = res1[*]
   endfor

   ;----------
   ; For each object, select the best-fit star

   res_best = replicate(res1[0], nobj)
   for iobj=0, nobj-1 do begin
      junk = min(res_all[iobj,*].rchi2, imin)
      res_best[iobj] = res_all[iobj,imin]
   endfor
; Get this info from the header ???
; FILENAME,OBJECT,TEFF,LOGG,SPTYPEB-V,HIERARCH [Fe/H]

   splog, 'Total time for ELODIE_BEST = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'

plot,zans.z*3e5,res_best.z*3e5,ps=7
for i=0,nstar-1 do oplot,zans.z*3e5,res_all[*,i].z*3e5,ps=3
oplot,[-200,200],[-200,200]

stop

junk=poly_array(2172,3)
synflux=fltarr(2172,nobj)
for i=0,26 do synflux[*,i]=junk#res_best[i].theta[1:3]
plot,djs_median(abs(synflux),1)/djs_median(objflux,1)

   return
end
