pro rline_loop, plate=plate, mjd=mjd

   primtarget = 32 ; BRG's

   ;----------
   ; Get a list of plates

   if (keyword_set(plate)) then begin
      nplate = n_elements(plate)
      if (NOT keyword_set(mjd)) then mjd = lonarr(nplate) ; zeros
      plist = replicate(create_struct('plate', 0L, 'mjd', 0L), nplate)
      plist.plate = plate
      plist.mjd = mjd
   endif else begin
      ; If PLATE is not specified, then determine a list of good plates...
      platelist, plist=plist
      ii = where(plist.nums[0] GT 0 $
;       AND plist.plate GE 300 AND plist.plate LE 349 $ ; ????
       AND plist.snvec[0] GT 15 AND plist.snvec[1] GT 15 $
       AND plist.snvec[2] GT 15 AND plist.snvec[3] GT 15)
      plist = plist[ii]
   endelse

   ;----------
   ; Loop through each plate

   nplate = n_elements(plist)
   for iplate=0, nplate-1 do begin

      splog, 'WORKING ON PLATE ', plist[iplate].plate, $
       ' MJD ', plist[iplate].mjd
      spec = rline_getplate(plist[iplate].plate, mjd=plist[iplate].mjd, $
       primtarget=primtarget, class='GALAXY', maxrchi2=2.0)

      if (keyword_set(spec)) then begin

         ; Accumulate ZANS and PLUG for the full sample
         sampzans = struct_append(sampzans, spec.zans)
         sampplug = struct_append(sampplug, spec.plug)

         rline_findew, spec,  ew, ewinv, fopt, finv
         pks = rline_findpeaks(ew, ewinv)

         ; Convert the peak position from pixel to log-wavelength
         npix = n_elements(spec[0].loglam)
         xvec = findgen(npix)
         linterp, xvec, spec[0].loglam, pks.x, xloglam

         if (keyword_set(pks)) then begin
            lines = rline_matchpeaks(spec, pks)
            for ipk=0, n_elements(lines)-1 do $
             allpks = struct_append( allpks, $
              create_struct('zans', spec[pks[ipk].fiber].zans, $
                            'plug', spec[pks[ipk].fiber].plug, $
                            pks[ipk], $
                            'xloglam', xloglam[ipk], $
                            lines[ipk]) )

         endif
      endif
   endfor

save,file='rline.ss'
stop

j = where(allpks.lambda EQ 0 AND allpks.xloglam LT alog10(7500))
goodpks = allpks[j]
readspec, goodpks.zans.plate, goodpks.zans.fiberid, mjd=goodpks.zans.mjd, $
 flux=flux, flerr=flerr, invvar=invvar, andmask=andmask, ormask=ormask, $
 loglam=loglam, wave=wave, zans=zans
synflux = 0 * flux
for i=0, n_elements(goodpks)-1 do $
 synflux[*,i] = synthspec(goodpks[i].zans, loglam=loglam[*,i])
save,file='goodpks.ss'

restore,'goodpks.ss'
k=1
thiswave = 10^goodpks[k].xloglam
ii=where(wave[*,k] NE 0)
splot,wave[ii,k],flux[ii,k], xrange=thiswave+[-200,200]
soplot,wave[ii,k],synflux[ii,k],color='blue'
soplot,wave[ii,k],flerr[ii,k],color='red'
soplot,[thiswave,thiswave],!y.crange,color='green


j = where(allpks.lambda EQ 0 AND allpks.xloglam LT alog10(7500))
splot,10^allpks.xloglam,allpks.sn,ps=3
soplot,10^allpks[j].xloglam,allpks[j].sn,ps=3,color='red'
k = j[where(allpks[j].sn GT 20)]

jj=j[0]
plotspec,allpks[jj].zans.plate,allpks[jj].zans.fiberid
soplot,10^allpks[jj].xloglam,1,ps=4,charsize=4,color='green'

end
