pro test_sbss_testplate

   plateid = 9001L
   tileid = 9001L
   setenv, 'HIPPARCOS_DIR=/u/schlegel/hipparcos'

   hipdat = hip_read(/small)
   tycdat = tyc_read(/small)
hipdat = [hipdat, tycdat]
   fname = '~schlegel/hoggpt/pro/hipparcos/plPlugMapT-9001.par'
   yanny_read, fname, allplug, hdr=hdr, enums=plugenum, structs=plugstruct
   ap = *allplug

   racen = dblarr(9)
   deccen = dblarr(9)
   for i=0, 8 do racen[i] = yanny_par(hdr, 'raCen'+string(i, form='(I1.1)'))
   for i=0, 8 do deccen[i] = yanny_par(hdr, 'decCen'+string(i, form='(I1.1)'))

   for ipoint=0, 8 do begin 
; rotate ra,dec to pointing ipoint and match to hipparcos
      plate_rotate, racen[0], deccen[0], racen[ipoint], deccen[ipoint], $
        ap.ra, ap.dec, ra, dec
      
      nmatch = djs_angle_match(ra, dec, hipdat.radeg, hipdat.dedeg, $
                               dtheta=1.5/3600.d, mcount=mcount, $
                               mindx=mindx, mdist=mdist, mmax=mmax)
      
      
      if nmatch GT 0 then begin
         wcoll = where(mcount GT 1, ncollide)
         if ncollide NE 0 then message, 'collision!!!'
         
         wmatch = where(mcount NE 0)
         nobj = total(ap[wmatch].holetype EQ 'OBJECT')
         npt = total(ap[wmatch].primtarget EQ ipoint)
         print, 'ipoint, nmatch, nobj, npt'
         print, ipoint, nmatch, nobj, npt
         
      endif
      
   endfor

  return
end
