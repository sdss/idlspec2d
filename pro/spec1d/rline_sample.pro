; This is a hacked version of RLINE_LOOP, for just getting the full
; sample selection of BRG's.
pro rline_sample, plate=plate, mjd=mjd

   primtarget = 32 ; BRG's

   ;----------
   ; Get a list of plates

   if (keyword_set(plate)) then begin
      nplate = n_elements(plate)
      if (NOT keyword_set(mjd)) then mjd = lonarr(nplate) ; zeros
      plist = replicate(create_struct('plateid', 0L, 'mjd', 0L), nplate)
      plist.plateid = plate
      plist.mjd = mjd
   endif else begin
      ; If PLATE is not specified, then determine a list of good plates...
      platelist, plist=plist
      ii = where(plist.nums[0] GT 0 $
;       AND plist.plateid GE 300 AND plist.plateid LE 349 $ ; ????
       AND plist.snvec[0] GT 15 AND plist.snvec[1] GT 15 $
       AND plist.snvec[2] GT 15 AND plist.snvec[3] GT 15)
      plist = plist[ii]
   endelse

   ;----------
   ; Loop through each plate

   nplate = n_elements(plist)
   for iplate=0, nplate-1 do begin

      splog, 'WORKING ON PLATE ', plist[iplate].plateid, $
       ' MJD ', plist[iplate].mjd
      spec = rline_getplate(plist[iplate].plateid, mjd=plist[iplate].mjd, $
       primtarget=primtarget, class='GALAXY', maxrchi2=2.0, /quick)

      if (keyword_set(spec)) then $
       sampzans = struct_append(sampzans, spec.zans)
      if (keyword_set(spec)) then $
       sampplug = struct_append(sampplug, spec.plug)

   endfor

save,file='rline_sample.ss'
stop

end
