;------------------------------------------------------------------------------
function linematch1, loglam, dloglam, zshift

   common com_linelist, linelist

   ;----------
   ; Read line lists and convert to vacuum

   if (NOT keyword_set(linelist)) then begin
      linefile = filepath('linematch.par', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      yanny_read, linefile, pdata
      linelist = *pdata[0]
      yanny_free, pdata
      vaclambda = linelist.lambda
      airtovac, vaclambda
      linelist.lambda = vaclambda
   endif

   ;----------
   ; Make a line list where the emission line wavelengths are redshifted
   ; to that of the galaxy.

   thislist = linelist.lambda * (1 + zshift * (linelist.type EQ 'EMISSION'))
;for i=0,n_elements(thislist)-1 do $
; print,linelist[i].type,linelist[i].lambda,thislist[i]

   ;----------
   ; Identify this feature with the nearest line.

   logdist = min(abs(loglam - alog10(thislist)), imin)
   if (logdist LE dloglam) then return, linelist[imin]

   ; Return no match
   retval = linelist[0]
   struct_assign, {junk:0}, retval ; Zero-out all elements
   return, retval
end
;------------------------------------------------------------------------------
function rline_matchpeaks, specstruct, pkstruct

   ;----------
   ; Loop through each possible line

   npix = n_elements(specstruct[0].loglam)
   xvec = findgen(npix)

   dloglam = 4.0e-4 ; Sky+object line exclusion width of 4 pixels

   npeak = n_elements(pkstruct)
   for ipeak=0, npeak-1 do begin
      thisspec = specstruct[ pkstruct[ipeak].fiber ]
      linterp, xvec, thisspec.loglam, pkstruct[ipeak].x, thisloglam
      thisline = linematch1(thisloglam, dloglam, thisspec.zans.z)
      if (ipeak EQ 0) then allline = thisline $
       else allline = [allline, thisline]
   endfor

   return, allline
end
;------------------------------------------------------------------------------
