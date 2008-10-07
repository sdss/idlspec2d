;+
; NAME:
;   plug2tsobj
;
; PURPOSE:
;   Construct a tsObj structure that matches all entries in a plugmap structure
;
; CALLING SEQUENCE:
;   tsobj = plug2tsobj(plateid, [ra, dec, plugmap=, dmin=, /silent ])
;
; INPUTS:
;   plateid    - Plate number; this can be either a scalar, in which case
;                the same plate is used for all objects, or a vector.
;
; OPTIONAL INPUTS:
;   ra         - Array of right ascension (degrees)
;   dec        - Array of declination (degrees)
;   plugmap    - Plug map structure, which must contain RA, DEC.
;                This must be set if RA and DEC are not set.
;   dmin       - Minimum separation between input position and position
;                of the closest match using MATCHRA,MATCHDEC in the calibObj
;                files; default to 1.0 arcsec.
;
; OUTPUTS:
;   tsobj      - tsObj structure, sorted such that each entry corresponds
;                to each entry in the PLUGMAP structure; return 0 if the
;                tsObj file was not found on disk.
;
; OPTIONAL OUTPUTS:
;   silent     - Make the call to MRDFITS silent, but still log warning
;                messages if the calibObj file is not found or is empty.
;
; COMMENTS:
;   The calibObj files are assumed to be in the directory $SPECTRO_DATA/calibobj.
;   These files were to have only the 640 objects for each plate.
;   But since plates can be re-plugged, we must re-sort these
;   files to match the object ordering in the plug-map structure.
;
;   If the TSOBJMAPROOT environment variable is set, then attempt to
;   read the survey tsObj-$PLATE files instead of the calibPlateP
;   files. 
;
;   Print a warning message for non-matched objects only if they are not skies,
;   according to the PLUGMAP structure.
;
; EXAMPLES:
;   Read the plug-map for plate 306, fibers 1 to 10, then construct the
;   calibObj structure:
;   > readspec, 306, indgen(10)+1, plug=plug
;   > tsobj = plug2tsobj(306,plugmap=plug)
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   fits_read
;   mrdfits
;   splog
;
; REVISION HISTORY:
;   25-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function plug2tsobj, plateid, ra1, dec1, plugmap=plugmap, dmin=dmin, $
 silent=silent

    root_dir = getenv('SPECTRO_DATA')
   tsroot_dir = getenv('TSOBJMAPROOT')
   if (NOT keyword_set(root_dir) and NOT keyword_set(tsroot_dir)) then $
    message, 'Environment variable SPECTRO_DATA or TSROOT_DIR must be set!'

   if (keyword_set(plugmap)) then begin
      if (n_elements(ra1) EQ 0) then ra = plugmap.ra
      if (n_elements(dec1) EQ 0) then dec = plugmap.dec
   endif
   if (n_elements(ra) EQ 0) then ra = ra1
   if (n_elements(dec) EQ 0) then dec = dec1

   if (NOT keyword_set(dmin)) then dmin = 2.0

   ;----------
   ; If PLATEID is a vector, then sort by plate number and call this routine
   ; iteratively, once for each plate number.

   if (n_elements(plateid) GT 1) then begin
      platenums = plateid[ uniq(plateid, sort(plateid)) ]
      nplate = n_elements(platenums)
      for iplate=0, nplate-1 do begin
         indx = where(plateid EQ platenums[iplate])
         tsobj1 = plug2tsobj( platenums[iplate], ra[indx], dec[indx] )

         if (iplate EQ 0) then begin
            tsobj = replicate(tsobj1[0], n_elements(ra))
            tmpobj = tsobj[0]
         endif

         for i=0, n_elements(tsobj1)-1 do begin
            struct_assign, tsobj1[i], tmpobj
            tsobj[indx[i]] = tmpobj
         endfor
      endfor
      return, tsobj
   endif
   
                                ; Prefer the tsObj files. 
   sortedmap = 0
   if keyword_set(tsroot_dir) then begin
      platestr = string(plateid[0],format='(i4.4)')
      filename = 'tsObj-' + platestr + '.fit*'
      filename = (findfile(filepath(filename, root_dir=tsroot_dir)))[0]
      if (keyword_set(filename)) then begin
         sortedmap = 1
      endif else begin
         filename = 'tsObj*-*' + platestr + '.fit*'
         filename = (findfile(filepath(filename, root_dir=root_dir)))[0]
      end
      hdu = 2                   ; Go for BEST.
   end else begin
      platestr = string(plateid[0],format='(i4.4)')
      filename = 'calibPlateP-' + platestr + '.fits'
      filename = (findfile(filepath(filename, root_dir=root_dir, $
                                    subdirectory='calibobj')))[0]
      hdu = 1
   end

   if (NOT keyword_set(filename)) then begin
      splog, 'WARNING: calibObj file not found for plate ' + platestr
      return, 0
   endif
   splog,'reading object calibration file: ', filename

   ; Make certain that the file exists and is valid
   message = 0
   fits_read, filename, junk, /no_abort, message=message
   if (keyword_set(message)) then tstemp = 0 $ ; File is invalid FITS file
    else tstemp = mrdfits(filename, hdu, silent=silent)
   if (NOT keyword_set(tstemp)) then begin
      splog, 'WARNING: calibObj file is empty: ' + filename
      return, 0
   endif

   tsobj1 = tstemp[0]
   struct_assign, {junk:0}, tsobj1 ; Zero-out this new structure
   tsobj = replicate(tsobj1, n_elements(ra))

   ;----------
   ; Find the tsObj-file entry for each plug-map entry by matching
   ; the RA,DEC positions on the sky.  Insist that they agree to 2 arcsec.

   if sortedmap then begin
                                ; The plugmap and the tsObj file
                                ; should match 1-to-1. But I'm
                                ; not certain about how, so will keep
                                ; the spherematch
      pmi1 = intarr(n_elements(ra))
      pmi2 = intarr(n_elements(ra))
      spherematch, ra, dec, tstemp.ra, tstemp.dec, dmin/3600, m1,m2,d12
      pmi1[m1] = 1
      pmi2[m2] = 1
      tsobj[m1] = tstemp[m2]

      duds = where(pmi1 eq 0 and not strmatch(plugmap.objtype,'SKY*'), dudcnt)
      for i=0,dudcnt-1 do begin
         splog, 'WARNING: Unmatched fiber ', strtrim(plugmap[duds[i]].objtype,2), duds[i],  $
                + ' (', ra[duds[i]], dec[duds[i]], ') r=', plugmap[duds[i]].mag[2]
      end
      duds = where(pmi2 eq 0 and tstemp.ra NE 0 and tstemp.dec NE 0, dudcnt)
      for i=0,dudcnt-1 do begin
         splog, 'WARNING: Unmatched non-0 calibration id ', duds[i],  $
                + ' at (', (tstemp.ra)[duds[i]], (tstemp.dec)[duds[i]], ')'
      end
   end else begin
      for iplug=0, n_elements(ra)-1 do begin
                                ; Assume that this object is non-existent if RA=DEC=0
         if (ra[iplug] NE 0 AND dec[iplug] NE 0) then begin
            adist = djs_diff_angle(tstemp.matchra, tstemp.matchdec, $
                                   ra[iplug], dec[iplug])
            thismin = min(adist, imin)
            if (thismin GT dmin/3600.) then begin
               if (keyword_set(plugmap)) then begin
                  if (strmatch(plugmap[iplug].objtype,'SKY*') EQ 0) then $
                     splog, 'Warning: Unmatched OBJTYPE=' $
                            + strtrim(plugmap[iplug].objtype,2) $
                            + ' at RA=', ra[iplug], ',DEC=', dec[iplug]
               endif else begin
                  splog, 'Warning: Unmatched ' $
                         + ' at RA=', ra[iplug], ',DEC=', dec[iplug]
               endelse
            endif else begin
               tsobj[iplug] = tstemp[imin]
            endelse
         endif
      endfor
   end
                                ; Extract the columns from Survey
                                ; tsObj-$plate files which we need for
                                ; flux calibration. Explicity do not
                                ; return any other columns: I want any
                                ; other consumers to _explode_.
   if keyword_set(tsroot_dir) then begin
      mytsobj = replicate(create_struct('PSFFLUX', fltarr(5), $
                                        'PSFFLUX_IVAR', fltarr(5), $
                                        'FIBERFLUX', fltarr(5), $
                                        'FIBERFLUX_IVAR', fltarr(5), $
                                        'FLAGS', lonarr(5), $
                                        'FLAGS2', lonarr(5)), n_elements(ra))
      
      for f=0,4 do begin
                                ; PSF fluxes first, then fiberfluxes.
                                ; -9999s have been mentioned. [?]
         duds = where(tsobj.psfcounts[f] lt -9998 or $
                tsobj.psfcountserr[f] lt -9998 or $
                tsobj.psfcounts[f] eq 0 or $
                tsobj.psfcountserr[f] eq 0, dudcnt, $
                complement=good, ncomplement=goodcnt)

         if goodcnt gt 0 then begin
            mytsobj[good].psfflux[f] = 10^((tsobj[good].psfcounts[f] - 22.5)/(-2.5))
            mytsobj[good].psfflux_ivar[f] = $
               1./(mytsobj[good].psfflux[f] * (1-10^(tsobj[good].psfcountserr[f]/(-2.5))))^2
         end

         if dudcnt gt 0 then begin
            mytsobj[duds].psfflux[f] = 0.
            mytsobj[duds].psfflux_ivar[f] = 0.
         end

         duds = where(tsobj.fibercounts[f] lt -9998 or tsobj.fibercountserr[f] lt -9998 or $
                      tsobj.fibercounts[f] eq 0 or tsobj.fibercountserr[f] eq 0, $
                      complement=good, ncomplement=goodcnt)

         if goodcnt gt 0 then begin
            mytsobj[good].fiberflux[f] = 10^((tsobj[good].fibercounts[f] - 22.5)/(-2.5))
            mytsobj[good].fiberflux_ivar[f] = $
               1./(mytsobj[good].fiberflux[f] * (1-10^(tsobj[good].fibercountserr[f]/(-2.5))))^2
         end

         if dudcnt gt 0 then begin
            mytsobj[duds].fiberflux[f] = 0.
            mytsobj[duds].fiberflux_ivar[f] = 0.
         end

         mytsobj.flags = tsobj.flags
         mytsobj.flags2 = tsobj.flags2
      end
      
      return, mytsobj
   end

   return, tsobj
end
;------------------------------------------------------------------------------
