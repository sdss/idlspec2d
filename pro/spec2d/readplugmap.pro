;+
; NAME:
;   readplugmap
;
; PURPOSE:
;   Read plugmap file and append tags as requested
;
; CALLING SEQUENCE:
;   plugmap = readplugmap( plugfile, [ plugdir=, /apotags, /deredden, $
;    /calibobj ] )
;
; INPUTS:
;   plugfile  - Name of Yanny-parameter plugmap file
;
; OPTIONAL KEYWORDS:
;   plugdir   - Directory for PLUGFILE
;   apotags   - If set, then add a number of tags to the output structure
;               constructed from the Yanny header.  These tags are:
;               CARTID, PLATEID, TILEID, RAPLATE, DECPLATE, REDDEN_MED.
;               Also add the tags FIBERSN[3], SYTHMAG[3] which are used
;               by the on-the-mountain reductions.
;   deredden  - If set, then deredden the MAG fields using the median
;               reddening value for the entire plate as found in the
;               Yanny header of the plugmap file; this is done for the
;               on-the-mountain reductions.
;   calibobj  - If set, then add a CALIBFLUX,CALIBFLUX_IVAR entries based upon
;               the calibObj files.  For stellar objects, this contains the
;               PSF fluxes in nMgy.  For galaxies, it contains the fiber fluxes
;               multiplied by the median (PSF/fiber) flux ratio for stars.
;               The MAG fields are left unchanged.
;               For objects with no calibObj entry, simply set these fields as:
;                 CALIBFLUX = 22.5 - 2.5*alog10(MAG), CALIBFLUX_IVAR = 0.
;               We apply the putative AB corrections to these fluxes
;               (but not to the MAG values).
;
; OUTPUTS:
;   plugmap   - Plugmap structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Do not use the calibObj structure if more than 10% of the non-sky
;   objects do not have fluxes.
;
; EXAMPLES:
;
; BUGS:
;   The AB corrections are hard-wired to be the same as in the photoop
;   product as of 18 Feb 2004.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   splog
;   struct_addtags()
;   yanny_readone()
;
; REVISION HISTORY:
;   29-Jan-2001  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function readplugmap, plugfile, plugdir=plugdir, $
 apotags=apotags, deredden=deredden, calibobj=calibobj

   ; The correction vector is here --- adjust this as necessary.
   ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
   correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

   thisfile = (findfile(djs_filepath(plugfile, root_dir=plugdir), $
    count=ct))[0]
   if (ct NE 1) then begin
      print, 'WARNING: Cannot find plugmap file ' + plugfile
      return, 0
   endif

   plugmap = yanny_readone(thisfile, 'PLUGMAPOBJ', hdr=hdr, /anonymous)
   if (NOT keyword_set(plugmap)) then begin
      print, 'WARNING: Invalid plugmap file ' + thisfile
      return, 0
   endif

   nplug = n_elements(plugmap)
   plateid = (yanny_par(hdr, 'plateId'))[0]
   redden_med = yanny_par(hdr, 'reddeningMed')

   if (keyword_set(apotags)) then begin
      addtags = { $
       cartid   : long((yanny_par(hdr, 'cartridgeId'))[0]), $
       plateid  : long(plateid), $
       tileid   : long((yanny_par(hdr, 'tileId'))[0]), $
       raplate  : float((yanny_par(hdr, 'raCen'))[0]), $
       decplate : float((yanny_par(hdr, 'decCen'))[0]), $
       redden_med : float(redden_med), $
       fibersn    : fltarr(3), $
       synthmag   : fltarr(3) }
      plugmap = struct_addtags(plugmap, replicate(addtags, nplug))
   endif

   if (keyword_set(calibobj)) then begin
      splog, 'Adding fields from calibObj file'
      addtags = replicate(create_struct( $
       'CALIBFLUX', fltarr(5), 'CALIBFLUX_IVAR', fltarr(5)), nplug)
      plugmap = struct_addtags(plugmap, addtags)

      iobj = where(strmatch(plugmap.holetype,'OBJECT*'))
      tsobj = plug2tsobj(plateid, plugmap=plugmap[iobj])

      ; Do not use the calibObj structure if more than 10% of the non-sky
      ; objects do not have fluxes.
      if (keyword_set(tsobj)) then begin
         qexist = tsobj.psfflux[2] NE 0
         qsky = strmatch(plugmap[iobj].objtype,'SKY*')
         if (total((qsky EQ 0) AND qexist) LT 0.9*total(sky EQ 0)) then begin
            splog, 'Discarding calibObj structure because < 90% matches'
            tsobj = 0
         endif
      endif

      if (keyword_set(tsobj)) then begin

         ; Assume that all objects not called a 'GALAXY' are stellar objects
         qstar = strmatch(plugmap[iobj].objtype, 'GALAXY*') EQ 0
         istar = where(qstar AND qexist, nstar)
         igal = where(qstar EQ 0 AND qexist, ngal)
         pratio = fltarr(5) + 1
         if (nstar GT 0) then begin
            plugmap[iobj[istar]].calibflux = tsobj[istar].psfflux
            plugmap[iobj[istar]].calibflux_ivar = tsobj[istar].psfflux_ivar
            ; Compute the ratio of PSF/FIBER flux for stars in each filter,
            ; using only stars that are brighter than 30 nMgy (= 18.8 mag).
            ; If no such stars, then this ratio is set to unity.
            for ifilt=0, 4 do begin
               v1 = tsobj[istar].psfflux[ifilt]
               v2 = tsobj[istar].fiberflux[ifilt]
               jj = where(v1 GT 30 AND v2 GT 30, ct)
               if (ct GT 0) then pratio[ifilt] = median([ v1[jj] / v2[jj] ])
            endfor

         endif
         splog, 'PSF/fiber flux ratios = ', pratio
         if (ngal GT 0) then begin
            for ifilt=0, 4 do begin
               plugmap[iobj[igal]].calibflux[ifilt] = $
                tsobj[igal].fiberflux[ifilt] * pratio[ifilt]
               plugmap[iobj[igal]].calibflux_ivar[ifilt] = $
                tsobj[igal].fiberflux_ivar[ifilt] / (pratio[ifilt])^2
            endfor
         endif

         ; Reject any fluxes based upon suspect PHOTO measurements,
         ; as indicated by the PHOTO flags.
         badbits2 = sdss_flagval('OBJECT2','SATUR_CENTER') $
          OR sdss_flagval('OBJECT2','INTERP_CENTER') $
          OR sdss_flagval('OBJECT2','PSF_FLUX_INTERP')
         qgoodphot = (tsobj.flags2 AND badbits2) EQ 0
         plugmap[iobj].calibflux = plugmap[iobj].calibflux * qgoodphot
         plugmap[iobj].calibflux_ivar = plugmap[iobj].calibflux_ivar * qgoodphot
      endif else begin
         splog, 'WARNING: No calibObj structure found for plate ', plateid
      endelse

      ;----------
      ; For any objects that do not have photometry from the calibObj
      ; structure, simply translate the flux from the plugmap MAG values.

      ibad = where(plugmap[iobj].calibflux EQ 0, nbad)
      if (nbad GT 0) then begin
         splog, 'Using plug-map fluxes for ', nbad, ' objects'
         plugmap[iobj[ibad]].calibflux = $
          10.^((22.5 - plugmap[iobj[ibad]].mag) / 2.5)
      endif

      ;----------
      ; Apply AB corrections to the CALIBFLUX values (but not to MAG)

      factor = exp(-correction/2.5 * alog(10))
      for j=0,4 do plugmap.calibflux[j] = plugmap.calibflux[j] * factor
      for j=0,4 do $
       plugmap.calibflux_ivar[j] = plugmap.calibflux_ivar[j] / factor^2
   endif

   if (keyword_set(deredden)) then begin
      splog, 'Applying reddening vector ', redden_med
      for ifilt=0, 4 do $
       plugmap[iobj].mag[ifilt] = plugmap[iobj].mag[ifilt] - redden_med[ifilt]
   endif

   return, plugmap
end
;------------------------------------------------------------------------------
