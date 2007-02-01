;+
; NAME:
;   readplugmap
;
; PURPOSE:
;   Read plugmap file and append tags as requested
;
; CALLING SEQUENCE:
;   plugmap = readplugmap( plugfile, [ plugdir=, /apotags, /deredden, $
;    exptime=, hdr=, /calibobj ] )
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
;   exptime   - Default exposure time SCI_EXPTIME to add to the output;
;               if there are multiple pointings and EXPTIME is set, then
;               the exposure time in each pointing is scaled such that
;               their sum is EXPTIME.
;   calibobj  - If set, then add a CALIBFLUX,CALIBFLUX_IVAR entries based upon
;               the calibObj files.  For stellar objects, this contains the
;               PSF fluxes in nMgy.  For galaxies, it contains the fiber fluxes
;               multiplied by the median (PSF/fiber) flux ratio for stars.
;               The MAG fields are left unchanged.
;               For objects with no calibObj entry, simply set these fields as:
;                 CALIBFLUX = 22.5 - 2.5*alog10(MAG), CALIBFLUX_IVAR = 0.
;               We apply the putative AB corrections to these fluxes
;               (but not to the MAG values).
;               Also add the SFD reddening values as SFD_EBV.
;
; OUTPUTS:
;   plugmap   - Plugmap structure
;
; OPTIONAL OUTPUTS:
;   hdr       - Header from Yanny-formatted plugmap file
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
;   dust_getval()
;   euler
;   plug2tsobj()
;   sdss_flagval()
;   splog
;   struct_addtags()
;   yanny_par()
;   yanny_read
;
; REVISION HISTORY:
;   29-Jan-2001  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function readplugmap, plugfile, plugdir=plugdir, $
 apotags=apotags, deredden=deredden, exptime=exptime, calibobj=calibobj, $
 hdr=hdr

   hdr = 0 ; Default return value

   ; The correction vector is here --- adjust this as necessary.
   ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
   correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

   thisfile = (findfile(djs_filepath(plugfile, root_dir=plugdir), $
    count=ct))[0]
   if (ct NE 1) then begin
      splog, 'WARNING: Cannot find plugmap file ' + plugfile
      return, 0
   endif

   yanny_read, thisfile, pstruct, hdr=hdr, stnames=stnames, /anonymous
   if (NOT keyword_set(pstruct)) then begin
      splog, 'WARNING: Invalid plugmap file ' + thisfile
      return, 0
   endif
   plugmap = *pstruct[(where(stnames EQ 'PLUGMAPOBJ'))[0]]
   nplug = n_elements(plugmap)

   ; Add the tags OFFSETID and SCI_EXPTIME for 
   plugmap = struct_addtags(plugmap, $
    replicate(create_struct('OFFSETID', 0L, 'SCI_EXPTIME', 0.), nplug))
   i = (where(stnames EQ 'PLUGMAPPOINT', ct))[0]
   if (ct GT 0) then begin
      splog, 'Using OFFSETID and SCI_EXPTIME from PLUGMAPPOINT structure'
      plugpoint = *pstruct[i]
      for j=0L, n_elements(plugpoint)-1L do begin
         k = where(abs(plugmap.xfocal - plugpoint[j].xfocal) LT 0.0001 $
          AND abs(plugmap.yfocal - plugpoint[j].yfocal) LT 0.0001, ct)
         if (ct GT 0) then begin
            plugmap[k[0]].offsetid = plugpoint[j].offsetid
            plugmap[k[0]].sci_exptime = plugpoint[j].sci_exptime
         endif
      endfor
   endif else begin
      ; Use default values
      plugmap.offsetid = 1
      sci_exptime = 1
   endelse
   if (keyword_set(exptime)) then begin
      iuniq = uniq(plugmap.offsetid, sort(plugmap.offsetid))
      exptot = total(plugmap[iuniq].sci_exptime)
      if (exptot GT 0) then begin
         splog, 'Rescaling SCI_EXPTIME values by ', exptime/exptot
         plugmap.sci_exptime = plugmap.sci_exptime * exptime/exptot
      endif
   endif

   plateid = (yanny_par(hdr, 'plateId'))[0]
   redden_med = yanny_par(hdr, 'reddeningMed')
   if (n_elements(redden_med) NE 5) then begin
      splog, 'WARNING: Wrong number of elements for reddeningMed'
      redden_med = fltarr(5)
   endif

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
       'CALIBFLUX', fltarr(5), $
       'CALIBFLUX_IVAR', fltarr(5), $
       'SFD_EBV', 0.), nplug)
      plugmap = struct_addtags(plugmap, addtags)

      iobj = where(strmatch(plugmap.holetype,'OBJECT*'))

      ;----------
      ; Read the SFD dust maps

      euler, plugmap[iobj].ra, plugmap[iobj].dec, ll, bb, 1
      plugmap[iobj].sfd_ebv = dust_getval(ll, bb, /interp)

      ;----------
      ; Attempt to read the calibObj photometry data

      tsobj = plug2tsobj(plateid, plugmap=plugmap[iobj], plugfile=plugfile)

      ; Do not use the calibObj structure if more than 10% of the non-sky
      ; objects do not have fluxes.
      if (keyword_set(tsobj)) then begin
         qexist = tsobj.psfflux[2] NE 0
         qsky = strmatch(plugmap[iobj].objtype,'SKY*')
         splog, 'Matched ', fix(total(qsky EQ 0 AND qexist)), $
          ' of ', fix(total(qsky EQ 0)), ' non-SKY objects'
         if (total((qsky EQ 0) AND qexist) LT 0.80*total(qsky EQ 0)) then begin
            splog, 'Discarding calibObj structure because < 80% matches'
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
      ; structure, simply translate the flux from the plugmap MAG values
      ; (as long as those values are in the range 0 < MAG < +50).

      for ifilt=0, 4 do begin
         ibad = where(plugmap[iobj].calibflux[ifilt] EQ 0 $
          AND plugmap[iobj].mag[ifilt] GT 0 $
          AND plugmap[iobj].mag[ifilt] LT 50, nbad)
         if (nbad GT 0) then begin
            splog, 'Using plug-map fluxes for ', nbad, $
             ' values in filter ', ifilt
            plugmap[iobj[ibad]].calibflux[ifilt] = $
             10.^((22.5 - plugmap[iobj[ibad]].mag[ifilt]) / 2.5)
            plugmap[iobj[ibad]].calibflux_ivar[ifilt] = 0
         endif
      endfor

      ;----------
      ; Apply AB corrections to the CALIBFLUX values (but not to MAG)

      factor = exp(-correction/2.5 * alog(10))
      for j=0,4 do plugmap.calibflux[j] = plugmap.calibflux[j] * factor[j]
      for j=0,4 do $
       plugmap.calibflux_ivar[j] = plugmap.calibflux_ivar[j] / factor[j]^2
   endif

   if (keyword_set(deredden)) then begin
      splog, 'Applying reddening vector ', redden_med
      iobj = where(strmatch(plugmap.holetype,'OBJECT*'))
      for ifilt=0, 4 do $
       plugmap[iobj].mag[ifilt] = plugmap[iobj].mag[ifilt] - redden_med[ifilt]
   endif

   return, plugmap
end
;------------------------------------------------------------------------------
