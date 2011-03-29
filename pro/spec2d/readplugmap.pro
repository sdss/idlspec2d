;+
; NAME:
;   readplugmap
;
; PURPOSE:
;   Read plugmap file and append tags as requested
;
; CALLING SEQUENCE:
;   plugmap = readplugmap( plugfile, [ spectrographid, plugdir=, $
;    /apotags, /deredden, /calibobj, exptime=, hdr=, fibermask=, _EXTRA= ] )
;
; INPUTS:
;   plugfile  - Name of Yanny-parameter plugmap file
;
; OPTIONAL INPUTS:
;   spectrographid  - The spectrograph number, either 1 or 2;
;                     if not set (or 0), then return all object fibers
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
;               the calibObj (deprecated) or photoPlate files.
;               For stellar objects, this contains the
;               PSF fluxes in nMgy.  For galaxies, it contains the fiber fluxes
;               multiplied by the median (PSF/fiber) flux ratio for stars.
;               The MAG fields are left unchanged.
;               For objects with no calibObj entry, simply set these fields as:
;                 CALIBFLUX = 22.5 - 2.5*alog10(MAG), CALIBFLUX_IVAR = 0.
;               We apply the putative AB corrections to these fluxes
;               (but not to the MAG values).
;               Also add the SFD reddening values as SFD_EBV.
;   _EXTRA    - Keywords for PLUG2TSOBJ(), such as MJD,INDIR
;
; OUTPUTS:
;   plugmap   - Plugmap structure
;
; OPTIONAL OUTPUTS:
;   hdr       - Header from Yanny-formatted plugmap file
;   fibermask - Byte array with bits set for unknown fibers
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
function readplugmap_sort, plugmap, fibermask=fibermask

   qobj = strmatch(plugmap.holetype,'OBJECT')
   indx = where(qobj, nfiber)
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) $
   else if (n_elements(fibermask) NE nfiber) then $
    message, 'Number of elements in FIBERMASK do not match NFIBER'

   blankmap = plugmap[0]
   struct_assign, {junk:0}, blankmap
   plugsort = replicate(blankmap, nfiber)
   plugsort.holetype = 'OBJECT'
   plugsort.objtype = 'NA'
   plugsort.fiberid = -1

   igood = where(qobj AND plugmap.fiberid GT 0, ngood)
   if (ngood EQ 0) then $
    message, 'No fibers found in plugmap!'
   iplace = plugmap[igood].fiberid - 1
   plugsort[iplace] = plugmap[igood]

   ; Set the appropriate fibermask bit if a fiber not found in plugmap file.
   ; Do this by first setting all bits to 1, then unsetting the good ones.
   fibermask = fibermask OR fibermask_bits('NOPLUG')
   fibermask[iplace] = fibermask[iplace] - fibermask_bits('NOPLUG')

   ; Fill in unplugged fibers with arbitrary entries, and assign
   ; them a FIBERID.  After this, plugsort.fiberid should run from 1...nfiber
   imissing = where(plugsort.fiberid LE 0, nmissing)
   splog, 'Number of missing fibers: ', nmissing
   if (nmissing GT 0) then begin
      ifill = where(qobj AND plugmap.fiberid LE 0, nfill)
      plugsort[imissing] = plugmap[ifill]
      plugsort[imissing].fiberid = imissing + 1
   endif

   return, plugsort
end
;------------------------------------------------------------------------------
function readplugmap, plugfile, spectrographid, plugdir=plugdir, $
 apotags=apotags, deredden=deredden, exptime=exptime, calibobj=calibobj, $
 hdr=hdr, fibermask=fibermask, _EXTRA=KeywordsForPhoto

   hdr = 0 ; Default return value
   if (keyword_set(fibermask)) then $
    message, 'FIBERMASK is already set!'

   ; The correction vector is here --- adjust this as necessary.
   ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
   correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

   ;----------
   ; Read the file

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

   plugmap.ra = (360d0 + plugmap.ra) MOD 360d0

   ;----------
   ; Trim to object fibers only, sort them, and trim to spectrographid

   plugmap = readplugmap_sort(plugmap, fibermask=fibermask)

   ;----------
   ; Add the tags OFFSETID and SCI_EXPTIME for 

   plugmap = struct_addtags(plugmap, $
    replicate(create_struct('OFFSETID', 0L, 'SCI_EXPTIME', 0.), $
    n_elements(plugmap)))
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

   ;----------
   ; Append some information from the plateHoles file

   platelist_dir = getenv('PLATELIST_DIR')
   platefile = 'plateHoles-' + string(plateid,format='(i6.6)') + '.par'
   if (keyword_set(platelist_dir)) then begin
      thisfile = (findfile(djs_filepath(platefile, $
       root_dir=platelist_dir, subdir=['plates','*','*']), count=ct))[0]
      if (ct GT 0) then begin
         plateholes = yanny_readone(thisfile, /anonymous)
         iobj = where(strmatch(plateholes.holetype,'BOSS*'))
         plateholes = plateholes[iobj]
         isort = lonarr(n_elements(plugmap)) - 1
         for i=0L, n_elements(plugmap)-1 do $
          isort[i] = where(plateholes.xfocal EQ plugmap[i].xfocal $
           AND plateholes.yfocal EQ plugmap[i].yfocal)
         plateholes = plateholes[isort>0]
         blankhole = plateholes[0]
         struct_assign, {junk: 0}, blankhole
         ibad = where(iobj EQ -1, nbad)
         for i=0L, nbad-1 do plateholes[ibad[i]] = blankhole
         htags = ['SOURCETYPE','LAMBDA_EFF','ZOFFSET','BLUEFIBER', $
          'BOSS_TARGET*','ANCILLARY_TARGET*', $
          'RUN','RERUN','CAMCOL','FIELD','ID']
         plugmap = struct_addtags(plugmap, $
          struct_selecttags(plateholes, select_tags=htags))
      endif
   endif

   ;----------
   ; Optionally add tags for SOS

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
      plugmap = struct_addtags(plugmap, replicate(addtags, n_elements(plugmap)))
   endif

   ;----------
   ; Read calibObj or photoPlate photometry data

   if (keyword_set(calibobj)) then begin
      splog, 'Adding fields from calibObj file'
      addtags = replicate(create_struct( $
       'CALIBFLUX', fltarr(5), $
       'CALIBFLUX_IVAR', fltarr(5), $
       'SFD_EBV', 0.), n_elements(plugmap))
      plugmap = struct_addtags(plugmap, addtags)

      ;----------
      ; Read the SFD dust maps

      euler, plugmap.ra, plugmap.dec, ll, bb, 1
      plugmap.sfd_ebv = dust_getval(ll, bb, /interp)

      ;----------
      ; Attempt to read the calibObj photometry data

      tsobj = plug2tsobj(plateid, _EXTRA=KeywordsForPhoto)

      ; Do not use the calibObj structure if more than 20% of the non-sky
      ; objects do not have fluxes.
      if (keyword_set(tsobj)) then begin
         qexist = tsobj.psfflux[2] NE 0
         qsky = strmatch(plugmap.objtype,'SKY*')
         splog, 'Matched ', fix(total(qsky EQ 0 AND qexist)), $
          ' of ', fix(total(qsky EQ 0)), ' non-SKY objects'
         if (total((qsky EQ 0) AND qexist) LT 0.80*total(qsky EQ 0)) then begin
            splog, 'Discarding calibObj structure because < 80% matches'
            tsobj = 0
         endif
      endif

      if (keyword_set(tsobj)) then begin

         ; Assume that all objects not called a 'GALAXY' are stellar objects
         qstar = strmatch(plugmap.objtype, 'GALAXY*') EQ 0
         istar = where(qstar AND qexist, nstar)
         igal = where(qstar EQ 0 AND qexist, ngal)
         if (tag_exist(tsobj,'FIBER2FLUX')) then begin
            fiberflux = transpose(tsobj.fiber2flux)
            fiberflux_ivar = transpose(tsobj.fiber2flux_ivar)
            pratio = [2.085, 2.085, 2.116, 2.134, 3.135]
         endif else begin
            fiberflux = transpose(tsobj.fiberflux)
            fiberflux_ivar = transpose(tsobj.fiberflux_ivar)
            pratio = [1.343, 1.336, 1.354, 1.363, 1.367]
         endelse
         if (nstar GT 0) then begin
            plugmap[istar].calibflux = tsobj[istar].psfflux
            plugmap[istar].calibflux_ivar = tsobj[istar].psfflux_ivar
;            ; Compute the ratio of PSF/FIBER flux for stars in each filter,
;            ; using only stars that are brighter than 30 nMgy (= 18.8 mag).
;            ; If no such stars, then this ratio is set to unity.
;            for ifilt=0, 4 do begin
;               v1 = tsobj[istar].psfflux[ifilt]
;               v2 = fiberflux[istar,ifilt]
;               jj = where(v1 GT 30 AND v2 GT 30, ct)
;               if (ct GT 0) then pratio[ifilt] = median([ v1[jj] / v2[jj] ])
;            endfor
         endif
         splog, 'PSF/fiber flux ratios = ', pratio
         if (ngal GT 0) then begin
            for ifilt=0, 4 do begin
               plugmap[igal].calibflux[ifilt] = $
                fiberflux[igal,ifilt] * pratio[ifilt]
               plugmap[igal].calibflux_ivar[ifilt] = $
                fiberflux_ivar[igal,ifilt] / (pratio[ifilt])^2
            endfor
         endif

         ; Reject any fluxes based upon suspect PHOTO measurements,
         ; as indicated by the PHOTO flags.
         badbits2 = sdss_flagval('OBJECT2','SATUR_CENTER') $
          OR sdss_flagval('OBJECT2','INTERP_CENTER') $
          OR sdss_flagval('OBJECT2','PSF_FLUX_INTERP')
         qgoodphot = (tsobj.flags2 AND badbits2) EQ 0
         plugmap.calibflux = plugmap.calibflux * qgoodphot
         plugmap.calibflux_ivar = plugmap.calibflux_ivar * qgoodphot
      endif else begin
         splog, 'WARNING: No calibObj structure found for plate ', plateid
      endelse

      ;----------
      ; For any objects that do not have photometry from the calibObj
      ; structure, simply translate the flux from the plugmap MAG values
      ; (as long as those values are in the range 0 < MAG < +50).

      for ifilt=0, 4 do begin
         ibad = where(plugmap.calibflux[ifilt] EQ 0 $
          AND plugmap.mag[ifilt] GT 0 $
          AND plugmap.mag[ifilt] LT 50, nbad)
         if (nbad GT 0) then begin
            splog, 'Using plug-map fluxes for ', nbad, $
             ' values in filter ', ifilt
            plugmap[ibad].calibflux[ifilt] = $
             10.^((22.5 - plugmap[ibad].mag[ifilt]) / 2.5)
            plugmap[ibad].calibflux_ivar[ifilt] = 0
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
      for ifilt=0, 4 do $
       plugmap.mag[ifilt] = plugmap.mag[ifilt] - redden_med[ifilt]
   endif

   ; Optionally trim to selected spectrograph
   nfiber = n_elements(plugmap)
   if (keyword_set(spectrographid)) then begin
      indx = (spectrographid-1)*nfiber/2 + lindgen(nfiber/2)
      plugmap = plugmap[indx]
      fibermask = fibermask[indx]
   endif

   return, plugmap
end
;------------------------------------------------------------------------------
