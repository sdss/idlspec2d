;+
; NAME:
;   readspec_footprint
;
; PURPOSE:
;   Routine for reading 1D spectro outputs on a given plate footprint
;
; CALLING SEQUENCE:
;   readspec, plate, [mjd=, topdir=, /best, /silent, $
;    zans=, plugmap=, tsobj= ]
;
; INPUTS:
;   plate      - Plate number
;
; OPTIONAL INPUTS:
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
;   best       - If set, then select the best observation of each object,
;                using the same algorithm to choose the primary observation
;                in the PLATEMERGE routine (for spAll files).
;   silent     - If set, then call MRDFITS with /SILENT.
;   topdir     - Top-level directory for data; default to the environment
;                variable $BOSS_SPECTRO_REDUX.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mjd        - If not specified, then this returns the MJD of PLATE
;   zans       - Redshift output structure [NFIBER]
;   plugmap    - Plug-map entries [NFIBER]
;   tsobj      - tsObj-structure output [NFIBER]
;
; COMMENTS:
;   Read all of the SDSS spectroscopic observations within 1.49 degrees
;   of the center of the specified plate.  If /BEST is set, then trim
;   to the best observation of each object.
;
; EXAMPLES:
;   Read all the SDSS spectroscopic observations on the SDSS-3 test
;   plate 2634, trimming to the best observations of each object:
;     IDL> readspec_footprint, 2634, zans=zans, /best
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   platelist
;   readspec
;   spheregroup()
;   sxpar()
;
; REVISION HISTORY:
;   24-Apr-2007  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro readspec_footprint, platenum, mjd=mjd, best=best, silent=silent, $
 zans=zans, plugmap=plug, tsobj=tsobj

   if (n_elements(platenum) EQ 0) then $
    message, 'PLATE is not specified or is not a scalar'

   ; Find all plates that may overlap this first plate
   readspec, platenum, mjd=mjd, objhdr=hdr, silent=silent
   if (NOT keyword_set(hdr)) then $
    message, 'Plate does not exist ' + string(platenum)
   ra = sxpar(hdr, 'RA')
   dec = sxpar(hdr, 'DEC')
   platelist, plist=plist
   plist = plist[where(strmatch(plist.status1d,'Done*'))]
;   plist = plist[where(strmatch(plist.status1d,'Done*') $
;    AND (strmatch(plist.platequality,'good*') $
;    OR strmatch(plist.platequality,'marginal*')))]
   adist = djs_diff_angle(plist.ra, plist.dec, ra, dec)
   plist = plist[where(adist LE 3.)]

   ; Read all these plates
   readspec, plist.plate, mjd=plist.mjd, zans=zans, silent=silent
   if (arg_present(plug)) then $
    readspec, plist.plate, mjd=plist.mjd, plug=plug, silent=silent
   if (arg_present(tsobj)) then $
    readspec, plist.plate, mjd=plist.mjd, tsobj=tsobj, silent=silent

   ; Store the plate quality for each object
   pquality = strarr(n_elements(zans))
   for iplate=0L, n_elements(plist)-1L do begin
      indx = where(zans.plate EQ plist[iplate].plate $
       AND zans.mjd EQ plist[iplate].mjd)
      pquality[indx] = plist[iplate].platequality
   endfor

   ; Trim to only those within the footprint of the requested plate
   adist = djs_diff_angle(zans.plug_ra, zans.plug_dec, ra, dec)
   indx = where(adist LT 1.49)
   zans = zans[indx]
   if (arg_present(plug)) then plug = plug[indx]
   if (arg_present(tsobj)) then tsobj = tsobj[indx]
   pquality = pquality[indx]

   ; Trim to the best observation of each object, using the same
   ; scoring algorithm as in PLATEMERGE
   if (keyword_set(best)) then begin
      specprimary = bytarr(n_elements(zans))
      score = 4 * (zans.sn_median GT 0) $
       + 2 * (strmatch(pquality,'good*') EQ 1) $
       + 1 * (zans.zwarning EQ 0) $
       + (zans.sn_median>0) / max(zans.sn_median+1.)
      dtheta = 2.0 / 3600.
      ingroup = spheregroup(zans.plug_ra, zans.plug_dec, dtheta, $
       multgroup=multgroup, firstgroup=firstgroup, nextgroup=nextgroup)
      for j=0L, n_elements(firstgroup)-1L do begin
         if (firstgroup[j] NE -1) then begin 
            if (multgroup[j] EQ 1) then begin
               specprimary[firstgroup[j]] = 1
            endif else begin
               indx = lonarr(multgroup[j])
               indx[0] = firstgroup[j]
               for k=0L, multgroup[j]-2L do indx[k+1] = nextgroup[indx[k]]
               foo = max(score[indx], ibest)
               specprimary[indx[ibest]] = 1
            endelse
         endif
      endfor
      indx = where(specprimary)
      zans = zans[indx]
      if (arg_present(plug)) then plug = plug[indx]
      if (arg_present(tsobj)) then tsobj = tsobj[indx]
   endif

   return
end
;------------------------------------------------------------------------------
