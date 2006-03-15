;+
; NAME:
;   plotspecqa
;
; PURPOSE:
;   Wrapper for PLOTSPEC to plot standard stars and/or sky spectra.
;
; CALLING SEQUENCE:
;   plotspecqa, [ /standards, /skies, _EXTRA= ]
;
; INPUTS:
;   plate      - Plate number(s)
;
; OPTIONAL INPUTS:
;   standards  - Select standard stars for plotting
;   skies      - Select sky spectra for plotting
;   _EXTRA     - Kewords for PLOTSPEC
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   plotspec
;   readspec
;   splog
;
; REVISION HISTORY:
;   15-Mar-2006  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro plotspecqa, plate, fiberid, mjd=mjd, standards=standards, skies=skies, $
 _EXTRA=EXTRA

   if (keyword_set(standards) OR keyword_set(skies)) then begin
      for iplate=0L, n_elements(plate)-1L do begin
         if (keyword_set(mjd)) then thismjd = mjd[iplate] $
          else thismjd = 0
         readspec, plate[iplate], fiberid, mjd=thismjd, plugmap=plug
         if (NOT keyword_set(plug)) then begin
            splot, 'No plugmap structure found for the selected plate'
            return
         endif
         qstandards = keyword_set(standards) $
          AND (strmatch(plug.objtype,'SPECTROPHOTO_STD*') $
          OR strmatch(plug.objtype,'REDDEN_STD*'))
         qskies = keyword_set(skies) AND strmatch(plug.objtype,'SKY*')
         iplot = where(qstandards OR qskies, nplot)
         if (nplot EQ 0) then begin
            splog, 'No standard stars or skies selected'
            return
         endif
         if (keyword_set(standards)) then $
          splog, 'Selected ', fix(total(qstandards)), $
           ' standards on plate ', plate[iplate]
         if (keyword_set(skies)) then $
          splog, 'Selected ', fix(total(qskies)), $
           ' skies on plate ', plate[iplate]
         plotspec, plate[iplate], plug[iplot].fiberid, mjd=thismjd, $
          _EXTRA=EXTRA
      endfor
   endif else begin
      plotspec, plate, fiberid, mjd=mjd, _EXTRA=EXTRA
   endelse

   return
end
;------------------------------------------------------------------------------
