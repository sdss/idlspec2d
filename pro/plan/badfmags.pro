;+
; NAME:
;   badfmags
;
; PURPOSE:
;   Find any F stars where the magnitudes in the calibObj files are grossly
;   discrepent from the plug-map files.
;
; CALLING SEQUENCE:
;   badfmags, [ maxdiff= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   maxdiff - Maximum r-band magnitude difference; default to 0.4
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine looks for calibration F stars that have very
;   discrepent magnitudes in the plug-map file vs. the calibObj files.
;   Only the r-band magnitudes are compared, and they must differ by
;   at least MAXDIFF in order to trigger a warning.
;
;   Before comparing the magnitudes, first offset all the F star
;   magnitudes by their median difference.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   platelist
;   readspec
;   splog
;
; REVISION HISTORY:
;   10-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro badfmags, maxdiff=maxdiff

   if (NOT keyword_set(maxdiff)) then maxdiff = 0.4

   splog, file='badfmag.log'

   platelist, plist=plist
   plist = plist[where(strmatch(plist.statuscombine,'Done*'))]

   for i=0L, n_elements(plist)-1 do begin
      readspec, plist[i].plate, mjd=plist[i].mjd, plug=plug, tsobj=calibobj, $
       /silent
      iphoto = where((plug.sectarget AND 34) NE 0, nphoto)
      if (keyword_set(calibobj) AND nphoto GT 1) then begin
         oldmag = plug[iphoto].mag
         newmag = 22.5 - 2.5*alog10(calibobj[iphoto].psfflux > 0.1)
         mdiff = newmag - oldmag
         for k=0, 4 do mdiff[k,*] = mdiff[k,*] - median(mdiff[k,*])
         for j=0, nphoto-1 do begin
            if (abs(mdiff[2,j]) GT maxdiff) then $
             splog, plist[i].plate, plist[i].mjd, iphoto[j]+1, $
              oldmag[2,j], newmag[2,j], mdiff[2,j]
         endfor
      endif
   endfor

   splog, /close

   return
end
;------------------------------------------------------------------------------
