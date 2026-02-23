;+
; NAME:
;   apo_plotbias
;
; PURPOSE:
;   Wrapper for sos_plotbias (to keep apo's procedures valid)
;
; CALLING SEQUENCE:
;   apo_plotbias, expnum, [ plotfile= ]
;
; INPUTS:
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   plotfile   - Plot file; if set, then send plot to this PostScript file
;                rather than to the default (X) display.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The histogram of bias values is plotted for all (4) camera files that
;   match the given exposure number.
;
;   A fiducial line is drawn as a thick blue line.  This line approximates
;   what we expect to see for each camera.
;
;   If $BOSS_SPECTRO_DATA is not set, then it is assumed to be
;     /data/spectro
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;    sos_plotbias
; REVISION HISTORY:
;   13-Jan-2025 Sean Morrison (UIUC)
;-
;------------------------------------------------------------------------------


pro apo_plotbias, expnum, plotfile=plotfile

   if (n_params() LT 1) then begin
      print, 'Syntax - apoplotbias, expnum, [plotfile= ]'
      return
   endif

    sos_plotbias, expnum, plotfile=plotfile
end
