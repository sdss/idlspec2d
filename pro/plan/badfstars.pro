;+
; NAME:
;   badfstars
;
; PURPOSE:
;   Find bad F stars used for spectro-photometry in the Spectro-2D reductions
;
; CALLING SEQUENCE:
;   badfstars, [ /doplot ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   doplot     - If set, then use PLOTSPEC to plot all the bad spectra
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This procedure looks at all the objects labelled as either
;   OBJTYPE='SPECTROPHOTO_STD' or 'REDDEN_STD' in the spAll.fits file.
;   These are considered to be bad F stars if the following conditions
;   are satisfied:
;     WCOVERAGE > 0.10
;     |cz| > 500 km/sec OR CLASS != 'STAR'
;   These objects are then listed in the file 'badfstars.log', and then
;   plotted if /DOPLOT is set.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   hogg_mrdfits()
;   plotspec
;
; REVISION HISTORY:
;   06-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro badfstars, doplot=doplot

   spfile = filepath('spAll.fits', root_dir=getenv('SPECTRO_DATA'))
   columns = ['PLATE','FIBERID','MJD', $
    'PRIMTARGET','SECTARGET','OBJTYPE', $
    'RA','DEC','CLASS','Z','ZWARNING','WCOVERAGE']

   spall = hogg_mrdfits(spfile, 1, columns=columns, nrowchunk=10000L)

   qfstar = strmatch(spall.objtype,'SPECTROPHOTO_STD*') $
    OR strmatch(spall.objtype,'REDDEN_STD*')
   ibad = where(qfstar AND spall.wcoverage GT 0.10 $
    AND (abs(spall.z) GT 500./3e5 OR strmatch(spall.class,'STAR*') EQ 0), nbad)

   splog, file='badfstars.log'
   for j=0L, nbad-1 do $
    splog, spall[ibad[j]].plate, spall[ibad[j]].mjd, spall[ibad[j]].fiberid, $
     format='(I5,I6,I4)', /noname
   splog, /close

   if (keyword_set(doplot)) then $
    plotspec,spall[ibad].plate,spall[ibad].fiberid,mjd=spall[ibad].mjd

   return
end
;------------------------------------------------------------------------------
