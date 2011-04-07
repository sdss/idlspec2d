;+
; NAME:
;   writespec
;
; PURPOSE:
;   Routine for writing Princeton-1D spectro outputs to an ASCII file
;
; CALLING SEQUENCE:
;   writespec, plate, fiberid, [ mjd=, filename= ]
;
; INPUTS:
;   plate      - Plate number(s)
;   fiber      - Fiber number(s), 1-indexed
;
; OPTIONAL INPUTS:
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mjd        - If not specified, then this is returned as an array of one
;                MJD per object.
;   filename   - Output file name; default to spec-pppp-mmmmm-fff.dat,
;                where pppp=plate number, mmmmm=MJD, fff=fiber ID.
;
; COMMENTS:
;   The data are read with READSPEC.  See the documentation for that
;   routine to see how to set environment variables that describe where
;   the data files are.
;
; EXAMPLES:
;   Write the spectrum of plate 401, fiber #100, to an ASCII file as follows: 
;     IDL> writespec, 401, 100
;   The default file name will be "spec-0401-51788-100.dat". This can be
;   changed by setting the FILENAME keyword: 
;     IDL> writespec, 401, 100, filename='junk.dat'
;
; BUGS:
;
; PROCEDURES CALLED:
;   readspec
;   sdss_flagname()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   25-Sep-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro writespec, plate, fiberid, mjd=mjd, filename=filename

   readspec, plate, fiberid, mjd=mjd, flux=objflux, flerr=objerr, $
    loglam=loglam, plug=plug
   wave = 10d^loglam

   primtarget = sdss_flagname('TARGET', plug.primtarget, /concat)
   sectarget = sdss_flagname('TTARGET', plug.sectarget, /concat)

   platestr = string(plate, format='(i4.4)')
   mjdstr = string(mjd, format='(i5.5)')
   fibstr = string(fiberid, format='(i4.4)')

   if (NOT keyword_set(filename)) then $
    filename = 'spec-' + platestr + '-' + mjdstr + '-' + fibstr + '.dat'

   openw, olun, filename, /get_lun

   printf, olun, '# Plate ' + platestr
   printf, olun, '# MJD ' + mjdstr
   printf, olun, '# Fiber ' + fibstr
   printf, olun, '# PRIMTARGET = ' + primtarget
   printf, olun, '# SECTARGET = ' + sectarget
   printf, olun, '#'
   printf, olun, '# Wavelength Flux       Error'
   printf, olun, '# [Ang]      [10^(-17) erg/cm/s/Ang]  [10^(-17) erg/cm/s/Ang]'
   for i=0, n_elements(objflux)-1 do begin
      printf, olun, wave[i], objflux[i], objerr[i], $
       format='(f10.3, e12.4, e12.4)'
   endfor

   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------
