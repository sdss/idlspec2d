;+
; NAME:
;   reduce_plate
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   reduce_plate, platenum, wave, template, [result, first= ]
;
; INPUTS:
;   platenum   - Plate number
;   wave       - wavelength of template  ([*,N] for N templates)
;   template   - template spectrum       ([*,N] for N templates)
;
; OPTIONAL INPUTS:
;   first      - only do the first "first" spectra on plate
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   result     - structure (see veldisp_struc) for info
;
; COMMENTS:
;   
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   readspec
;   veldisp
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   28-Jun-2000  Written - Doug Finkbeiner, Berkeley
;-
; D. Finkbeiner 28 Jun 2000  Master program to reduce 1 plate (1-d)


FUNCTION generate_filename, mjd, plate
  
  fname = string('spVelDisp-', mjd, '-', plate, $
                 format='(A,I5.5,A,I4.4)')
  return, fname
END 



FUNCTION get_specpath

  dum = findfile('/deep2/dfink', count=deepcount)
  IF (deepcount GT 0) THEN BEGIN   ; check if in Berkeley
      specpath = '/deepscr0/dfink/spectro/'
  ENDIF ELSE BEGIN 
      specpath = '/data/spectro/'
  ENDELSE
  print, 'SPEC path: ', specpath

  return, specpath
END



PRO reduce_plate, platenum, wave, template, result, first=first

  IF NOT keyword_set(platenum) THEN platenum = 306
  zap = 1  ; zap 5577

; Read in target list

  listfile = filepath('regress1d_all.dat', $
                      root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

  readcol, listfile[0], plt, mjd, fiber, zchic, class, targ, $
    format='(I,L,I,F,A,L)'

; define structure to hold results - store info from regress file
  result = veldisp_struc(n_elements(plt))
  result.plate = plt
  result.mjd   = mjd
  result.fiber = fiber
  result.zchic = zchic
  result.class = class
  result.primtarget = targ

; Consider only galaxies
  isgalaxy = (class EQ 'GALAXY') OR ((targ AND 96) NE 0)
  w = where((plt EQ platenum) AND isgalaxy, ct)
  print, ct, ' potential galaxies found'

  IF keyword_set(first) THEN BEGIN 
      w = w[0:first-1]
      ct = first
      print, 'List trimmed to ', ct
  ENDIF 

; Trim list
  IF ct EQ 1 THEN w = w[0]
  IF ct LT 1 THEN BEGIN 
      print, 'No objects meet criteria - skipping plate', platenum, '.'
      return
  END 
  result = result[w]

; Get spectro data path
  specpath = get_specpath()

; Read data
  readspec, platenum, result.fiber, flux=galflux, flerr=galsig, wave=galwave, plug=galplug, /silent, root_dir=specpath+'2d_3c/'

  result.run    = galplug.objid[0]
  result.rerun  = galplug.objid[1]
  result.camcol = galplug.objid[2]
  result.field  = galplug.objid[3]
  result.id     = galplug.objid[4]

; fix up template
  templatesig = template*0+1
  lambda_match, galwave[*, 0], wave, template
  lambda_match, galwave[*, 0], wave, templatesig
  lambda_match, galwave[*, 0], wave, wave

; Remove blank spectra from list
  bad = bytarr(ct)
  FOR i=0, ct-1 DO BEGIN 
      IF stdev(galflux[*, i]) EQ 0 THEN bad[i] = 1
  ENDFOR 
  w = where(bad EQ 0, nobj)
  print, 'Using ', nobj
  galflux = galflux[*, w]
  galwave = galwave[*, w]
  galsig  = galsig[*, w]
  galplug = galplug[w]
  result = result[w]

; Zap 5577
  IF keyword_set(zap) THEN BEGIN 
      linemask = (galwave LT 5586) AND (galwave GE 5573)
      galflux = djs_maskinterp(galflux, linemask, iaxis=0)
  ENDIF 

; Restrict wavelength range to "keep" (in Angstroms)
  keep = [3500, 6100]

; Call veldisp
  veldisp, galflux, galsig, galwave, template, templatesig, wave, result, $
    sigmast=0.05, maxsig=6, /nodif, keep=keep

  z     = result.z
  zchic = result.zchic

  badz = where(abs(z-zchic) GT 0.002, ct)
  IF ct EQ 0 THEN BEGIN 
      print, 'ALL z values agree with Chicago.'
  ENDIF ELSE BEGIN 
      print, 'Bad z fiber numbers: ', fiber[badz]
      print, 'Agree with Chicago', (nobj-ct), ' out of', nobj, ' = ', $
        (nobj-ct)/float(nobj)*100, '%'
  ENDELSE 

; Write veldisp FITS file
  fname = generate_filename(result[0].mjd, result[0].plate)
  mwrfits, result, specpath+'veldisp/brgtemplate/'+fname+'.fits', /create

; Write trimmed tsobj file
  suffix= '-'+string(result[0].mjd,format='(I5.5)')+ $
    '-'+string(platenum[0],format='(I4.4)')+'.fits'

  readspec, platenum, tsobj=tsobj
  help, tsobj
  tsobjtrim = tsobj[(result.fiber-1)]
  mwrfits, tsobjtrim,specpath+'veldisp/brgtemplate/tsObj-trim'+suffix, /create

  print
  print, 'Plate', platenum, ' finished.  ', systime()
  print

  return

END
