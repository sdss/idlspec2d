; D. Finkbeiner 28 Jun 2000  Master program to reduce 1 plate (1-d)


FUNCTION generate_filename, mjd, plate
  
  fname = string('spVelDisp-', mjd, '-', plate, '.fits', $
                 format='(A,I5.5,A,I4.4,A)')
  return, fname
END 

PRO reduce_plate, platenum, first=first

  IF NOT keyword_set(platenum) THEN platenum = 306
  zap = 1  ; zap 5577

; Read in target list
  listpath = '~/idlspec2d/etc/'
  readcol, listpath+'regress1d_all.dat', plt, mjd, fiber, zchic, class, targ, $
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

; Read template
  readcol, '~/brg.08.gk.01', wave, template, /silent

  dum = findfile('/deep2/dfink', count=deepcount)
  IF (deepcount GT 0) THEN BEGIN   ; check if in Berkeley
      path = '/deepscr0/dfink/spectro/'+string(platenum, format='(I4.4)')+'/'
      flist = findfile(path+'*.fits')
  ENDIF ELSE BEGIN 
      path = '/data/spectro/2d_3c/'+string(platenum, format='(I4.4)')+'/'
      flist = findfile(path+'2dnew/*.fits')
  ENDELSE 
  fname = flist[0]

; Get log lambda
  if (NOT keyword_set(hdr)) then hdr = headfits(fname)
  naxis1 = sxpar(hdr, 'NAXIS1')
  coeff0 = sxpar(hdr, 'COEFF0')
  coeff1 = sxpar(hdr, 'COEFF1')
  loglam = coeff0 + coeff1 * findgen(naxis1)

; Read data
  galflux = (mrdfits(fname, 0))[*, result.fiber]
  galsig  = (mrdfits(fname, 1))[*, result.fiber]
  galwave = (10.^loglam)#fltarr(n_elements(result))

;  readspec, result.plate, 0, flux=galflux, wave=galwave, $
;    flerr=galsig, /silent

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
;  galplug = galplug[*, w]
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

  fname = generate_filename(result[0].mjd, result[0].plate)
  mwrfits, result, '/data/spectro/veldisp/brgtemplate/'+fname

  print
  print, 'Plate', platenum, ' finished.  ', systime()
  print
  return

END
