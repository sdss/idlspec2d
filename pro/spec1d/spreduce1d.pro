;+
; NAME:
;   spreduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   spreduce1d, [ platefile, fiberid= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platefile  - Plate file(s) from spectro-2D; default to all files
;                matching 'spPlate*.fits'
;   fiberid    - If specified, then only reduce these fiber numbers;
;                this must be a vector with unique values between 1 and 640.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Names of output files are derived from PLATEFILE.
;   For example, if PLATEFILE='spPlate-0306-51690.fits', then
;     ZALLFILE = 'spZall-0306-51690.fits'
;     ZBESTFILE = 'spZbest-0306-51690.fits'
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   mrdfits()
;   mwrfits
;   splog
;   sxaddpar
;   sxdelpar
;   sxpar()
;   zfind()
;   zrefind()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro spreduce1d, platefile, fiberid=fiberid

   if (NOT keyword_set(platefile)) then platefile = findfile('spPlate*.fits')

   ;----------
   ; If multiple plate files exist, then call this script recursively
   ; for each such plate file.

   if (n_elements(platefile) GT 1) then begin
      for i=0, n_elements(platefile)-1 do begin
         spreduce1d, platefile[i]
      endfor
      return
   endif else begin
      platefile = platefile[0]
   endelse

   ;----------
   ; Determine names of output files

   platemjd = strmid(platefile, 8, 10)

   zallfile = 'spZall-' + platemjd + '.fits'
   zbestfile = 'spZbest-' + platemjd + '.fits'
   if (NOT keyword_set(logfile)) then $
    logfile = 'spDiag1d-' + platemjd + '.log'

   splog, filename=logfile
   splog, 'Log file ', logfile, ' opened ', systime()
   splog, 'IDL version: ', string(!version,format='(99(a," "))')
   stime0 = systime(1)

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   npixobj = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
   objivar = mrdfits(platefile,1)
   andmask = mrdfits(platefile,2)
;   ormask = mrdfits(platefile,3)
;   dispmap = mrdfits(platefile,4)
   plugmap = mrdfits(platefile,5)

   objivar = skymask(objivar, andmask)
andmask = 0 ; Free memory

   ;----------
   ; Trim to specified fibers if FIBERID is set

   if (keyword_set(fiberid)) then begin
      objflux = objflux[*,fiberid-1]
      objivar = objivar[*,fiberid-1]
      plugmap = plugmap[fiberid-1]
      nobj = n_elements(fiberid)
   endif else begin
      fiberid = lindgen(nobj) + 1
   endelse

   ;----------
   ; Look for where the S/N is unreasonably large

   for iobj=0L, nobj-1 do begin
      junk = where(objflux[*,iobj] * sqrt(objivar[*,iobj]) GT 200., ct)
      if (ct GT 0) then $
       splog, 'WARNING: Fiber #', iobj+1, ' has ', ct, ' pixels with S/N > 200'
   endfor

   ;----------
   ; Find GALAXY redshifts

   npoly = 3
   zmin = -0.01 ; -3000 km/sec
   zmax = 0.60 ; Max z for a rest-frame template to 2300 Ang to cover 3700 Ang
   pspace = 2
   nfind = 5

   eigenfile = 'spEigenGal*.fits'

   splog, 'Compute GALAXY redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_gal = zfind(objflux, objivar, hdr=hdr, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=pspace*3)
   splog, 'CPU time to compute GALAXY redshifts = ', systime(1)-t0

   splog, 'Locally re-fitting GALAXY redshifts'
   t0 = systime(1)
   res_gal = zrefind(objflux, objivar, hdr=hdr, $
    pwidth=5, pspace=1, width=5, zold=res_gal)
   splog, 'CPU time to re-fit GALAXY redshifts = ', systime(1)-t0

   res_gal.class = 'GALAXY'
   res_gal.subclass = ' '

   res_all = res_gal ; Append results

   ;----------
   ; Find QSO redshifts

   npoly = 3
   zmin = 0.0033 ; +1000 km/sec
   zmax = 6.00 ; Max range to use for now, with the template starting at
               ; 525 Ang (rest), which corresponds to 3700 Ang at this z.
   pspace = 4
   nfind = 5

   eigenfile = 'spEigenQSO*.fits'

   splog, 'Compute QSO redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_qso = zfind(objflux, objivar, hdr=hdr, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=pspace*5)
   splog, 'CPU time to compute QSO redshifts = ', systime(1)-t0

   splog, 'Locally re-fitting QSO redshifts'
   t0 = systime(1)
   res_qso = zrefind(objflux, objivar, hdr=hdr, $
    pwidth=11, pspace=1, width=11, zold=res_qso)
   splog, 'CPU time to re-fit QSO redshifts = ', systime(1)-t0

   res_qso.class = 'QSO'
   res_qso.subclass = ' '

   res_all = [res_all, res_qso] ; Append results

   ;----------
   ; Find STAR redshifts

   npoly = 4 ; With only 1 eigen-template, fit more polynomial terms for stars.
   zmin = -0.002 ; -600 km/sec
   zmax = 0.002 ; +600 km/sec
   pspace = 1
   nfind = 1

   eigenfile = 'spEigenStar*.fits'

   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   shdr = headfits( djs_filepath(eigenfile, root_dir=eigendir) )
   nstar = sxpar(shdr, 'NAXIS2') > 1

   for istar=0, nstar-1 do begin
      subclass = sxpar(shdr, 'NAME'+strtrim(string(istar),2))

      splog, 'Compute STAR (' + subclass + ') redshifts:', $
       ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
      t0 = systime(1)
      res_star = zfind(objflux, objivar, hdr=hdr, $
       eigenfile=eigenfile, columns=istar, npoly=npoly, $
       zmin=zmin, zmax=zmax, pspace=1, $
       nfind=nfind, width=pspace*3)
      splog, 'CPU time to compute STAR redshifts = ', systime(1)-t0

      res_star.class = 'STAR'
      res_star.subclass = subclass

      res_all = [res_all, res_star] ; Append results
   endfor

   ;----------
   nper = (size(res_all,/dimens))[0]

   ;----------
   ; Sort results for each object by ascending order in chi^2/DOF,
   ; but putting any results with zero degrees-of-freedom at the end.

   minvdiff = 1000.0 ; km/s
   cspeed = 2.9979e5

   for iobj=0, nobj-1 do begin
      res1 = res_all[*,iobj]

      rchi2 = res1.rchi2

      isort = sort(rchi2 + (res1.dof EQ 0)*max(rchi2))
      for ii=0, nper-1 do begin
         res_all[ii,iobj] = res1[isort[ii]]
      endfor

      ; Find the difference in reduced chi^2 between each result and the next
      res1 = res_all[*,iobj]
      rchi2 = res1.rchi2
      for ii=0, nper-2 do begin
         inext = (where( $
          abs(res1[ii+1:nper-1].z - res1[ii].z) GT minvdiff/cspeed $
          AND res1[ii+1:nper-1].dof GT 0))[0]
         if (inext NE -1) then $
          res_all[ii,iobj].rchi2diff = rchi2[ii+1+inext] - rchi2[ii]
      endfor
   endfor

   ;----------
   ; Insist that all SKY fibers are identified as SKY

   for iobj=0, nobj-1 do begin
      if (strtrim(plugmap[iobj].objtype,2) EQ 'SKY') then begin
         if (nper GT 1) then $
          res_all[1:nper-1,iobj] = res_all[0:nper-2,iobj]
         rcopy = { class: 'SKY', $
                   subclass: ' ', $
                   tfile: ' ' }
         res1 = res_all[0,iobj]
         struct_assign, rcopy, res1
         res_all[0,iobj] = res1
      endif
   endfor

   ;----------
   ; Re-classify some objects as UNKNOWN.  Do this for any non-SKY objects
   ; where:
   ; (1) The reduced chi^2 difference with the next best identification
   ;     is very small.
   ; (2) The classification is a QSO where the first eigen-component
   ;     is negative.
   ; (2) The classification is a STAR where the first eigen-component
   ;     is negative.

   minrchi2diff = 0.01

   qlowconf = bytarr(nobj)
   qlowconf = qlowconf OR (res_all[0,*].rchi2diff LT minchi2diff)
   qlowconf = qlowconf OR (strtrim(res_all[0,*].class) EQ 'STAR' $
    AND res_all[0,*].theta[0] LT 0)
   qlowconf = qlowconf OR (strtrim(res_all[0,*].class) EQ 'QSO' $
    AND res_all[0,*].theta[0] LT 0)

   for iobj=0, nobj-1 do begin
      if (qlowconf[0] EQ 0 AND res_all[0,iobj].class NE 'SKY') then begin
         if (nper GT 1) then $
          res_all[1:nper-1,iobj] = res_all[0:nper-2,iobj]
         rcopy = { class: 'UNKNOWN', $
                   subclass: ' ', $
                   tfile: ' ' }
         res1 = res_all[0,iobj]
         struct_assign, rcopy, res1
         res_all[0,iobj] = res1
      endif
   endfor

   ;----------
   ; Add other fields to the output structure

   res1 = { plate:    long(sxpar(hdr, 'PLATEID')), $
            tile:     long(sxpar(hdr, 'TILEID')), $
            mjd:      long(sxpar(hdr, 'MJD')), $
            fiberid:  0L        , $
            objid:    lindgen(5), $
            objtype:  ' '       , $
            plug_ra:  0.0d      , $
            plug_dec: 0.0d      }
   res_prepend = make_array(value=res1, dimension=size(res_all,/dimens))
   res_all = struct_addtags(res_prepend, res_all)

   for iobj=0, nobj-1 do begin
      res_all[*,iobj].fiberid = fiberid[iobj]
      res_all[*,iobj].objid = plugmap[iobj].objid
      res_all[*,iobj].objtype = plugmap[iobj].objtype
      res_all[*,iobj].plug_ra = plugmap[iobj].ra
      res_all[*,iobj].plug_dec = plugmap[iobj].dec
   endfor

   res1 = { wavemin:   0.0, $
            wavemax:   0.0, $
            wcoverage: 0.0, $
            spec1_g:   sxpar(hdr, 'SPEC1_G'), $
            spec1_r:   sxpar(hdr, 'SPEC1_R'), $
            spec1_i:   sxpar(hdr, 'SPEC1_I'), $
            spec2_g:   sxpar(hdr, 'SPEC2_G'), $
            spec2_r:   sxpar(hdr, 'SPEC2_R'), $
            spec2_i:   sxpar(hdr, 'SPEC2_I') }
   res_append = make_array(value=res1, dimension=size(res_all,/dimens))
   res_all = struct_addtags(res_all, res_append)

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')
   for iobj=0, nobj-1 do begin
      igood = where(objivar[*,iobj] NE 0, ngood)
      res_all[*,iobj].wavemin = $
       10^(objloglam0 + (igood[0]>0)*objdloglam) * (ngood NE 0)
      res_all[*,iobj].wavemax = $
       10^(objloglam0 + (igood[(ngood-1)>0])*objdloglam) * (ngood NE 0)
      res_all[*,iobj].wcoverage = ngood * objdloglam
   endfor

   ;----------
   ; Write the output files

   sxaddpar, hdr, 'NAXIS', 0
   sxdelpar, hdr, 'NAXIS1'
   sxdelpar, hdr, 'NAXIS2'
   sxaddpar, hdr, 'EXTEND', 'T', after='NAXIS'
   sxaddpar, hdr, 'VERS1D', idlspec2d_version(), $
    'Version of idlspec2d for 1D reduction', after='VERSCOMB'

   mwrfits, 0, zallfile, hdr, /create ; Retain the original header in first HDU
   mwrfits, res_all, zallfile

   mwrfits, 0, zbestfile, hdr, /create ; Retain the original header in first HDU
   mwrfits, (res_all[0,*])[*], zbestfile

   ;----------
   ; Close log file

   splog, 'Total time for SPREDUCE1D = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE1D at ', systime()
   splog, /close

   return
end
;------------------------------------------------------------------------------
