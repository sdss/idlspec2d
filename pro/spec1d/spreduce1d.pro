;+
; NAME:
;   spreduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   spreduce1d, platefile, [outfile]
;
; INPUTS:
;   platefile  - Plate file(s) from spectro-2D; default to all files
;                matching 'spPlate*.fits'
;
; OPTIONAL INPUTS:
;   outfile    - Name of output file; default to a name derived from PLATEFILE.
;                For example, if PLATEFILE='spPlate-0306-51690.fits', then
;                OUTFILE='spZ-0306-51690.fits'.
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
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   mrdfits()
;   mwrfits
;   splog
;   sxpar()
;   writefits
;   zfind()
;   zrefind()
;
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro spreduce1d, platefile, outfile

   if (NOT keyword_set(platefile)) then platefile = findfile('spPlate*.fits')

   ;----------
   ; If multiple plate files exist, then call this script recursively
   ; for each such plate file.

   if (n_elements(platefile) GT 1) then begin
      for i=0, n_elements(platefile)-1 do begin
         if (keyword_set(outfile)) then $
          spreduce1d, platefile[i], outfile[i] $
         else $
          spreduce1d, platefile[i]
      endfor
      return
   endif else begin
      platefile = platefile[0]
   endelse

   ;----------
   ; Determine name of output file

   if (NOT keyword_set(outfile)) then $
    outfile = 'spZ' + strmid(platefile,strpos(platefile,'spPlate')+7)

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

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   skymask = (andmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for iobj=0L, nobj-1 do $
    skymask[*,iobj] = smooth(float(skymask[*,iobj]),5) GT 0
   ibad = where(skymask)
andmask = 0 ; Free memory
   if (ibad[0] NE -1) then objivar[ibad] = 0

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
   zmax = 0.5 ; Limit for now until template goes further into the UV ???
   pspace = 2
   nfind = 5

   eigenfile = 'spEigenGals.fits'

   splog, 'Compute GALAXY redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_gal = zfind(objflux, objivar, hdr=hdr, fiberid=plugmap.fiberid, $
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
   zmax = 3.75 ; Max range to use for now, with the template starting at
               ; 800 Ang (rest), which corresponds to 3800 Ang at this z. ???
   pspace = 4
   nfind = 5

   eigenfile = 'spEigenQSO.fits'

   splog, 'Compute QSO redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_qso = zfind(objflux, objivar, hdr=hdr, fiberid=plugmap.fiberid, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=pspace*3)
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

   eigenfile = 'spEigenStar.fits'

   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   shdr = headfits( djs_filepath(eigenfile, root_dir=eigendir) )
   nstar = sxpar(shdr, 'NAXIS2') > 1

   for istar=0, nstar-1 do begin
      subclass = sxpar(shdr, 'NAME'+strtrim(string(istar),2))

      splog, 'Compute STAR (' + subclass + ') redshifts:', $
       ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
      t0 = systime(1)
      res_star = zfind(objflux, objivar, hdr=hdr, fiberid=plugmap.fiberid, $
       eigenfile=eigenfile, columns=istar, npoly=npoly, $
       zmin=zmin, zmax=zmax, pspace=1, $
       nfind=nfind, width=pspace*3)
      splog, 'CPU time to compute STAR redshifts = ', systime(1)-t0

      res_star.class = 'STAR'
      res_star.subclass = subclass

      res_all = [res_all, res_star] ; Append results
   endfor

   ;----------
   ; Sort results for each object by ascending order in chi^2/DOF,
   ; but putting any results with zero degrees-of-freedom at the end.

   for iobj=0, nobj-1 do begin
      res1 = res_all[*,iobj]
;      sortval = res1.chi2 / (res1.dof > 1)
      sortval = (res1.chi2 + (res1.dof EQ 0) * max(res1.chi2)) $
       / (res1.dof + (res1.dof EQ 0))
      isort = sort(sortval)
      for ii=0, n_elements(res1)-1 do begin
         res_all[ii,iobj] = res1[isort[ii]]
      endfor
   endfor

   ;----------
   ; Insist that all SKY fibers are identified as SKY

   nper = (size(res_all,/dimens))[0]
   for iobj=0, nobj-1 do begin
      if (strtrim(plugmap[iobj].objtype,2) EQ 'SKY') then begin
         if (nper GT 1) then $
          res_all[1:nper-1,iobj] = res_all[0:nper-2,iobj]
         rcopy = { plate: res_all[0,iobj].plate, $
                   mjd: res_all[0,iobj].mjd, $
                   fiberid: res_all[0,iobj].fiberid, $
                   class: 'SKY', $
                   subclass: ' ', $
                   tfile: ' ' }
         res1 = res_all[0,iobj]
         struct_assign, rcopy, res1
         res_all[0,iobj] = res1
      endif
   endfor

   splog, 'Total time for SPREDUCE1D = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE1D at ', systime()

   ;----------
   ; Write the output file

   sxaddpar, hdr, 'NAXIS', 0
   writefits, outfile, 0, hdr ; Retain the original header in the first HDU
   mwrfits, res_all, outfile

   return
end
;------------------------------------------------------------------------------
