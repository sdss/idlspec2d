;+
; NAME:
;   spreduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   spreduce1d, platefile
;
; INPUTS:
;   platefile  - Plate file from spectro-2D
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Presently, just writes a file 'z-'+PLATEFILE.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   mwrfits
;   splog
;   sxpar()
;   zfind()
;
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro spreduce1d, platefile

platefile = 'spPlate-0306-51690.fits'
outfile = 'z-' + platefile


djs_readcol, '/home/schlegel/idlspec2d/etc/regress1d_all.dat', $
 chicplate, junk, chicfiberid, chicz, chicclass, format='(L,L,L,F,A)'
ii=where(chicplate EQ 306)
chicz=chicz[ii]
chicclass=chicclass[ii]

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

; Trim to QSO's only...
;jj = where(chicclass EQ 'QSO')
;jj = where(chicclass EQ 'GALAXY')
;jj = jj[0:4] ; Trim to first few objects
jj = lindgen(6)
chicz = chicz[jj]
chicclass = chicclass[jj]
objflux = objflux[*,jj]
objivar = objivar[*,jj]
plugmap = plugmap[jj]
nobj=n_elements(jj)

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
   zmin = -0.01
   zmax = 0.6
zmax = 0.25 ; ???
   pspace = 2
   nfind = 5

   eigenfile = filepath('spEigenGals.fits', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

   splog, 'Compute GALAXY redshifts: ZMIN=', zmin, ' ZMAX=', zmax, $
    ' PSPACE=', pspace
   t0 = systime(1)
   res_gal = zfind(objflux, objivar, hdr=hdr, fiberid=plugmap.fiberid, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=pspace*3)
   splog, 'CPU time to compute GALAXY redshifts = ', systime(1)-t0

   res_gal.class = 'GALAXY'

   ;----------
   ; Find QSO redshifts

   npoly = 3
   zmin = -0.01
   zmax = 5.0
zmax = 2.2 ; ???
   pspace = 4
   nfind = 5

   eigenfile = filepath('spEigenQSO.fits', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

   splog, 'Compute QSO redshifts: ZMIN=', zmin, ' ZMAX=', zmax, $
    ' PSPACE=', pspace
   t0 = systime(1)
   res_qso = zfind(objflux, objivar, hdr=hdr, fiberid=plugmap.fiberid, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=pspace*3)
   splog, 'CPU time to compute QSO redshifts = ', systime(1)-t0

   res_qso.class = 'QSO'

   ;----------
   ; Append and sort results for each object by chi^2/DOF

   res_all = [res_gal, res_qso]

   for iobj=0, nobj-1 do begin
      res1 = res_all[*,iobj]
      sortval = res1.chi2 / (res1.dof > 1)
      isort = sort(sortval)
      for ii=0, n_elements(res1)-1 do begin
         res_all[ii,iobj] = res1[isort[ii]]
      endfor
   endfor

print,[transpose(chicz),res_all[0,*].z]

;   vres = veldisp(objflux, objivar, starflux, starmask, $
;    zoffset=zoffset)

   ;----------
   ; Write the output file

stop

set_plot,'ps'
device,file='qsoz-306a.ps'
djs_plot, chicz, davez,ps=4,xr=[0,3],yr=[0,3], charsize=2, $
 xtitle='Chicago-z',ytitle='z from \chi^2-min',title='QSOs on plate 306'
oplot,[0,3],[0,3]
device,/close
set_plot,'x'

set_plot,'ps'
device,file='qsoz-306b.ps'
djs_plot, chicz, ((davez-chicz)*3e5 > (-5800))<5800, $
 ps=4,yr=[-1,1]*6000,xr=[0,3], charsize=2, $
 xtitle='Chicago-z',ytitle='\Delta z from \chi^2-min',title='QSOs on plate 306'
oplot,[0,3],[0,0]
device,/close
set_plot,'x'

djs_plot, chicz, ((davez-chicz)*3e5 > (-580))<580, $
 ps=4,yr=[-1,1]*600,xr=[0,1], charsize=2
j=where( abs(davez-chicz)*3e5 LT 400 )

j=where( abs(davez-chicz)*3e5 LT 4000 )
print,mean((davez-chicz)[j]*3e5),djsig((davez-chicz)[j]*3e5)
print,n_elements(j)/float(nobj)

   mwrfits, result, outfile, /create

   return
end
;------------------------------------------------------------------------------
