;+
; NAME:
;   spreduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   spreduce1d, [ platefile, fiberid=, /doplot, /debug ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platefile  - Plate file(s) from spectro-2D; default to all files
;                matching 'spPlate*.fits'
;   fiberid    - If specified, then only reduce these fiber numbers;
;                this must be a vector with unique values between 1 and
;                the number of rows in the plate file (typically 640).
;   doplot     - If set, then generate plots.  Send plots to a PostScript
;                file spDiagDebug1d-$PLATE-$MJD.ps unless /DEBUG is set.
;   debug      - If set, then send plots to the X display and wait for
;                a keystroke after each plot; setting /DEBUG forces /DOPLOT.
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
;   cpbackup
;   dfpsclose
;   dfpsplot
;   elodie_best()
;   filter_thru()
;   mrdfits()
;   mwrfits
;   qaplot_fcalibvec
;   splog
;   skymask()
;   speclinefit
;   struct_addtags()
;   sxaddpar
;   sxdelpar
;   sxpar()
;   synthspec()
;   vdispfit
;   zfind()
;   zrefind()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro spreduce1d, platefile, fiberid=fiberid, doplot=doplot, debug=debug

   if (NOT keyword_set(platefile)) then begin
      platefile = findfile('spPlate*.fits*', count=nplate)
   endif else begin
      if (size(platefile,/tname) NE 'STRING') then $
       message, 'PLATEFILE must be a file name'
      if (keyword_set(platefile)) then nplate = n_elements(platefile) $
       else nplate = 0
   endelse
   if (keyword_set(debug)) then doplot = 1

   ;----------
   ; If multiple plate files exist, then call this script recursively
   ; for each such plate file.

   if (nplate EQ 0) then begin
      splog, 'No plate files specified or found'
      return
   endif else if (nplate EQ 1) then begin
      platefile = platefile[0]
   endif else begin
      for i=0, nplate-1 do begin
         spreduce1d, platefile[i], fiberid=fiberid, doplot=doplot, debug=debug
      endfor
      return
   endelse

   ;----------
   ; Determine names of output files

   platemjd = strmid(platefile, 8, 10)

   zallfile = 'spZall-' + platemjd + '.fits'
   zbestfile = 'spZbest-' + platemjd + '.fits'
   zlinefile = 'spZline-' + platemjd + '.fits'
   if (NOT keyword_set(logfile)) then $
    logfile = 'spDiag1d-' + platemjd + '.log'
   plotfile = 'spDiag1d-' + platemjd + '.ps'

   if (keyword_set(doplot) AND NOT keyword_set(debug)) then begin
      debugfile = 'spDiagDebug1d-' + platemjd + '.ps'
      cpbackup, debugfile
      dfpsplot, debugfile, /color
   endif

   stime0 = systime(1)

   if (keyword_set(logfile)) then begin
      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + systime()
   endif
   if (keyword_set(plotfile)) then $
    splog, 'Plot file ' + plotfile
   if (keyword_set(debugfile)) then $
    splog, 'Debug plot file ' + debugfile
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   if (NOT keyword_set(hdr)) then $
    message, 'Plate file not valid: ' + platefile
   npixobj = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
   objivar = mrdfits(platefile,1)
   andmask = mrdfits(platefile,2)
   ormask = mrdfits(platefile,3)
;   dispmap = mrdfits(platefile,4)
   plugmap = mrdfits(platefile,5)

   anyandmask = transpose(andmask[0,*])
   anyormask = transpose(ormask[0,*])
   for ipix=1, npixobj-1 do $
    anyandmask = anyandmask OR transpose(andmask[ipix,*])
   for ipix=1, npixobj-1 do $
    anyormask = anyormask OR transpose(ormask[ipix,*])

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')

   ;----------
   ; Trim to specified fibers if FIBERID is set

   if (keyword_set(fiberid)) then begin
      if (min(fiberid) LE 0 OR max(fiberid) GT nobj) then $
       message, 'Invalid value for FIBERID: must be between 0 and '+string(nobj)
      objflux = objflux[*,fiberid-1]
      objivar = objivar[*,fiberid-1]
      anyandmask = anyandmask[fiberid-1]
      anyormask = anyormask[fiberid-1]
      plugmap = plugmap[fiberid-1]
      nobj = n_elements(fiberid)
   endif else begin
      fiberid = lindgen(nobj) + 1
   endelse
   splog, 'Number of fibers = ', nobj

   ;----------
   ; Look for where the S/N is unreasonably large
   ; or where flux is unphysically negative.

   for iobj=0L, nobj-1 do begin
      junk = where(abs(objflux[*,iobj]) * sqrt(objivar[*,iobj]) GT 200., ct)
      if (ct GT 0) then $
       splog, 'WARNING: Fiber #', fiberid[iobj], $
        ' has ', ct, ' pixels with S/N > 200'

      junk = where(objflux[*,iobj] * sqrt(objivar[*,iobj]) LE -10., ct)
      if (ct GT 0) then $
       splog, 'WARNING: Fiber #', fiberid[iobj], $
        ' has ', ct, ' pixels with Flux < -10*Noise'
   endfor

   ;----------
   ; Mask out points that are unphysically negative (10-sigma negatives),
   ; and mask the neighboring 2 pixels in each direction.

   for iobj=0L, nobj-1 do begin
      thismask = objflux[*,iobj] * sqrt(objivar[*,iobj]) LE -10.
      thismask = smooth(float(thismask),5) GT 0
      objivar[*,iobj] = objivar[*,iobj] * (1 - thismask)
   endfor

   ;----------
   ; Find GALAXY redshifts

   npoly = 3
   zmin = -0.01 ; -3000 km/sec
   zmax = 0.60 ; Max z for a rest-frame template to 2300 Ang to cover 3700 Ang
   pspace = 2
   nfind = 5
   plottitle = 'Galaxy Redshift'

   eigenfile = 'spEigenGal-*.fits'

   splog, 'Compute GALAXY redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_gal = zfind(objflux, objivar, hdr=hdr, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=5*pspace, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to compute GALAXY redshifts = ', systime(1)-t0

   splog, 'Locally re-fitting GALAXY redshifts'
   t0 = systime(1)
   res_gal = zrefind(objflux, objivar, hdr=hdr, $
    pwidth=5, pspace=1, width=5, zold=res_gal, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to re-fit GALAXY redshifts = ', systime(1)-t0

   ; Only solve for velocity dispersions for the best-fit
   splog, 'Find velocity dispersions for galaxies'
   t0 = systime(1)
   ifind = 0
   vdans = vdispfit(objflux, objivar, hdr=hdr, zobj=res_gal[ifind,*].z, $
    eigenfile='spEigenElodie.fits', columns=lindgen(24), yfit=dispflux)
   res_gal[ifind,*].vdisp = reform([vdans.vdisp],1,nobj)
   res_gal[ifind,*].vdisp_err = reform([vdans.vdisp_err],1,nobj)
   res_gal[ifind,*].vdispchi2 = reform([vdans.vdispchi2],1,nobj)
   res_gal[ifind,*].vdispnpix = reform([vdans.vdispnpix],1,nobj)
   res_gal[ifind,*].vdispdof = reform([vdans.vdispdof],1,nobj)
   splog, 'CPU time to fit GALAXY velocity dispersions = ', systime(1)-t0

   res_gal.class = 'GALAXY'
   res_gal.subclass = ' '

   res_all = res_gal ; Append results

   ;----------
   ; Find QSO redshifts

   npoly = 3
   zmin = 0.0033 ; +1000 km/sec
   zmax = 7.00 ; Max range to use for now, with the template starting at
               ; 460 Ang (rest), which corresponds to 3680 Ang at this z.
   pspace = 4
   nfind = 5
   plottitle = 'QSO Redshift'

   eigenfile = 'spEigenQSO-*.fits'

   splog, 'Compute QSO redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_qso = zfind(objflux, objivar, hdr=hdr, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=7*pspace, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to compute QSO redshifts = ', systime(1)-t0

   splog, 'Locally re-fitting QSO redshifts'
   t0 = systime(1)
   res_qso = zrefind(objflux, objivar, hdr=hdr, $
    pwidth=11, pspace=1, width=11, zold=res_qso, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to re-fit QSO redshifts = ', systime(1)-t0

   res_qso.class = 'QSO'
   res_qso.subclass = ' '

   res_all = [res_all, res_qso] ; Append results

   ;----------
   ; Find STAR redshifts

   npoly = 4 ; With only 1 eigen-template, fit more polynomial terms for stars.
   zmin = -0.004 ; -1200 km/sec
   zmax = 0.004 ; +1200 km/sec
   pspace = 1
   nfind = 1

   eigenfile = 'spEigenStar-*.fits'

   ; Select the stars eigen-file here to detemine how many templates are in it
   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   if (ct EQ 0) then $
    message, 'Unable to find EIGENFILE matching '+eigenfile
   eigenfile = fileandpath(allfiles[ (reverse(sort(allfiles)))[0] ])
   shdr = headfits(djs_filepath(eigenfile, root_dir=eigendir))
   nstar = sxpar(shdr, 'NAXIS2') > 1

   for istar=0, nstar-1 do begin
      subclass = strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)
      plottitle = subclass + '-Star Redshift'

      splog, 'Compute STAR (' + subclass + ') redshifts:', $
       ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
      t0 = systime(1)
      res_star = zfind(objflux, objivar, hdr=hdr, $
       eigenfile=eigenfile, columns=istar, npoly=npoly, $
       zmin=zmin, zmax=zmax, pspace=1, nfind=nfind, width=5*pspace, $
       plottitle=plottitle, doplot=doplot, debug=debug)
      splog, 'CPU time to compute STAR redshifts = ', systime(1)-t0

      res_star.class = 'STAR'
      res_star.subclass = subclass

      res_all = [res_all, res_star] ; Append results
   endfor

   ;----------
   ; Find CV STAR redshifts

   npoly = 3
   zmin = -0.0033 ; -1000 km/sec
   zmax = 0.0033 ; +1000 km/sec
   pspace = 1
   nfind = 1

   eigenfile = 'spEigenCVstar-*.fits'

   subclass = 'CV'
   plottitle = subclass + '-Star Redshift'

   splog, 'Compute STAR (' + subclass + ') redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_cvstar = zfind(objflux, objivar, hdr=hdr, $
    eigenfile=eigenfile, npoly=npoly, $
    zmin=zmin, zmax=zmax, pspace=1, nfind=nfind, width=5*pspace, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to compute STAR redshifts = ', systime(1)-t0

   res_cvstar.class = 'STAR'
   res_cvstar.subclass = subclass

   res_all = [res_all, res_cvstar] ; Append results

   ;----------
   nper = (size(res_all,/dimens))[0]

   ;----------
   ; Sort results for each object by ascending order in chi^2/DOF,
   ; but putting any results with zero degrees-of-freedom at the end.

   minvdiff = 1000.0 ; km/s
   cspeed = 2.99792458e5

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
   ; Generate the synthetic spectra, and count the fraction of points
   ; that deviate more than N sigma (where N goes from 1 to NFSIG).

   t0 = systime(1)
   nfsig = 10
   chi68p = fltarr(nper,nobj)
   fracnsigma = fltarr(nfsig,nper,nobj)
   fracnsighi = fltarr(nfsig,nper,nobj)
   fracnsiglo = fltarr(nfsig,nper,nobj)
   counts_spectro = fltarr(5,nper,nobj)
   counts_synth = fltarr(5,nper,nobj)

   objloglam = objloglam0 + lindgen(npixobj) * objdloglam
   wavevec = 10d^objloglam
   flambda2fnu = wavevec^2 / 2.99792e18

   fthru = filter_thru(objflux * rebin(flambda2fnu,npixobj,nobj), $
    waveimg=wavevec, mask=(objivar EQ 0))
   counts_spectro[*,0,*] = transpose(fthru) * 10^((48.6 - 2.5*17.)/2.5)

   ; Loop in reverse order, so that we look at the best-fit spectra last,
   ; and keep those spectra around for later.

; Save time for now and only look at best fit, since SYNTHSPEC and
; FILTER_THRU are so slow ???
;   for iper=nper-1, 0, -1 do begin
   for iper=0, 0, -1 do begin
      ; Copy this for all fits, since the measured magnitudes are the same
      counts_spectro[*,iper,*] = counts_spectro[*,0,*]

      synflux = synthspec(res_all[iper,*], loglam=objloglam)

      for iobj=0, nobj-1 do begin
         igood = where(objivar[*,iobj] GT 0, ngood)
         if (ngood GT 0) then begin
            chivec = (objflux[igood,iobj] - synflux[igood,iobj]) $
             * sqrt(objivar[igood,iobj])
            abschivec = abs(chivec)
            chi68p[iper,iobj] = (abschivec[sort(abschivec)])[floor(0.68*ngood)]
            for isig=0, nfsig-1 do begin
               fracnsigma[isig,iper,iobj] = total(abschivec GT isig+1) / ngood
               fracnsighi[isig,iper,iobj] = total(chivec GT isig+1) / ngood
               fracnsiglo[isig,iper,iobj] = total(chivec LT (-isig-1)) / ngood
            endfor
         endif
      endfor

      fthru = filter_thru(synflux * rebin(flambda2fnu,npixobj,nobj), $
       waveimg=wavevec)
      counts_synth[*,iper,*] = transpose(fthru) * 10^((48.6 - 2.5*17.)/2.5)
   endfor

   splog, 'CPU time to generate chi^2 statistics = ', systime(1)-t0

   ;----------
   ; Zero-out the dispersion template if the best-fit was not a galaxy.

   for iobj=0, nobj-1 do begin
      if (strtrim(res_gal[iobj].class,2) NE 'GALAXY') then $
       dispflux[*,iobj] = 0
   endfor

   ;----------
   ; Add other fields to the output structure

   splog, 'Adding other fields to output structure'
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
            zwarning:  0L, $
            sn_median: 0.0, $
            chi68p: 0.0, $
            fracnsigma: fltarr(nfsig), $
            fracnsighi: fltarr(nfsig), $
            fracnsiglo: fltarr(nfsig), $
            counts_spectro: fltarr(5), $
            counts_synth: fltarr(5), $
            counts_sky: fltarr(5), $
            anyandmask: 0L, $
            anyormask:  0L, $
            spec1_g:   sxpar(hdr, 'SPEC1_G'), $
            spec1_r:   sxpar(hdr, 'SPEC1_R'), $
            spec1_i:   sxpar(hdr, 'SPEC1_I'), $
            spec2_g:   sxpar(hdr, 'SPEC2_G'), $
            spec2_r:   sxpar(hdr, 'SPEC2_R'), $
            spec2_i:   sxpar(hdr, 'SPEC2_I') }
   res_append = make_array(value=res1, dimension=size(res_all,/dimens))
   res_all = struct_addtags(res_all, res_append)

   for iobj=0, nobj-1 do begin
      igood = where(objivar[*,iobj] NE 0, ngood)
      res_all[*,iobj].wavemin = $
       10^(objloglam0 + (igood[0]>0)*objdloglam) * (ngood NE 0)
      res_all[*,iobj].wavemax = $
       10^(objloglam0 + (igood[(ngood-1)>0])*objdloglam) * (ngood NE 0)
      res_all[*,iobj].wcoverage = ngood * objdloglam
      res_all[*,iobj].anyandmask = anyandmask[iobj]
      res_all[*,iobj].anyormask = anyormask[iobj]
      if (ngood GT 0) then $
       res_all[*,iobj].sn_median = $
        median( sqrt(objivar[igood,iobj]) * abs(objflux[igood,iobj]))
   endfor

   res_all.chi68p = chi68p
   res_all.fracnsigma = fracnsigma
   res_all.fracnsighi = fracnsighi
   res_all.fracnsiglo = fracnsiglo
   res_all.counts_spectro = counts_spectro
   res_all.counts_synth = counts_synth

   ;----------
   ; Generate output headers for spZbest, spZall, spZline files.

   sxaddpar, hdr, 'NAXIS', 0
   sxdelpar, hdr, 'NAXIS1'
   sxdelpar, hdr, 'NAXIS2'
   sxaddpar, hdr, 'EXTEND', 'T', after='NAXIS'
   sxaddpar, hdr, 'VERS1D', idlspec2d_version(), $
    ' Version of idlspec2d for 1D reduction', after='VERSCOMB'
   spawn, 'uname -n', uname
   sxaddpar, hdr, 'UNAME', uname[0]
   ww = strsplit(uname[0], '.', /extract)
   if (ww[1<(n_elements(ww)-1)] EQ 'fnal') then return

   ;----------
   ; Call the line-fitting code for this plate

   splog, 'Call line-fitting code'

; Should be equivalent ???
;   speclinefit, platefile, fiberid=fiberid, $
;    zhdr=hdr, zans=(res_all[0,*])[*], synflux=synflux, dispflux=dispflux, $
;    zline=zline, doplot=doplot, debug=debug

   speclinefit, fiberid=fiberid, $
    hdr=hdr, objflux=objflux, objivar=objivar, $
    zhdr=hdr, zans=(res_all[0,*])[*], synflux=synflux, dispflux=dispflux, $
    outfile=zlinefile, zline=zline, doplot=doplot, debug=debug

   ;----------
   ; Classify galaxies and QSO's based upon emission lines:
   ;   log10(OIII/Hbeta) > 0.7 - 1.2 * (log10(NII/Halpha) - 0.4)  AGN
   ;                     <                                        STARFORMING
   ; If the H_alpha E.W. > 500 Ang, then upgrade STARFORMING -> STARBURST.
   ; If any galaxies or quasars have lines detected at the 10-sigma level
   ;   with sigmas > 200 km/sec at the 5-sigma level, call them BROADLINE.

   nline = (size(zline, /dimens))[0]
   i5007 = where(strtrim(zline.linename,2) EQ '[O_III] 5007')
   ihbeta = where(strtrim(zline.linename,2) EQ 'H_beta')
   ihalpha = where(strtrim(zline.linename,2) EQ 'H_alpha')
   i6583 = where(strtrim(zline.linename,2) EQ '[N_II] 6583')

   q_good = zline[i5007].linearea_err GT 0 $
    AND zline[ihbeta].linearea_err GT 0 $
    AND zline[ihalpha].linearea_err GT 0 $
    AND zline[i6583].linearea_err GT 0
   q_good = q_good $
    AND  zline[i5007].linearea GT 3 * zline[i5007].linearea_err $
    AND zline[ihbeta].linearea GT 3 * zline[ihbeta].linearea_err $
    AND zline[ihalpha].linearea GT 3 * zline[ihalpha].linearea_err $
    AND zline[i6583].linearea GT 3 * zline[i6583].linearea_err
   q_agn = zline[i5007].linearea * (zline[i6583].linearea)^(1.2) $
    GT 10^(0.22) * zline[ihbeta].linearea * (zline[ihalpha].linearea)^(1.2)
   q_obj = strtrim((res_all[0,*].class)[*],2) EQ 'GALAXY' $
    OR strtrim((res_all[0,*].class)[*],2) EQ 'QSO'
   q_stronghalpha = zline[ihalpha].lineew GT 50 $
    AND zline[ihalpha].lineew_err GT 0 $
    AND zline[ihalpha].lineew GT 3 * zline[ihalpha].lineew_err

   ; Find the maximum of (sigma - 5*sigma_err) for all lines of each object
   ; Insist that the lines be detected at the 10-sigma level.
   maxsigma = fltarr(nobj)
   for iobj=0, nobj-1 do $
    for iline=0, nline-1 do $
     if (strtrim(zline[iline,iobj].linename,2) NE 'Ly_alpha' $
      AND zline[iline,iobj].linearea GT 10*zline[iline,iobj].linearea_err $
      AND zline[iline,iobj].linesigma_err GT 0) then $
       maxsigma[iobj] = maxsigma[iobj] > $
        (zline[iline,iobj].linesigma - 5*zline[iline,iobj].linesigma_err)

   indx = where(q_good AND q_obj AND q_agn)
   if (indx[0] NE -1) then res_all[0,indx].subclass $
    = strtrim(res_all[0,indx].subclass + ' AGN', 2)

   indx = where(q_good AND q_obj AND (q_agn EQ 0) AND (q_stronghalpha EQ 0))
   if (indx[0] NE -1) then res_all[0,indx].subclass $
    = strtrim(res_all[0,indx].subclass + ' STARFORMING', 2)

   indx = where(q_good AND q_obj AND (q_agn EQ 0) AND (q_stronghalpha EQ 1))
   if (indx[0] NE -1) then res_all[0,indx].subclass $
    = strtrim(res_all[0,indx].subclass + ' STARBURST', 2)

   indx = where(q_obj AND maxsigma GT 200.)
   if (indx[0] NE -1) then res_all[0,indx].subclass $
    = strtrim(res_all[0,indx].subclass + ' BROADLINE', 2)

   ;----------
   ; Find the best-fit Elodie star for all objects classified as stars

   fitindx = where(strtrim(res_all[0,*].class,2) EQ 'STAR', nfit)
   splog, 'Fitting to Elodie spectra for ', nfit, ' stars'
   t0 = systime(1)

   res_elodie = elodie_best(objflux, objivar, hdr=hdr, fitindx=fitindx)

   splog, 'CPU time to fit to Elodie = ', systime(1)-t0

   ;----------
   ; Set ZWARNING flags.

   splog, 'Setting flags'
   zwarning = lonarr(nper,nobj)

   ; Warning: Sky fiber.
   for iobj=0, nobj-1 do begin
      if (strtrim(plugmap[iobj].objtype,2) EQ 'SKY') then $
       zwarning[*,iobj] = zwarning[*,iobj] OR sdss_flagval('ZWARNING', 'SKY')
   endfor

   ; Warning: too little wavelength coverage.
   qflag = res_all.wcoverage LT 0.18
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'LITTLE_COVERAGE')

   ; Warning: delta-chi^2 is too small as compared to the next best ID.
   minrchi2diff = 0.01
   qflag = res_all.rchi2diff LT minrchi2diff $
    OR res_all.rchi2diff LT minrchi2diff * res_all.rchi2
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'SMALL_DELTA_CHI2')

   ; Warning: synthetic spectrum is negative (for STAR only).
   qflag = (strtrim(res_all.class) EQ 'STAR' $
    AND strtrim(res_all.subclass) NE 'CV' $
    AND res_all.theta[0] LE 0)
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'NEGATIVE_MODEL')

   ; Warning: Fraction of points above 5 sigma is too large (> 5%),
   ; except for QSO's where we just look at the fraction of high outliers
   ; since we expect absorption lines that could give many low outliers.
   qflag = (strtrim(res_all.class) NE 'QSO' AND fracnsigma[4,*,*] GT 0.05) $
    OR (strtrim(res_all.class) EQ 'QSO' AND fracnsighi[4,*,*] GT 0.05)
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'MANY_OUTLIERS')

   ; Warning: Redshift-error warning flag set to -1, which means that
   ; the chi^2 minimum was at the edge of the redshift-fitting range.
;   qflag = res_all.z_err EQ -1
;   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'Z_FITLIMIT')

   ; Warning: For QSOs, if C_IV, CIII], Mg_II, H_beta or H_alpha are negative
   ; and have at least a few pixels in the fit (DOF > 2).
   for iobj=0, nobj-1 do begin
      if (strtrim(res_all[0,iobj].class,2) EQ 'QSO') then begin
         indx = where(zline.fiberid EQ res_all[0,iobj].fiberid $
          AND (zline.linename EQ 'C_IV 1549' $
            OR zline.linename EQ 'C_III] 1908' $
            OR zline.linename EQ 'Mg_II 2799' $
            OR zline.linename EQ 'H_beta' $
            OR zline.linename EQ 'H_alpha') )
         if (indx[0] NE -1) then begin
            qflag = total(zline[indx].linearea LT 0 $
             AND zline[indx].linearea_err GT 0 $
             AND zline[indx].linedof GT 2) NE 0
            zwarning[0,iobj] = zwarning[0,iobj] $
             OR qflag * sdss_flagval('ZWARNING', 'NEGATIVE_EMISSION')
         endif
      endif
   endfor

   res_all.zwarning = zwarning

   ;----------
   ; Write the output files

   splog, 'Writing output files'
   sxaddpar, hdr, 'NAXIS', 0
   sxdelpar, hdr, 'NAXIS1'
   sxdelpar, hdr, 'NAXIS2'
   sxaddpar, hdr, 'EXTEND', 'T', after='NAXIS'
   sxaddpar, hdr, 'VERS1D', idlspec2d_version(), $
    'Version of idlspec2d for 1D reduction', after='VERSCOMB'
   spawn, 'uname -n', uname
   sxaddpar, hdr, 'UNAME', uname[0]
   ww = strsplit(uname[0], '.', /extract)
   if (ww[1<(n_elements(ww)-1)] EQ 'fnal') then return

   zans = struct_addtags((res_all[0,*])[*], res_elodie)
   mwrfits, 0, zbestfile, hdr, /create ; Retain the original header in first HDU
   mwrfits, zans, zbestfile
   mwrfits, synflux, zbestfile
   mwrfits, dispflux, zbestfile

   sxaddpar, hdr, 'DIMS0', nper, ' Number of fits per objects'
   sxaddpar, hdr, 'DIMS1', nobj, ' Number of objects'
   mwrfits, 0, zallfile, hdr, /create ; Retain the original header in first HDU
   mwrfits, res_all, zallfile

   if (keyword_set(debugfile)) then dfpsclose

   ;----------
   ; Generate final QA plots

   splog, 'Generating QA plots'

   if (keyword_set(plotfile)) then begin
      cpbackup, plotfile
      dfpsplot, plotfile, /color
   endif

   plottitle = string(zans[0].plate, zans[0].mjd, $
    format='("Flux-Calibration Errors Plate=", i4, " MJD=", i5)')
   qaplot_fcalibvec, objloglam, objflux, objivar, synflux, plugmap, zans, $
    plottitle=plottitle

   if (keyword_set(plotfile)) then dfpsclose

   ;----------
   ; Close log file

   splog, 'Total time for SPREDUCE1D = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE1D at ' + systime()
   if (keyword_set(logfile)) then splog, /close

   return
end
;------------------------------------------------------------------------------
