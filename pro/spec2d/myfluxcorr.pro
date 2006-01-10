pro myfluxcorr, bsmearfile, rsmearfile, bscifile, rscifile, corrfile, $
 adderr=adderr, rset=rset, bset=bset

   ;----------
   ; Read the plug-map file (for identifying sky fibers)

   bsmearflux = mrdfits(bsmearfile,0)
   bsmearivar = mrdfits(bsmearfile,1)
   bsmearmask = mrdfits(bsmearfile,2)
   bsmearset  = mrdfits(bsmearfile,3)
   plugmap    = mrdfits(bsmearfile,5)
   traceset2xy, bsmearset, xx, bsmearloglam
   correct_dlam, bsmearflux, bsmearivar, bsmearset, dlam=1.0e-4
   if keyword_set(bset) then bsmearfluxcalib = bspline_valu(bsmearloglam, bset)
   
   rsmearflux = mrdfits(rsmearfile,0)
   rsmearivar = mrdfits(rsmearfile,1)
   rsmearmask = mrdfits(rsmearfile,2)
   rsmearset  = mrdfits(rsmearfile,3)
   traceset2xy, rsmearset, xx, rsmearloglam
   correct_dlam, rsmearflux, rsmearivar, rsmearset, dlam=1.0e-4
   if keyword_set(rset) then rsmearfluxcalib = bspline_valu(rsmearloglam, rset)
 
   ;---------------------------------------------------------
   ;  now put blue and red together
   ;
      smearwave = [temporary(bsmearloglam),temporary(rsmearloglam)]
      smearflux = [temporary(bsmearflux),temporary(rsmearflux)]
      smearivar = [temporary(bsmearivar),temporary(rsmearivar)]

      if keyword_set(rset) then begin
          smearfluxcalib = [temporary(bsmearfluxcalib),$
                            temporary(rsmearfluxcalib)]
          divideflat, smearflux, invvar=smearivar, smearfluxcalib, minval=0.0001
      endif

      ;--------------------------------------------------------------------
      ;    bkpts from 3700 to 9300
      ;
      wavemin = alog10(3700.0)
      wavemax = alog10(9300.0)
      waverange = wavemax - wavemin
      wavemid   = 0.5*(wavemin + wavemax)
      bkpt = findgen(500)*8.0e-4 + alog10(3700.0)
      ncoeff = 3 ; ???

      nfiber = (size(smearwave,/dimens))[1]

      smearset = ptrarr(nfiber)

      for ifiber=0, nfiber-1 do begin

      print, ifiber

      ; Schlegel counter of step number...
      ; print, format='("Step ",i5," of ",i5,a1,$)', $
      ;   ifiber, nfiber, string(13b)

        sset = 0
        good = where(smearivar[*,ifiber] GT 0, ngood)
        if ngood GT 10 then begin
           inside = where(bkpt GT min(smearwave[good,ifiber]) AND $
                          bkpt LT max(smearwave[good,ifiber]),ninside)

           if ninside GT 1 then begin
             outside = [inside[0]-1,inside,inside[ninside-1]+1]

             sset = bspline_iterfit(smearwave[*,ifiber], smearflux[*,ifiber], $
               invvar=smearivar[*,ifiber], bkpt=bkpt[outside], $
               lower=5, upper=5, requiren=2)

           endif
        endif
        smearset[ifiber] = ptr_new(sset)
      endfor 



   nfiles = n_elements(corrfile)

   for ifile = 0, nfiles - 1 do begin
 
   if (bsmearfile EQ bscifile[ifile] AND $
       rsmearfile EQ rscifile[ifile]) then begin
      ;----------
      ; Special case: If the science image and smear image is the same,
      ; then force their ratio to be unity.

      fitimg = smearwave*0.0 + 1.0
      xy2traceset, smearwave, fitimg, corrset, ncoeff=ncoeff

      corrset.coeff[0,*] = 1.0
      corrset.coeff[1:*,*] = 0.0

   endif else begin

     bsciflux = mrdfits(bscifile[ifile],0)
     bsciivar = mrdfits(bscifile[ifile],1)
     bscimask = mrdfits(bscifile[ifile],2)
     bsciset  = mrdfits(bscifile[ifile],3)
     traceset2xy, bsciset, xx, bsciloglam
     correct_dlam, bsciflux, bsciivar, bsciset, dlam=1.0e-4
     if keyword_set(bset) then bscifluxcalib = bspline_valu(bsciloglam, bset)

     rsciflux = mrdfits(rscifile[ifile],0)
     rsciivar = mrdfits(rscifile[ifile],1)
     rscimask = mrdfits(rscifile[ifile],2)
     rsciset  = mrdfits(rscifile[ifile],3)
     traceset2xy, rsciset, xx, rsciloglam
     correct_dlam, rsciflux, rsciivar, rsciset, dlam=1.0e-4
     if keyword_set(rset) then rscifluxcalib = bspline_valu(rsciloglam, rset)

      sciwave = [temporary(bsciloglam),temporary(rsciloglam)]
      sciflux = [temporary(bsciflux),temporary(rsciflux)] 
      sciivar = [temporary(bsciivar),temporary(rsciivar)]

      if keyword_set(rset) then begin
          scifluxcalib = [temporary(bscifluxcalib),temporary(rscifluxcalib)]
          divideflat, sciflux, invvar=sciivar, scifluxcalib, minval=0.0001
      endif

     
     smearfit = sciwave * 0.0

     for ifiber=0,nfiber - 1 do $
       smearfit[*,ifiber] = bspline_valu(sciwave[*,ifiber],*(smearset[ifiber])) 

     ;----------
     ; Convert these fits into a trace set

     xy2traceset, sciwave, sciflux, corrset, $
        invvar=sciivar, xmin=wavemin, xmax=wavemax, $
        ncoeff=ncoeff, yfit=crazyimg, inputfunc=smearfit

     traceset2xy, corrset, sciwave, fitimg

   ;----------
   ; Which fits are considered bad?
   ; Require the median S/N per pix for the science exposure to be > 2.0.
   ; Require the median S/N per pix for the smear exposure to be > 2.0.
   ; Require the fits for a fiber are never less than 0.2 * median.
   ; Require the fits for a fiber are never greater than 5 * median.
   ; Require the fiber isn't SKY.

     scisnmed = djs_median(sciflux * sqrt(sciivar), 1)
     smearsnmed = djs_median(smearflux * sqrt(smearivar), 1)

     medfit = median(fitimg)
     splog, 'Median flux-correction factor = ', medfit

     qgood = scisnmed GT 3.0 $
      AND strtrim(plugmap.objtype,2) NE 'SKY'

     igood = where(qgood EQ 1, ngood)
     splog, 'Total of', ngood, ' good flux-correction vectors'
     if (igood[0] EQ -1) then $
      message, 'Trouble in flux-correction fits'

     ;----------
     ; Make a mean fit for the good fits
 
     meanfit = total(fitimg[*,igood],2) / ngood
     meanwave = total(sciwave[*,igood],2) / ngood
     xy2traceset, meanwave, meanfit, cset1, ncoeff=ncoeff

     ;----------
     ; Replace bad fits with the mean fit

     ibad = where(qgood EQ 0, nbad)
     splog, 'Replacing', nbad, ' bad flux-correction vectors'
     if nbad GT 0 then $
         corrset.coeff[*,ibad] = cset1.coeff # replicate(1,nbad)

   endelse

   mwrfits, corrset, corrfile[ifile], /create

   endfor

   ptr_free,smearset
 
   return
end
