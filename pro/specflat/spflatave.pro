;+
; NAME:
;   spflatave
;
; PURPOSE:
;   Average together a set (or all!) 2D pixel flats.
;
; CALLING SEQUENCE:
;   spflatave, [ mjd=, mjstart=, mjend=, mjout=, indir=, outdir=, docam= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - Valid MJD's for input pixel flats.
;   mjstart    - Valid starting MJD for input pixel flats.
;   mjend      - Valid ending MJD for input pixel flats.
;   mjout      - MJD for name of output average pixel flat; default to 0,
;                resulting in file names like 'pixflatave-00000-b1.fits'.
;   indir      - Input directory; default to current directory.
;   outdir     - Output directory; default to same as INDIR.
;   docam      - Camera names; default to all cameras: ['b1', 'b2', 'r1', 'r2']
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Some sigma-clipping is done before combining each pixel, clipping
;   at 2-sigma if more than 7 frames, and a lower sigma for fewer frames.
;
;   The output file has two HDU's, the first being the average flat,
;   the second being the standard deviation at each pixel.
;     --> Comment this out for the time being!!???
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   djs_iterstat
;   fileandpath()
;   headfits()
;   mrdfits()
;   splog
;   sxaddpar
;   sxpar()
;   writefits
;
; REVISION HISTORY:
;   17-Jul-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro spflatave, mjd=mjd, mjstart=mjstart, mjend=mjend, mjout=mjout, $
 indir=indir, outdir=outdir, docam=docam

   if (NOT keyword_set(indir)) then indir = ''
   if (NOT keyword_set(outdir)) then outdir = indir
   if (NOT keyword_set(mjout)) then mjout = 0L

   if (keyword_set(docam)) then camnames = docam $
    else camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   for icam=0, ncam-1 do begin

      ; Find all input pixel flats for this camera that match the
      ; specified input MJD's.
      files = findfile(djs_filepath('pixflat-*-'+camnames[icam]+'.fits', $
       root_dir=indir), count=nfile)
      if (nfile GT 0) then begin
         thismjd = long(strmid(fileandpath(files),8,5))
         qkeep = bytarr(nfile) + 1
         if (keyword_set(mjstart)) then $
          qkeep = qkeep AND (thismjd GE mjstart)
         if (keyword_set(mjend)) then $
          qkeep = qkeep AND (thismjd LE mjend)
         if (keyword_set(mjd)) then $
          for ifile=0, nfile-1 do $
           qkeep[ifile] = qkeep[ifile] $
            AND (total(thismjd[ifile] EQ long(mjd)) NE 0)
         ikeep = where(qkeep, nfile)
         if (nfile GT 0) then files = files[ikeep]
      endif

      splog, 'Found ' + string(nfile) + ' files for camera ' + camnames[icam]

      if (nfile GT 0) then begin
         ;----------
         ; Read all the images into a single array

         hdr = headfits(files[0])
         naxis1 = sxpar(hdr,'NAXIS1')
         naxis2 = sxpar(hdr,'NAXIS2')
         npix = naxis1 * naxis2
         pixflatarr = fltarr(naxis1, naxis2, nfile)
         for ifile=0, nfile-1 do $
          pixflatarr[*,*,ifile] = mrdfits(files[ifile])

         ;----------
         ; For each MJD's pixel flat, remove the first two and last two
         ; nonzero pixels in each row.

         for ifile=0, nfile-1 do begin
            for iy=0, naxis2-1 do begin
               igood = where(pixflatarr[*,iy,ifile] GT 0, ngood)
               if (ngood GT 0) then begin
                  pixflatarr[igood[0]:(igood[0]+1)<(naxis1-1),iy,ifile] = 0
                  pixflatarr[(igood[ngood-1]-1)>0:igood[ngood-1],iy,ifile] = 0
               endif
            endfor
         endfor

         ;----------
         ; Generate a map of the sigma at each pixel (doing some rejection).
         ; This is a horrible loop over each pixel, but shouldn't take more
         ; than a few minutes.

         if (nfile LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 flats
          else if (nfile EQ 3) then sigrej = 1.1 $
          else if (nfile EQ 4) then sigrej = 1.3 $
          else if (nfile EQ 5) then sigrej = 1.6 $
          else if (nfile EQ 6) then sigrej = 1.9 $
          else sigrej = 2.0
         maxiter = 2

         aveimg = fltarr(naxis1, naxis2)
         sigimg = fltarr(naxis1, naxis2)
         for ipix=0L, npix-1 do begin
            ; Only consider pixels with values > 0
            ; The SPFLATTEN2 routine will have set values = 0 where there
            ; were not enough counts to determine the flat-field value.
            vals = pixflatarr[lindgen(nfile)*npix+ipix]
            ii = where(vals GT 0)
            if (ii[0] NE -1) then begin
               djs_iterstat, vals[ii], $
                sigrej=sigrej, maxiter=maxiter, sigma=sigma1, mean=mean1
               aveimg[ipix] = mean1
               sigimg[ipix] = sigma1
            endif
         endfor

         ;----------
         ; Reject pixels in the average flat where the dispersion between
         ; the input flats was large

;         maskimg = sigimg LT 0.05
;         junk = where(maskimg EQ 0, nbad)
;         splog, 'Reject ', nbad, ' pixels'
;
;         aveimg = aveimg * maskimg

         ;----------
         ; Keep pixels near the left and right edge that didn't have
         ; the minimum counts in any of the pixel flats, and set those
         ; pixels equal to unity.

         ; MASKSUM will equal zero for any pixels with no good flats
         if (nfile EQ 1) then $
          masksum = (pixflatarr GT 0) $
         else $
          masksum = total(pixflatarr GT 0, 3)
pixflatarr = 0 ; Clear memory

         for iy=0, naxis2-1 do begin
            igood = where(masksum[*,iy] GT 0 AND aveimg[*,iy] GT 0, ngood)
            ; Don't do anything if the whole row is bad
            if (ngood GT 0) then begin
               ; Set left-most pixels to unity
               if (igood[0] NE 0) then $
                aveimg[0:igood[0]-1,iy] = 1.0
               ; Set right-most pixels to unity
               if (igood[ngood-1] NE naxis1-1) then $
                aveimg[igood[ngood-1]+1:naxis1-1,iy] = 1.0
            endif
         endfor

         ;----------
         ; Append comments to the header

         for ifile=0, nfile-1 do $
          sxaddpar, hdr, 'COMMENT', 'Include file ' + fileandpath(files[ifile])

         ;----------
         ; Write the output file

         outfile = djs_filepath( string(mjout, camnames[icam], $
          format='("pixflatave-",i5.5,"-",a2,".fits")'), root_dir=outdir)

         splog, 'Writing file ' + outfile
         writefits, outfile, aveimg, hdr
;         mwrfits, sigimg, outfile

      endif

   endfor

   return
end
;------------------------------------------------------------------------------
