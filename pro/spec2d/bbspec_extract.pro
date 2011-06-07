;+
; NAME:
;   bbspec_extract
;
; PURPOSE:
;   Run bbspec 2D extraction code
;
; CALLING SEQUENCE:
;   bbspec_extract, image, invvar, flux, fluxivar, basisfile=, $
;    [ ximg=, frange=, yrange=, ymodel=, tmproot=, outfile=, /batch ]
;
; INPUTS:
;   image      - Image [NX,NY]
;   innvar     - Inverse variance image corresponding to IMAGE [NX,NY]
;   basisfile  - File with PSF model for bbspec
;
; OPTIONAL INPUTS:
;   ximg       - X centroids of fibers on IMAGE, to replace those entries
;                in the PSF file [NY,NFIBER]
;   frange     - Fiber number range (0-indexed); default to all fibers
;                represented in the PSF file; extract all 20-fiber bundles
;                spanned by this range, so specifying [25,30] would extract
;                fiber numbers [20,39]
;   yrange     - 0-indexed ange of rows to extract; default to all rows
;                that contain any unmasked pixels
;   tmproot    - Root file name for temporary files; default to 'tmp-';
;                necessary to be a unique string identifier if multiple
;                instances of this procedure running in the same directory,
;                for example those spawned by BBSPECT_TEST.
;   outfile    - If set, then output FITS file with 3 HDUs containing FLUX,
;                FLUXIVAR, YMODEL
;   batch      - If set, then batch each bundle of 20 fibers to a different
;                PBS job; FRANGE keyword cannot also be set
;
; OUTPUTS:
;   flux       - Extracted flux vectors [NY,NFIBER]
;   fluxivar   - Extracted inverse variance vectors [NY,NFIBER]
;
; OPTIONAL OUTPUTS:
;   ymodel     - Model-fit image
;
; COMMENTS:
;   Currently hard-wired to extract in chunks of 20 fibers by 100 rows,
;   with the 15 rows on the top and bottom to be used as padding and
;   discarded.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-May-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro bbspec_extract_readpsf, basisfile1, pbasis, phdr

   ; Determine the number of HDUs in the PSF file
   nhdu = 0
   basisfile = findfile(basisfile1, count=ct)
   if (ct EQ 0) then $
    message, 'PSF file not found '+string(basisfile1)
   while (size(headfits(basisfile,exten=nhdu,/silent),/tname) EQ 'STRING') do $
    nhdu++
   if (nhdu EQ 0) then $
    message, 'Error reading PSF file '+string(basisfile)

   phdr = ptrarr(nhdu)
   pbasis = ptrarr(nhdu)
   for ihdu=0, nhdu-1 do begin
      pbasis[ihdu] = ptr_new(mrdfits(basisfile, ihdu, hdr1))
      phdr[ihdu] = ptr_new(hdr1)
   endfor

   return
end
;------------------------------------------------------------------------------
pro bbspec_extract_shiftpsf, pbasis, phdr, psffile, ximg=ximg

   psftype = strtrim(sxpar(*phdr[0],'PSFTYPE'),2)
   ny = sxpar(*phdr[0],'NAXIS1')
   nfiber = sxpar(*phdr[0],'NAXIS2')

   for ihdu=0, n_elements(pbasis)-1 do begin
      basis1 = *pbasis[ihdu]
      hdr1 = *phdr[ihdu]

      ; Replace with the X centroids shifted, and trim to only the first entries
      ; if the PSF is only solved for the first fibers in the first rows
      if (ihdu EQ 0 AND keyword_set(ximg)) then basis1 = ximg[0:ny-1,0:nfiber-1]
 
      mwrfits, basis1, psffile, hdr1, create=(ihdu EQ 0), /silent
   endfor

   return
end
;------------------------------------------------------------------------------
function bbspec_cmd_batch, cmd, pbsfile=pbsfile

   if (NOT keyword_set(pbsfile)) then $
    message, 'PBSFILE not set'
   openw, olun, pbsfile, /get_lun
   printf, olun, '# Auto-generated batch file '+systime()
   printf, olun, '#PBS -l nodes=1'
   printf, olun, '#PBS -l walltime=48:00:00'
   printf, olun, '#PBS -W umask=0022'
   printf, olun, '#PBS -V'
   printf, olun, '#PBS -j oe'
   printf, olun, 'cd $PBS_O_WORKDIR'
   printf, olun, 'set -o verbose'
   printf, olun, cmd
   close, olun
   free_lun, olun

   splog, 'Submitting file '+pbsfile
   spawn, 'qsub '+pbsfile, jobid

   return, jobid
end
;------------------------------------------------------------------------------
pro bbspec_extract, image, invvar, flux, fluxivar, basisfile=basisfile, $
 ximg=ximg, frange=frange1, yrange=yrange1, tmproot=tmproot1, ymodel=ymodel, $
 outfile=outfile, batch=batch

   if (n_params() NE 4 OR keyword_set(basisfile) EQ 0) then $
    message, 'Parameters not set'
   if (keyword_set(yrange1)) then yrange = yrange1 $
    else yrange = minmax(where(total(invvar,1) GT 0))
   if (keyword_set(tmproot1)) then tmproot = tmproot1 $
    else tmproot = 'tmp-'
   if (keyword_set(batch) AND keyword_set(frange1)) then $
    message, 'Cannot specify both BATCH and FRANGE'
   imgfile = tmproot+'img.fits'

   stime0 = systime(1)

   ; Write the image file
   mwrfits, image, imgfile, /create, /silent
   mwrfits, invvar, imgfile, /silent

   ; If XIMG is set, then shift the PSF traces
   if (keyword_set(ximg)) then begin
      psffile = tmproot+'psf.fits'
      bbspec_extract_readpsf, basisfile, pbasis, phdr
      bbspec_extract_shiftpsf, pbasis, phdr, psffile, ximg=ximg
   endif else begin
      psffile = basisfile
   endelse
   psfhdr = headfits(psffile)
   ny = sxpar(psfhdr,'NAXIS1')
   nfiber = sxpar(psfhdr,'NAXIS2')
   if (keyword_set(frange1)) then frange = frange1 $
    else frange = [0,nfiber-1]
   brange = [frange[0]/20,ceil(frange[1]/20)] ; range of bundle numbers

   pyfile = djs_filepath('extract_spectra.py', root_dir=getenv('BBSPEC_DIR'), $
    subdir='bin')
   fstr = ' -f '+strtrim(yrange[0],2)+','+strtrim(yrange[1],2)+',100'

   ; Extract the spectra
   if (NOT keyword_set(batch)) then begin
      fluxfile = tmproot+'flux.fits'
      cmd = pyfile+' -i '+imgfile+' -p '+psffile+' -o '+fluxfile+fstr
      if (keyword_set(brange)) then begin
         cmd += ' -b '+strtrim(brange[0],2)
         for j=brange[0]+1, brange[1] do cmd += ','+strtrim(j,2)
      endif
      splog, 'SPAWN ', cmd
      spawn, cmd, res, errcode
      flux = mrdfits(fluxfile,0)
      fluxivar = mrdfits(fluxfile,1)
      ymodel = mrdfits(fluxfile,5)
;      rmfile, fluxfile
   endif else begin
      njob = brange[1] - brange[0] + 1
      bvec = brange[0] + lindgen(njob)
      jobid = lonarr(njob)
      splog, 'Batching ', njob, ' jobs'
      fluxfile = tmproot+'flux'+'-'+string(bvec,format='(i2.2)')+'.fits'
      pbsfile = tmproot+'flux'+'-'+string(bvec,format='(i2.2)')
      cmd = pyfile+' -i '+imgfile+' -p '+psffile+' -o '+fluxfile+fstr $
       + ' -b '+strtrim(bvec,2)
      for ijob=0, njob-1 do begin
         jobid[ijob] = bbspec_cmd_batch(cmd[ijob], pbsfile=pbsfile[ijob])
      endfor
      bbspec_batch_wait, jobid
      ymodel = 0
      flux = 0
      fluxivar = 0
      for i=0, njob-1 do flux += mrdfits(fluxfile[i],0)
      for i=0, njob-1 do fluxivar += mrdfits(fluxfile[i],1)
      for i=0, njob-1 do ymodel += mrdfits(fluxfile[i],5)
;      for i=0, njob-1 do rmfile, fluxfile[i]
   endelse

   for i=0, n_elements(phdr)-1 do ptr_free, phdr[i]
   for i=0, n_elements(pbasis)-1 do ptr_free, pbasis[i]

   splog, 'Time to bbspec = ', systime(1)-stime0, ' seconds'

   return
end
;------------------------------------------------------------------------------
