;+
; NAME:
;   spreduce
;
; PURPOSE:
;   Extract, wavelength-calibrate, and flatten SDSS spectral frame(s).
;
; CALLING SEQUENCE:
;   spreduce, flatname, arcname, objname, [ run2d=, $
;    plugfile=, lampfile=, indir=, plugdir=, outdir=, $
;    ecalibfile=, plottitle=, /do_telluric, writeflatmodel=, writearcmodel=, $
;    /bbspec ]
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;   objname    - Name(s) of object SDSS image(s)
;
; REQUIRED KEYWORDS:
;   plugfile   - Name of plugmap file (Yanny parameter file)
;
; OPTIONAL KEYWORDS:
;   run2d      - 2D reduction name, to include in output headers
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in $IDLSPEC2D_DIR/etc.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   plugdir    - Input directory for PLUGFILE; default to '.'
;   outdir     - Directory for output files; default to '.'
;   ecalibfile - opECalib.par file for SDSSPROC
;   plottitle  - Prefix for titles in QA plots.
;   do_telluric- Passed to EXTRACT_OBJECT
;   writeflatmodel - passed to SPCALIB to write out flat image
;                    model info to file
;   writearcmodel  - passed to SPCALIB to write out arc image
;                    model info to file
;   bbspec         - use bbspec extraction code
;   splitsky       - split sky model between spatial halves
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
;   Should test that arcs and flats are valid images with CHECKFLAVOR.
;
; PROCEDURES CALLED:
;   extract_object
;   fibermask_bits()
;   get_tai
;   qaplot_arcline
;   qaplot_fflat
;   readplugmap()
;   reject_science()
;   select_arc()
;   select_flat()
;   sortplugmap
;   sdssproc
;   spcalib
;   splog
;   sxaddpar
;   sxpar()
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/skylines.dat
;
; REVISION HISTORY:
;   12-Oct-1999  Written by D. Schlegel & S. Burles, APO
;      Apr-2010  Added "write[flat,arc]model" pass-through (A. Bolton, Utah)
;-
;------------------------------------------------------------------------------
pro spreduce, flatname, arcname, objname, run2d=run2d, $
 plugfile=plugfile, lampfile=lampfile, $
 indir=indir, plugdir=plugdir, outdir=outdir, $
 ecalibfile=ecalibfile, plottitle=plottitle, do_telluric=do_telluric, $
 writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, bbspec=bbspec, $
 splitsky=splitsky

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = '.'

   stime0 = systime(1)

   ;---------------------------------------------------------------------------
   ; Locate skyline file for sky wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(skylinefile)) then begin
      fullskyfile = (findfile(skylinefile, count=ct))[0]
   endif else begin
      skydefault = filepath('skylines.dat', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      fullskyfile = (findfile(skydefault, count=ct))[0]
   endelse
   if (NOT keyword_set(fullskyfile)) then $
    message, 'No SKYLINEFILE found '+skylinefile

   ;---------------------------------------------------------------------------
   ; Determine spectrograph ID and color from first object file
   ;---------------------------------------------------------------------------

   sdssproc, objname[0], indir=indir, spectrographid=spectrographid, $
    color=color, ecalibfile=ecalibfile, hdr=objhdr

   ;---------------------------------------------------------------------------
   ; Read PLUGMAP file and sort
   ; Look for photoPlate file in directory OUTDIR
   ;---------------------------------------------------------------------------

   plugmap = readplugmap(plugfile, spectrographid, $
    plugdir=plugdir, /calibobj, mjd=sxpar(objhdr,'MJD'), indir=outdir, $
    exptime=sxpar(objhdr,'EXPTIME'), hdr=hdrplug, fibermask=fibermask)
   if (NOT keyword_set(plugmap)) then begin
      splog, 'ABORT: Plug map not found ' $
       + djs_filepath(plugfile, root_dir=plugdir)
      return
   endif
 
   ;---------------------------------------------------------------------------
   ; REDUCE CALIBRATION FRAMES
   ;---------------------------------------------------------------------------

   flatinfoname = filepath( $
       'spFlat-'+string(format='(a1,i1,a)',color,spectrographid, $
       '-'), root_dir=outdir)

   arcinfoname = filepath( $
       'spArc-'+string(format='(a1,i1,a)',color,spectrographid, $
       '-'), root_dir=outdir)

   heap_gc   ; Garbage collection for all lost pointers

   cartid = long(yanny_par(hdrplug, 'cartridgeId'))
   spcalib, flatname, arcname, fibermask=fibermask, cartid=cartid, $
    lampfile=lampfile, indir=indir, $
    ecalibfile=ecalibfile, plottitle=plottitle, $
    flatinfoname=flatinfoname, arcinfoname=arcinfoname, $
    arcstruct=arcstruct, flatstruct=flatstruct, $
    writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
    bbspec=bbspec

   ;----------
   ; Find the mid-point in time for all of the science observations

   for iobj=0, n_elements(objname)-1 do begin
      objhdr = sdsshead(objname[iobj], indir=indir)
      get_tai, objhdr, tai_beg, tai_mid, tai_end
      if (iobj EQ 0) then begin
         tai1 = tai_beg
         tai2 = tai_end
      endif else begin
         tai1 = tai1 < tai_beg
         tai2 = tai2 > tai_end
      endelse
   endfor
   tai_mid = 0.5 * (tai1 + tai2)

   ;----------
   ; Select one set of flats+arcs for all science exposures.
   ; We *could* use a different flat+arc for each exposure,
   ; which is something we have done in the past.
   ; Instead, now just choose the flat taken nearest the mid-point
   ; of all the science integrations, and take its corresponding arc.

   bestflat = select_flat(flatstruct, tai_mid)
   if (NOT keyword_set(bestflat)) then begin
      splog, 'ABORT: No good flats (saturated?)'
      return
   endif
   splog, 'Best flat = ', bestflat.name

;   bestarc = select_arc(arcstruct)
   bestarc = arcstruct[ bestflat.iarc ]

   if (NOT keyword_set(bestarc)) then begin
      splog, 'ABORT: No good arcs (saturated?)'
      return
   endif

   splog, 'Best arc = ', bestarc.name
   if ((color EQ 'blue' AND bestarc.bestcorr LT 0.5) $
    OR (color EQ 'red'  AND bestarc.bestcorr LT 0.5) ) then begin
      splog, 'ABORT: Best arc correlation = ', bestarc.bestcorr
      return
   endif else $
    if ((color EQ 'blue' AND bestarc.bestcorr LT 0.7) $
    OR (color EQ 'red'  AND bestarc.bestcorr LT 0.7) ) then begin
      splog, 'WARNING: Best arc correlation = ', bestarc.bestcorr
   endif

   lambda = *(bestarc.lambda)
   xpeak = *(bestarc.xpeak)
   wset = *(bestarc.wset)
   dispset = *(bestarc.dispset)

   qaplot_arcline, *(bestarc.xdif_tset), wset, lambda, $
    rejline=*(bestarc.rejline), $
    color=color, title=plottitle+' Arcline Fit for '+bestarc.name

   ; Generate the sdProc files for the arc images and the PSF models
   if (keyword_set(bbspec)) then begin
      sdssproc, bestarc.name, indir=indir, /outfile, $
       /applybias, /applypixflat, /applycrosstalk
      arcstr = strmid(bestarc.name,4,11)
      flatstr = strmid(bestflat.name,4,11)
      ; Assume files exist: sdProc-$arcstr spArc-$arcstr spFlat-$flatstr
      ;- Create PSF file
      pyfile = djs_filepath('make-my-psf.py', root_dir=getenv('BBSPEC_DIR'), $
       subdir='examples')
      cmd = 'python '+pyfile+' '+arcstr+' '+flatstr
      spawn, cmd, res, errcode
      if (keyword_set(errcode)) then begin
         splog, errcode
         message, 'Error calling '+cmd
      endif
      bbspec_pixpsf, arcstr, flatstr, /batch
   endif

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------

   for iobj=0, N_elements(objname)-1 do begin

      stimeobj = systime(1)
      splog, prelog=objname[iobj]

      ;----------
      ; Read object image
      ;   Minflat will mask all pixels with low 2d pixelflat values
 
      splog, 'Reading object ', objname[iobj]
      sdssproc, objname[iobj], image, invvar, indir=indir, hdr=objhdr, $
       /applybias, /applypixflat, spectrographid=spectrographid, color=color, $
       ecalibfile=ecalibfile, minflat=0.8, maxflat=1.2, $
       nsatrow=nsatrow, fbadpix=fbadpix, /applycrosstalk, ccdmask=ccdmask

      ;-----
      ; Decide if this science exposure is bad

      qbadsci = reject_science(image, objhdr, nsatrow=nsatrow, fbadpix=fbadpix)

      sxaddpar, objhdr, 'RUN2D', run2d, ' Spectro-2D reduction name'

      ; In case TAI-BEG,TAI-END were missing from the header, add them in.
      get_tai, objhdr, tai_beg, tai_mid, tai_end
      sxaddpar, objhdr, 'TAI-BEG', tai_beg
      sxaddpar, objhdr, 'TAI-END', tai_end

      sxaddpar, objhdr, 'FRAMESN2', 0.0
      sxaddpar, objhdr, 'TILEID', long(yanny_par(hdrplug, 'tileId')), $
       'Cartridge used in this plugging', after='PLATEID'
      sxaddpar, objhdr, 'CARTID', long(yanny_par(hdrplug, 'cartridgeId')), $
       'Cartridge used in this plugging', after='PLATEID'
      redden = float(yanny_par(hdrplug, 'reddeningMed'))
      if (n_elements(redden) NE 5) then redden = fltarr(5)
      for j=0, n_elements(redden)-1 do $
       sxaddpar, objhdr, string('REDDEN',j+1, format='(a6,i2.2)'), redden[j], $
       ' Median extinction in '+(['u','g','r','i','z'])[j]+'-band'

      if (qbadsci) then begin

         ; We will have already output an abort message in the REJECT_SCIENCE() proc.
         splog, 'Skipping reduction of this bad science exposure'

      endif else if (NOT keyword_set(bestflat)) then begin

         splog, 'ABORT: No good flats'

      endif else begin

         xsol = *(bestflat.xsol)
         fflat = *(bestflat.fflat)
         superflatset = *(bestflat.superflatset)
         widthset = *(bestflat.widthset)
         dispset = *(bestarc.dispset)
         proftype = bestflat.proftype

         sxaddpar, objhdr, 'XSIGMA', max(bestflat.medwidth)
         sxaddpar, objhdr, 'WSIGMA', max(bestarc.medwidth)

         qaplot_fflat, fflat, wset, $
          title=plottitle+'Fiber-Flats for '+bestflat.name

         ;----------
         ; Combine FIBERMASK bits from the plug-map file, best flat
         ; and best arc

         fibermask = fibermask $
          OR (*(bestflat.fibermask) AND fibermask_bits('BADTRACE')) $
          OR (*(bestflat.fibermask) AND fibermask_bits('BADFLAT')) $
          OR (*(bestarc.fibermask) AND fibermask_bits('BADARC'))

         ;----------
         ; Determine output file name and modify the object header

         framenum = sxpar(objhdr, 'EXPOSURE')
         outname = filepath( $
          'spFrame-'+string(format='(a1,i1,a,i8.8,a)',color,spectrographid, $
          '-',framenum,'.fits'), root_dir=outdir)

         sxaddpar, objhdr, 'PLUGFILE', fileandpath(plugfile)
         sxaddpar, objhdr, 'FLATFILE', fileandpath(bestflat.name)
         sxaddpar, objhdr, 'ARCFILE', fileandpath(bestarc.name)
         sxaddpar, objhdr, 'OBJFILE', fileandpath(objname[iobj])
         sxaddpar, objhdr, 'LAMPLIST', fileandpath(lampfile)
         sxaddpar, objhdr, 'SKYLIST', fileandpath(fullskyfile)

         ;-----
         ; Extract the object frame

         extract_object, outname, objhdr, image, invvar, plugmap, wset, $
          xpeak, lambda, xsol, fflat, fibermask, color=color, $
          proftype=proftype, superflatset=superflatset, $
          widthset=widthset, dispset=dispset, skylinefile=fullskyfile, $
          plottitle=plottitle, do_telluric=do_telluric, bbspec=bbspec, $
          splitsky=splitsky, ccdmask=ccdmask

         splog, 'Elapsed time = ', systime(1)-stimeobj, ' seconds', $
          format='(a,f6.0,a)' 
      endelse

   endfor

   heap_gc   ; Garbage collection for all lost pointers
   splog, 'Elapsed time = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)', prelog=''

   return
end
;------------------------------------------------------------------------------
