;+
; NAME:
;   spreduce
;
; PURPOSE:
;   Extract, wavelength-calibrate, and flatten SDSS spectral frame(s).
;
; CALLING SEQUENCE:
;   spreduce, flatname, arcname, objname, pixflatname=pixflatname, $
;    plugfile=plugfile, lampfile=lampfile, $
;    indir=indir, plugdir=plugdir, outdir=outdir, $
;    ecalibfile=ecalibfile, plottitle=plottitle, summarystruct=summarystruct
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
;   pixflatname- Name of pixel-to-pixel flat, produced with SPFLATTEN.
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in $IDLSPEC2D_DIR/etc.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   plugdir    - Input directory for PLUGFILE; default to '.'
;   outdir     - Directory for output files; default to '.'
;   ecalibfile - opECalib.par file for SDSSPROC
;   plottitle  - Prefix for titles in QA plots.
;   summarystruct - Good intentions ???
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
;   qaplot_arcline
;   qaplot_fflat
;   select_arc()
;   select_flat()
;   sortplugmap
;   sdssproc
;   spcalib
;   splog
;   yanny_free
;   yanny_read
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/skylines.dat
;
; REVISION HISTORY:
;   12-Oct-1999  Written by D. Schlegel & S. Burles, APO
;-
;------------------------------------------------------------------------------

pro spreduce, flatname, arcname, objname, pixflatname=pixflatname, $
 plugfile=plugfile, lampfile=lampfile, $
 indir=indir, plugdir=plugdir, outdir=outdir, $
 ecalibfile=ecalibfile, plottitle=plottitle, summarystruct=summarystruct

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = '.'

   stime0 = systime(1)

   ;---------------------------------------------------------------------------
   ; Locate skyline file for sky wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(skylinefile)) then begin
      skyfilenames = (findfile(skylinefile, count=ct))[0]
      if (ct EQ 0) then message, 'No SKYLINEFILE found '+skylinefile
   endif else begin
      skydefault = filepath('skylines.dat', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      skyfilenames = (findfile(skydefault, count=ct))[0]
      if (skyfilenames EQ '') then message, 'No SKYLINEFILE found '+skydefault
   endelse

   skylinefile = skyfilenames[0]

   ;---------------------------------------------------------------------------
   ; Determine spectrograph ID and color from first object file
   ;---------------------------------------------------------------------------

   sdssproc, objname[0], indir=indir, spectrographid=spectrographid, $
             color=color, ecalibfile=ecalibfile

   ;---------------------------------------------------------------------------
   ; Read PLUGMAP file and sort
   ;---------------------------------------------------------------------------
 
   plugpath = filepath(plugfile, root_dir=plugdir)
   plugfilename = (findfile(plugpath, count=ct))[0]
   if (ct NE 1) then $
    message, 'Cannot find plugMapFile ' + plugfile

   yanny_read, plugfilename, pstruct, hdr=hdrplug
   plugmap = *pstruct[0]
   yanny_free, pstruct

   ;-------------------------------------------------------------------------
   ; Plugsort will also return a mask of unplugged fibers
   ;-------------------------------------------------------------------------

   plugsort = sortplugmap(plugmap, spectrographid, fibermask=fibermask)
 
   ;---------------------------------------------------------------------------
   ; REDUCE CALIBRATION FRAMES
   ;---------------------------------------------------------------------------

   flatinfoname = filepath( $
       'spFlatInfo-'+string(format='(a1,i1,a)',color,spectrographid, $
       '-'), root_dir=outdir)

   arcinfoname = filepath( $
       'spArcInfo-'+string(format='(a1,i1,a)',color,spectrographid, $
       '-'), root_dir=outdir)

   spcalib, flatname, arcname, pixflatname=pixflatname, fibermask=fibermask, $
    lampfile=lampfile, indir=indir, arcstruct, flatstruct, $
    ecalibfile=ecalibfile, plottitle=plottitle, $
    flatinfoname=flatinfoname, arcinfoname=arcinfoname

   ;----------
   ; Select the best arc

   bestarc = select_arc(arcstruct)

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
   wset =  *(bestarc.wset)

   qaplot_arcline, *(bestarc.xdif_tset), wset, lambda, $
    color=color, title=plottitle+' Arcline Fit for '+bestarc.name

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
       pixflatname=pixflatname, spectrographid=spectrographid, color=color, $
       ecalibfile=ecalibfile, minflat = 0.2

      ;----------
      ; Construct the best flat for this object from all of the reduced
      ; flat-field frames

      ; Set TAI equal to the time half-way through the exposure
      tai = sxpar(objhdr, 'TAI') + 0.5 * sxpar(objhdr, 'EXPTIME')
      bestflat = select_flat(flatstruct, tai)
      splog, 'Best flat = ', bestflat.name

      if (NOT keyword_set(bestflat)) then begin

         splog, 'ABORT: No good flats (saturated?)'

      endif else begin

         xsol = *(bestflat.xsol)
         fflat = *(bestflat.fflat)
         widthset = *(bestflat.widthset)

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
          'spSpec2d-'+string(format='(a1,i1,a,i8.8,a)',color,spectrographid, $
          '-',framenum,'.fits'), root_dir=outdir)

         sxaddpar, objhdr, 'PLUGMAPF', plugfilename
         sxaddpar, objhdr, 'FLATFILE', bestflat.name
         sxaddpar, objhdr, 'ARCFILE', bestarc.name
         sxaddpar, objhdr, 'OBJFILE', objname[iobj]
         sxaddpar, objhdr, 'LAMPLIST', lampfile
         sxaddpar, objhdr, 'SKYLIST', skylinefile
         sxaddpar, objhdr, 'PIXFLAT', pixflatname

         ;-----
         ; Extract the object frame

         extract_object, outname, objhdr, image, invvar, plugsort, wset, $
          xpeak, lambda, xsol, fflat, fibermask, widthset=widthset, $
          color=color, plottitle=plottitle

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
