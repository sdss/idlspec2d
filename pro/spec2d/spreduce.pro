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
;                default to the file 'lamphgcdne.dat' (APO) or
;                'lampHeNe.dat' (LCO) in $IDLSPEC2D_DIR/etc.
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
;   MWM_fluxer  - Utilize MWM optional settings (ie gaia reddening and different S/N cuts)
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
;-     Nov-2018  Adapted for the SDSS-V BHM (H.Ibarra)
;      Apr-2019  Backwards for SDSS-IV data (H.Ibarra)
;      Oct-2021  Merging spreduce_legacy back in and updating to current fps standards (S.Morrison)
;-
;------------------------------------------------------------------------------
pro spreduce, flatname, arcname, objname, run2d=run2d, plugfile=plugfile, $
    calobjobssfile=calobjobssfile, lampfile=lampfile, indir=indir, plugdir=plugdir,$
    outdir=outdir, ecalibfile=ecalibfile, plottitle=plottitle, do_telluric=do_telluric,$
    writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, bbspec=bbspec, $
    splitsky=splitsky, nitersky=nitersky, plates=plates, legacy=legacy, saveraw=saveraw,$
    gaiaext=gaiaext, map3d = map3d, MWM_fluxer=MWM_fluxer, no_db=no_db,debug=debug,$
    corrline=corrline, nbundles=nbundles, bundlefibers=bundlefibers, $
    noreject=noreject, force_arc2trace=force_arc2trace

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = '.'
   if (NOT keyword_set(nitersky)) then nitersky = 1

   stime0 = systime(1)

   ;---------------------------------------------------------------------------
   ; Locate skyline file for sky wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(skylinefile)) then begin
      fullskyfile = (findfile(skylinefile, count=ct))[0]
   endif else begin
      skydefault = filepath('skylines.dat', root_dir=getenv('IDLSPEC2D_DIR'), $
                            subdirectory='etc')
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
   ; Read PLUGMAP/confSummary file and sort
   ; Look for photoPlate file in directory OUTDIR
   ;---------------------------------------------------------------------------
   if keyword_set(plates) or keyword_set(legacy) then begin
        plugmap = readplugmap(plugfile, spectrographid, plugdir=plugdir, $
                              cartid=sxpar(objhdr,'CARTID'),$
                              /calibobj, mjd=sxpar(objhdr,'MJD'),indir=outdir, $
                              exptime=sxpar(objhdr,'EXPTIME'), hdr=hdrplug, $
                              fibermask=fibermask, plates=plates, gaiaext=gaiaext,$
                              map3d = map3d, MWM_fluxer=MWM_fluxer, no_db=no_db)
        if (NOT keyword_set(plugmap)) then begin
            splog, 'ABORT: Plug map not found ' $
                + djs_filepath(plugfile, root_dir=plugdir)
            return
        endif
        plugmap = *(plugmap[0])
        fibermap = *(fibermap[0])
        hdrplug  = *(hdrplug[0])
    endif else begin
        fps=1
        calobssum = readplugmap(plugfile, spectrographid, plugdir=plugdir,$
                                /calibobj, mjd=sxpar(objhdr,'MJD'), indir=outdir, $
                                exptime=sxpar(objhdr,'EXPTIME'), hdr=hdrcal, $
                                fibermask=fibermaskcal, gaiaext=gaiaext, map3d = map3d,$
                                MWM_fluxer=MWM_fluxer, no_db=no_db)
        fibermaskcal = *(fibermaskcal[0])
        hdrcal = *(hdrcal[0])
        if (NOT keyword_set(calobssum)) then begin
            for i=0, n_elements(plugfile)-1 do begin
                splog, 'ABORT: obsSummary file not found ' $
                    + djs_filepath(plugfile[i], root_dir=plugdir)
            endfor
            return
        endif
    endelse
 
 
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

cart = strtrim(sxpar(objhdr, 'CARTID'),2)
if cart eq 'FPS-S' then obs = 'LCO' else obs = 'APO'
tai=sxpar(objhdr,'TAI-BEG')+(sxpar(objhdr, 'EXPTIME')/2.0)
airmass = tai2airmass(sxpar(objhdr,'RADEG'),sxpar(objhdr,'DECDEG'), tai=tai, site=obs)
   if keyword_set(fps) then begin
        cartid = strtrim(yanny_par_fc(hdrcal, 'cartridgeId'),2)
        spcalib, flatname, arcname, fibermask=fibermaskcal, cartid=cartid, $
                lampfile=lampfile, indir=indir, ecalibfile=ecalibfile, noreject=noreject, $
                plottitle=plottitle, flatinfoname=flatinfoname, arcinfoname=arcinfoname, $
                arcstruct=arcstruct, flatstruct=flatstruct, writeflatmodel=writeflatmodel, $
                writearcmodel=writearcmodel, bbspec=bbspec, plates=plates, legacy=legacy, $
                nbundles=nbundles, bundlefibers=bundlefibers, saveraw=saveraw, debug=debug, $
                traceflat=traceflat, force_arc2trace=force_arc2trace
   endif else begin
        cartid = strtrim(yanny_par_fc(hdrplug, 'cartridgeId'),2)

        spcalib, flatname, arcname, fibermask=fibermask, cartid=cartid, $
                lampfile=lampfile, indir=indir, ecalibfile=ecalibfile, noreject=noreject, $
                plottitle=plottitle, flatinfoname=flatinfoname, arcinfoname=arcinfoname, $
                arcstruct=arcstruct, flatstruct=flatstruct, writeflatmodel=writeflatmodel, $
                writearcmodel=writearcmodel, bbspec=bbspec, plates=plates, legacy=legacy, $
                nbundles=nbundles, bundlefibers=bundlefibers, saveraw=saveraw, debug=debug, $
                traceflat=traceflat

   endelse

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
   print,bestflat.iarc
;   bestarc = select_arc(arcstruct)
   bestarc = arcstruct[ bestflat.iarc ]

   if 1 eq 0 then begin
   ;if keyword_set(traceflat) then begin
        ; Use spTraceFlats for fiber throughputs and superflats instead of associated flats
        ; This was turned off and replaced with padding the ends of the fiberflats using the scaled
        ; traceflats
        iarc = bestflat.iarc
        ccd = strtrim(sxpar(*bestarc.hdr, 'CAMERAS'),2)
        tmjd = strtrim(sxpar(flathdr, 'MJD'),2)
        tflat = filepath('spTraceFlat-'+ccd+'-*.fits.gz', $
                         root_dir=get_trace_dir(tmjd))
        tflat = file_search(tflat, /fold_case, count=ct)
        if ct gt 0 then begin
            tflat=tflat[0]
            t_fflat = mrdfits(tflat, 0, t_hdr,/silent)
            t_fmask = mrdfits(tflat, 2, /silent)
            t_wset  = mrdfits(tflat, 3, /silent)
            t_sfset = mrdfits(tflat, 4, /silent)
            t_xsol  = mrdfits(tflat, 5, /silent)
            t_MEDWIDTH = [sxpar(t_hdr, 'MEDWIDT0'), sxpar(t_hdr, 'MEDWIDT1'), $
                          sxpar(t_hdr, 'MEDWIDT2'), sxpar(t_hdr, 'MEDWIDT3')]
            bestflat = create_struct( name='FLAT_STRUCT', $
                                      'NAME', FILE_BASENAME(tflat), $
                                      'IARC', iarc, $
                                      'PROFTYPE', sxpar(t_hdr, 'PROFTYPE'), $
                                      'MEDWIDTH', t_MEDWIDTH, $
                                      'FIBERMASK', ptr_new(t_fmask), $
                                      'XSOL', ptr_new(t_xsol), $
                                      'WIDTHSET', ptr_new(t_wset), $
                                      'FFLAT', ptr_new(t_fflat), $
                                      'SUPERFLATSET', ptr_new(t_sfset), $
                                      'HDR', ptr_new(t_hdr))
        endif
   endif

  foreach arc, arcstruct do begin
    if arc.qbad eq 1 then continue
    lambda = *(arc.lambda)
    xpeak = *(arc.xpeak)
    wset = *(arc.wset)
    dispset = *(arc.dispset)
    reslset = *(arc.reslset)

    qaplot_arcline, *(arc.xdif_tset), wset, lambda, $
        rejline=*(arc.rejline), $
        color=color, title=plottitle+' Arcline Fit for '+arc.name

   endforeach



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
   reslset = *(bestarc.reslset)

   qaplot_arcline, *(bestarc.xdif_tset), wset, lambda, $
    rejline=*(bestarc.rejline), $
    color=color, title=plottitle+' Best Arcline Fit for '+bestarc.name

   ; Generate the sdProc files for the arc images and the PSF models
   if (keyword_set(bbspec)) then begin
      sdssproc, bestarc.name, indir=indir, /outfile, /applycrosstalk
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
      bbspec_pixpsf, arcstr, flatstr, /batch, nbundles=nbundles, bundlefibers=bundlefibers
   endif

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------
   if keyword_set(fps) then begin
         objobssum = readplugmap(plugfile, spectrographid, plugdir=plugdir,$
                                /calibobj, mjd=sxpar(objhdr,'MJD'), indir=outdir, $
                                exptime=sxpar(objhdr,'EXPTIME'), hdr=hdrobj, $
                                fibermask=fibermaskobj, gaiaext=gaiaext, map3d = map3d,$
                                MWM_fluxer=MWM_fluxer)
        if (NOT keyword_set(objobssum)) then begin
            for i=0, n_elements(objobssfile)-1 do begin
                splog, 'ABORT: obsSummary file not found ' $
                    + djs_filepath(objobssfile[i], root_dir=plugdir)
            endfor
            return
        endif
   endif
   for iobj=0, N_elements(objname)-1 do begin

      stimeobj = systime(1)
      splog, prelog=objname[iobj]

      ;----------
      ; Read object image
      ;   Minflat will mask all pixels with low 2d pixelflat values
 
      splog, 'Reading object ', objname[iobj]
      sdssproc, objname[iobj], image, invvar, rdnoiseimg=rdnoise, $
            indir=indir, hdr=objhdr, spectrographid=spectrographid, color=color, $
            ecalibfile=ecalibfile, minflat=0.8, maxflat=1.2, outfile=saveraw, $
            nsatrow=nsatrow, fbadpix=fbadpix, /applycrosstalk, ccdmask=ccdmask

      ;-----
      ; Decide if this science exposure is bad
      if keyword_set(plates) then begin
            srvymode=sxpar(objhdr, 'SRVYMODE')
            if strmatch(srvymode, 'MWM lead', /fold_case) eq 1 then begin
                threshold=2000.
            endif else begin
                threshold=1000.
            endelse
      endif else begin
            racen = sxpar(objhdr, 'RACEN')
            deccen = sxpar(objhdr, 'DECCEN')
            euler, racen, deccen, ll, bb, 1
            if abs(bb) lt 15. then threshold=2000. else threshold=1000.
      endelse
      qbadsci = reject_science(image, objhdr, nsatrow=nsatrow, fbadpix=fbadpix, threshold=threshold)

      sxaddpar, objhdr, 'RUN2D', run2d, ' Spectro-2D reduction name'

      ; In case TAI-BEG,TAI-END were missing from the header, add them in.
      get_tai, objhdr, tai_beg, tai_mid, tai_end
      sxaddpar, objhdr, 'TAI-BEG', tai_beg
      sxaddpar, objhdr, 'TAI-END', tai_end

      sxaddpar, objhdr, 'FRAMESN2', 0.0
      if keyword_set(fps) then begin
            hdrplug= *(hdrobj[iobj])  ;hdrobj[nt[iobj]+1:nt[iobj+1]-1]
      endif else begin
            sxaddpar, objhdr, 'TILEID', long(yanny_par_fc(hdrplug, 'tileId')), $
                    'Tile ID for SDSS BOSS plates', after='PLATEID'
      endelse
      sxaddpar, objhdr, 'CARTID', strtrim(yanny_par_fc(hdrplug, 'cartridgeId'),2), $
                'Cartridge used in this plugging', after='PLATEID'
      redden = float(yanny_par_fc(hdrplug, 'reddeningMed'))
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
         reslset = *(bestarc.reslset)
         superflat_minval = bestflat.superflat_minval

         sxaddpar, objhdr, 'XSIGMA', max(bestflat.medwidth)
         sxaddpar, objhdr, 'WSIGMA', max(bestarc.medwidth)
         sxaddpar, objhdr, 'WDISPR', max(bestarc.medresol)

         qaplot_fflat, fflat, wset, $
          title=plottitle+'Fiber-Flats for '+bestflat.name

         ;----------
         ; Combine FIBERMASK bits from the plug-map file, best flat
         ; and best arc
         if keyword_set(fps) then begin
                plugmap   = *(objobssum[iobj])
                fibermask = *(fibermaskobj[iobj])
         endif else begin
                superflat_minval = 0
         endelse
         fibermask = fibermask $
                    OR (*(bestflat.fibermask) AND fibermask_bits('BADTRACE')) $
                    OR (*(bestflat.fibermask) AND fibermask_bits('BADFLAT')) $
                    OR (*(bestarc.fibermask) AND fibermask_bits('BADARC'))
         plugmap.fibermask = fibermask
         ;----------
         ; Determine output file name and modify the object header

         framenum = sxpar(objhdr, 'EXPOSURE')
         outname = filepath( $
          'spFrame-'+string(format='(a1,i1,a,i8.8,a)',color,spectrographid, $
          '-',framenum,'.fits'), root_dir=outdir)

         if keyword_set(fps) then $
            sxaddpar, objhdr, 'confSFILE', fileandpath(plugfile[iobj]) $
            else sxaddpar, objhdr, 'PLUGFILE', fileandpath(plugfile)
         if keyword_set(traceflat) then $
            sxaddpar, objhdr, 'TRACFLAT', fileandpath(traceflat)
         sxaddpar, objhdr, 'FLATFILE', fileandpath(bestflat.name)
         sxaddpar, objhdr, 'ARCFILE', fileandpath(bestarc.name)
         sxaddpar, objhdr, 'OBJFILE', fileandpath(objname[iobj])
         sxaddpar, objhdr, 'LAMPLIST', fileandpath(lampfile)
         sxaddpar, objhdr, 'SKYLIST', fileandpath(fullskyfile)
         sxaddpar, objhdr, 'OBSMODE', strtrim(yanny_par_fc(hdrplug, 'OBSMODE'),2)
         ;-----
         ; Extract the object frame
        
        extract_object, outname, objhdr, image, invvar, rdnoise, plugmap, wset, $
                xpeak, lambda, xsol, fflat, fibermask, color=color, $
                proftype=proftype, superflatset=superflatset, reslset=reslset, $
                widthset=widthset, dispset=dispset, skylinefile=fullskyfile, $
                plottitle=plottitle, do_telluric=do_telluric, bbspec=bbspec, $
                splitsky=splitsky, ccdmask=ccdmask, nitersky=nitersky, corrline=corrline, $
                nbundles=nbundles, bundlefibers=bundlefibers, debug=debug, $
                superflat_minval = superflat_minval

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
