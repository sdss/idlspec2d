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
;    indir=indir, plugdir=plugdir, outdir=outdir
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
;                default to the file 'lamphgcdne.dat' in the IDL path.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   plugdir    - Input directory for PLUGFILE; default to '.'
;   outdir     - Directory for output files; default to '.'
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
;   djs_locate_file()
;   djs_median()
;   djs_plot
;   extract_image
;   extract_object
;   fiberflat()
;   fitarcimage
;   qaplot_arcline
;   qaplot_fflat
;   sortplugmap
;   sdssproc
;   splog
;   trace320crude
;   traceset2xy
;   xy2traceset
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   12-Oct-1999  Written by D. Schlegel & S. Burles, APO
;-
;------------------------------------------------------------------------------

pro spreduce, flatname, arcname, objname, pixflatname=pixflatname, $
 plugfile=plugfile, lampfile=lampfile, $
 indir=indir, plugdir=plugdir, outdir=outdir

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
      skydefault = getenv('IDLSPEC2D_DIR') + '/etc/skylines.dat'
      skyfilenames = (findfile(skydefault, count=ct))[0]
      if (skyfilenames EQ '') then message, 'No SKYLINEFILE found '+skydefault
   endelse

   skylinefile = skyfilenames[0]

   ;---------------------------------------------------------------------------
   ; Determine spectrograph ID and color from first object file
   ;---------------------------------------------------------------------------

   sdssproc, objname[0], indir=indir, spectrographid=spectrographid, color=color

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
   ; Plugsort will return mask of good (1) and bad (0) fibers too
   ;-------------------------------------------------------------------------

   plugsort = sortplugmap(plugmap, spectrographid, fibermask)
 
   ;---------------------------------------------------------------------------
   ; LOOP THROUGH FLAT+ARC IMAGE TO IDENTIFY THE BEST PAIR
   ;---------------------------------------------------------------------------

   nflat = N_elements(flatname)
   narc = N_elements(arcname)

   ibest = -1 ; Index number for best flat+arc pair
   bestcorr = -2.0 ; Any call to FITARCIMAGE returns values [-1,1]
   oldflatfile = ''

   for ifile=0, (nflat<narc)-1 do begin

      splog, ifile+1, (nflat<narc), $
       format='("Looping through flat+arc pair #",I3," of",I3)'

      ; Check if this is the same flat field as the last one read

      if (flatname[ifile] NE oldflatfile) then begin

         ;---------------------------------------------------------------------
         ; Read flat-field image
         ;---------------------------------------------------------------------

         splog, 'Reading flat ', flatname[ifile]
         sdssproc, flatname[ifile], image, invvar, indir=indir, $
          hdr=flathdr, pixflatname=pixflatname, nsatrow=nsatrow

         ;-----
         ; Decide if this flat is bad:
         ;   Reject if more than 1% of the pixels are marked as bad.
         ;   Reject if more than 10 rows are saturated.

         fbadpix = float(N_elements(where(invvar EQ 0))) / N_elements(invvar)

         qbadflat = 0
         if (fbadpix GT 0.01) then begin
            qbadflat = 1
            splog, 'Reject flat ' + flatname[ifile] + $
             ' (' + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels)'
         endif
         if (nsatrow GT 10) then begin
            qbadflat = 1
            splog, 'Reject flat ' + flatname[ifile] + $
             ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
         endif

         if (NOT qbadflat) then begin
            ;------------------------------------------------------------------
            ; Create spatial tracing from flat-field image
            ;------------------------------------------------------------------

            splog, 'Tracing 320 fibers in ',  flatname[ifile]
            tmp_xsol = trace320crude(image, invvar, yset=ycen, maxdev=0.15)

            splog, 'Fitting traces in ',  flatname[ifile]
            xy2traceset, ycen, tmp_xsol, tset, ncoeff=5, maxdev=0.1
            traceset2xy, tset, ycen, tmp_xsol
         endif

         oldflatfile = flatname[ifile]
      endif

      if (NOT qbadflat) then begin

         ;---------------------------------------------------------------------
         ; Read the arc
         ;---------------------------------------------------------------------

         splog, 'Reading arc ', arcname[ifile]
         sdssproc, arcname[ifile], image, invvar, indir=indir, $
          hdr=archdr, pixflatname=pixflatname, nsatrow=nsatrow

         ;-----
         ; Decide if this arc is bad:
         ;   Reject if more than 40 rows are saturated.

         qbadarc = 0
         if (nsatrow GT 40) then begin
            qbadarc = 1
            splog, 'Reject arc ' + arcname[ifile] + $
             ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
         endif

         if (NOT qbadarc) then begin
            ;------------------------------------------------------------------
            ; Extract the arc image
            ;------------------------------------------------------------------

            splog, 'Extracting arc image with simple gaussian'
            sigma = 1.0
            proftype = 1 ; Gaussian
            highrej = 15
            lowrej = 15
            npoly = 1 ; maybe more structure
            wfixed = [1,1] ; Just fit the first gaussian term

            extract_image, image, invvar, tmp_xsol, sigma, flux, fluxivar, $
             proftype=proftype, wfixed=wfixed, $
             highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

            ;-------------------------------------------------------------------
            ; Compute correlation coefficient for this arc image
            ;-------------------------------------------------------------------

            splog, 'Searching for wavelength solution'
            tmp_aset = 0
            fitarcimage, flux, fluxivar, aset=tmp_aset, $
             color=color, lampfile=lampfile, bestcorr=corr
; FOR NOW, REVERT TO THE OLD CODE! ???
;            fitarcimage_old, flux, fluxivar, aset=tmp_aset, $
;             color=color, lampfile=lampfile, bestcorr=corr

            ;-----
            ; Determine if this is the best flat+arc pair
            ; If so, then save the information that we need

            if (corr GT bestcorr) then begin
               ibest = ifile
               bestcorr = corr
               arcimg = flux
               arcivar = fluxivar
               xsol = tmp_xsol
               aset = tmp_aset
            endif
         endif
      endif

   endfor

   ;---------------------------------------------------------------------------
   ; Make sure that the best flat+arc pair is good enough
   ;---------------------------------------------------------------------------

   if (ibest EQ -1) then begin
      splog, 'ABORT: No good flats and/or no good arcs (saturated?)'
      return
   endif

   splog, 'Best flat = ', flatname[ibest]
   splog, 'Best arc = ', arcname[ibest]

   if ((color EQ 'blue' AND bestcorr LT 0.5) $
    OR (color EQ 'red'  AND bestcorr LT 0.5) ) then begin
      splog, 'ABORT: Best arc correlation = ', bestcorr
      return
   endif else $
    if ((color EQ 'blue' AND bestcorr LT 0.7) $
    OR (color EQ 'red'  AND bestcorr LT 0.7) ) then $
      splog, 'WARNING: Best arc correlation = ', bestcorr

   ;---------------------------------------------------------------------------
   ; Compute wavelength calibration for arc lamp only
   ;---------------------------------------------------------------------------

   arccoeff = 6

   splog, 'Searching for wavelength solution'
   fitarcimage, arcimg, arcivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
    aset=aset, fibermask=fibermask, wfirst=wfirst, $
    color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset
; FOR NOW, REVERT TO THE OLD CODE! ???
;   fitarcimage_old, arcimg, arcivar, xpeak, ypeak, wsetold, ncoeff=arccoeff, $
;    aset=aset, $
;    color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset

   if (NOT keyword_set(wset)) then begin
      splog, 'ABORT: Wavelength solution failed'
      return 
   endif

   qaplot_arcline, xdif_tset, lambda, filename=arcname[ibest], color=color

   ;---------------------------------------------------------------------
   ; Read best flat-field image (again)
   ;---------------------------------------------------------------------

   splog, 'Reading flat ', flatname[ibest]
   sdssproc, flatname[ibest], image, invvar, indir=indir, $
    hdr=flathdr, pixflatname=pixflatname

   ;---------------------------------------------------------------------------
   ; Extract the flat-field image
   ;---------------------------------------------------------------------------

   splog, 'Extracting flat-field image with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 15
   lowrej = 15
   npoly = 1  ; just fit flat background to each row
   wfixed = [1,1] ; Just fit the first gaussian term

      ; --->Kill<--- or adjust first and last column ???
      invvar[0,*] = 0.0
      invvar[2047,*] = 0.0

   extract_image, image, invvar, xsol, sigma, flat_flux, flat_fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

   junk = where(flat_flux GT 1.0e5, nbright)
   splog, 'Found ', nbright, ' bright pixels in extracted flat ', $
    flatname[ibest], format='(a,i7,a,a)'

   ;---------------------------------------------------------------------------
   ; Compute fiber-to-fiber flat-field variations
   ;---------------------------------------------------------------------------

;   fflat = fiberflat(flat_flux, flat_fluxivar, wset, fibermask=fibermask)
   fflat = fiberflat(flat_flux, flat_fluxivar, wset, fibermask=fibermask, $
    /dospline)

   qaplot_fflat, fflat, wset, filename=flatname[ibest]

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------

   for iobj=0, N_elements(objname)-1 do begin

      stimeobj = systime(1)
      splog, camname=objname[iobj]

      ;------------------
      ; Read object image

      splog, 'Reading object ', objname[iobj]
      sdssproc, objname[iobj], image, invvar, indir=indir, hdr=objhdr, $
       pixflatname=pixflatname, spectrographid=spectrographid, color=color

      framenum = sxpar(objhdr, 'EXPOSURE')
      outname = filepath( $
       'spSpec2d-'+string(format='(a1,i1,a,i8.8,a)',color,spectrographid, $
       '-',framenum,'.fits'), root_dir=outdir)

      sxaddpar, objhdr, 'PLUGMAPF', plugfilename
      sxaddpar, objhdr, 'FLATFILE', flatname[ibest]
      sxaddpar, objhdr, 'ARCFILE', arcname[ibest]
      sxaddpar, objhdr, 'OBJFILE', objname[iobj]
      sxaddpar, objhdr, 'LAMPLIST', lampfile
      sxaddpar, objhdr, 'SKYLIST', skylinefile
      sxaddpar, objhdr, 'PIXFLAT', pixflatname

      extract_object, outname, objhdr, image, invvar, plugsort, wset, $
               xpeak, lambda, xsol, fflat, fibermask, color=color

      splog, 'Elapsed time = ', systime(1)-stimeobj, ' seconds', $
       format='(a,f6.0,a)' 

   endfor

   heap_gc   ; Garbage collection for all lost pointers
   splog, 'Elapsed time = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)', camname=''

   return
end
;------------------------------------------------------------------------------
