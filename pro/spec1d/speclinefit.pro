;+
; NAME:
;   speclinefit
;
; PURPOSE:
;   Line-fitting calling script for entire plate(s)
;
; CALLING SEQUENCE:
;   speclinefit, [ platefile, fiberid=, $
;    hdr=, objflux=, objivar=, $
;    zhdr=, zans=, synflux=, dispflux=, $
;    outfile=, zline=, yfit=, /doplot, /debug ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platefile  - Plate file(s) from spectro-2D; default to all files
;                matching 'spPlate*.fits' if neither PLATEFILE nor OBJFLUX
;                are set.
;                If specified, then all inputs are read from this file
;                and its corresponding spZbest file.  The following inputs
;                are ignored: HDR,OBJFLUX,OBJIVAR,ZHDR,ZANS,SYNFLUX,DISPFLUX.
;                If set, then an output FITS file and a log file are created,
;                as is a PostScript plot file if /DOPLOT is set.
;   fiberid    - If specified, then only reduce these fiber numbers;
;                this must be a vector with unique values between 1 and
;                the number of rows in the plate file (typically 640).
;                This keyword must be set to only those fibers that exist
;                in the spZbest file if FIBERID was specified when
;                running SPREDUCE1D.
;   hdr        - Header from plate file.  The fields COEFF0,COEFF1 are
;                used for the wavelength mapping.
;   objflux    - Object spectra [NPIX,NOBJ]
;   objivar    - Object inverse variance [NPIX,NOBJ]
;   zhdr       - Optional header for output file.
;   zans       - Output structure from SPREDUCE1D.  The following fields
;                are used in this procedure: PLATE,MJD,FIBERID,CLASS,Z.
;   synflux    - Best-fit synthetic spectra; used for STAR background terms.
;   dispflux   - Best-fit dispersion spectra; used for GALAXY background terms.
;   outfile    - Name of output file (the name is generated from PLATEFILE
;                if that is set).
;   doplot     - If set, then generate plots.  Send plots to a PostScript
;                file if PLATEFILE is set and /DEBUG is not set.
;   debug      - If set, then send plots to the X display and wait for
;                a keystroke after each plot; setting /DEBUG forces /DOPLOT.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   zline      - Output structure with line fits [NLINE,NOBJ].
;                Entries are blank where nothing was measured.
;   yfit       - Best-fit background terms + line fits [NPIX,NOBJ]
;
; COMMENTS:
;   If PLATEFILE is set, then names of other input/output files are derived
;   from that name.  For example, if PLATEFILE='spPlate-0306-51690.fits', then
;     ZBESTFILE = 'spZbest-0306-51690.fits'   (input file)
;     ZLINEFILE = 'spZline-0306-51690.fits'   (input file)
;     LOGFILE   = 'spDiagLine-0306-51690.log' (output file)
;     PLOTFILE  = 'spDiagLine-0306-51690.ps'  (output file)
;   An output file is written only if PLATEFILE or OUTFILE are set.
;
;   The first background term is simply the best-fit dispersion template
;   for galaxies, nothing for QSOs, or the synthetic spectrum for other
;   types of objects.
;   At any wavelengths where this background model does not exist (at the
;   line center), then fit another term that is a linear term
;   with a width of 0.030 in log-wavelength (300 pix) for QSOs, or a
;   linear term with a width of 0.005 (50 pix) for other types of objects.
;
;   The initial guess for the dispersion is 2000 km/sec for QSO's
;   (ZGUESS=0.0030), and 105 kms/sec for all other objects (ZGUESS=0.00015).
;
;   The continuum level is evaluated at the line center for QSOs,
;   and otherwise is set to [300,600] km/s on either side of the line center.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/emlines.par
;
; PROCEDURES CALLED:
;   cpbackup
;   create_linestruct()
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   fileandpath()
;   headfits()
;   linebackfit()
;   mrdfits()
;   mwrfits
;   poly_array()
;   splog
;   skymask()
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   12-Feb-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro speclinefit, platefile, fiberid=fiberid, $
 hdr=hdr, objflux=objflux, objivar=objivar, $
 zhdr=zhdr, zans=zans, synflux=synflux, dispflux=dispflux, $
 outfile=outfile, zline=zline, yfit=yfit, doplot=doplot, debug=debug

   if (NOT keyword_set(platefile) AND NOT keyword_set(objflux)) then begin
      platefile = findfile('spPlate*.fits*', count=nfile)
   endif else begin
      if (size(platefile,/tname) EQ 'STRING') then $
       nfile = n_elements(platefile) $
      else $
       nfile = 0
   endelse
   if (keyword_set(debug)) then doplot = 1

   ;----------
   ; If multiple plate files exist, then call this script recursively
   ; for each such plate file.

   if (nfile GT 1) then begin
      for i=0, nfile-1 do begin
         speclinefit, platefile[i], fiberid=fiberid, $
          doplot=doplot, debug=debug
      endfor
      return
   endif else if (nfile EQ 1) then begin
      platefile = platefile[0] ; in case this is a 1-element array
   endif

   ;----------
   ; Determine names of output files

   if (keyword_set(platefile)) then begin
      thisname = fileandpath(platefile, path=thispath)

      zbestfile = djs_filepath(repstr(thisname,'spPlate','spZbest'), $
       root_dir=thispath)

      outfile = djs_filepath(repstr(thisname,'spPlate','spZline'), $
       root_dir=thispath)

      logfile = repstr(thisname,'spPlate','spDiagLine')
      logfile = djs_filepath(repstr(logfile,'fits','log'), root_dir=thispath)

      if (keyword_set(doplot) AND NOT keyword_set(debug)) then $
       plotfile = djs_filepath(repstr(logfile,'log','ps'), root_dir=thispath)

      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + systime()
      if (keyword_set(plotfile)) then begin
         splog, 'Plot file ' + plotfile
         cpbackup, plotfile
         dfpsplot, plotfile, /color
      endif
      splog, 'IDL version: ' + string(!version,format='(99(a," "))')
      spawn, 'uname -a', uname
      splog, 'UNAME: ' + uname[0]
      splog, 'idlspec2d version ' + idlspec2d_version()
      splog, 'idlutils version ' + idlutils_version()
   endif

   stime0 = systime(1)

   ;----------
   ; Read the 2D output file

   if (keyword_set(platefile)) then begin
      splog, 'Reading spPlate file ', platefile

      objflux = mrdfits(platefile,0,hdr)
      if (NOT keyword_set(hdr)) then begin
         splog, 'ABORT: Plate file not valid: ' + platefile, /close
         return
      endif
      objivar = mrdfits(platefile,1)
      andmask = mrdfits(platefile,2)
      ormask = mrdfits(platefile,3)

      objivar = skymask(objivar, andmask, ormask)
      andmask = 0 ; Free memory
      ormask = 0 ; Free memory

      splog, 'Reading spZbest file ', zbestfile

      zhdr = headfits(zbestfile)
      zans = mrdfits(zbestfile, 1)
      synflux = mrdfits(zbestfile, 2)
      dispflux = mrdfits(zbestfile, 3)
   endif

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]

   if (NOT keyword_set(hdr) $
    OR NOT keyword_set(objflux) $
    OR NOT keyword_set(objivar) $
    OR NOT keyword_set(zans) $
    OR NOT keyword_set(synflux) $
    OR NOT keyword_set(dispflux)) then begin
      splog, 'ABORT: Unable to read all input arrays', $
       close=keyword_set(platefile)
      return
   endif

   ;----------
   ; Construct the wavelength mapping

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')
   objloglam = objloglam0 + lindgen(npixobj) * objdloglam

   ;----------
   ; Trim to specified fibers if FIBERID is set.
   ; The following logic works if the spZbest file does not contain
   ; all fibers, for example if SPREDUCE1D was run with FIBERID set.

   if (keyword_set(fiberid)) then begin
      fiblist = fiberid
   endif else begin
      fiblist = lindgen(nobj) + 1
   endelse

   ; Find the object in the ZANS structure that matches each fiber ID.
   nfib = n_elements(fiblist)
   zindx = lonarr(nfib)
   for i=0, nfib-1 do zindx[i] = (where(zans.fiberid EQ fiblist[i]))[0]
   if ((where(zindx EQ -1))[0] NE -1) then $
    message, 'Some FIBERIDs do not exist in spZbest file'

   ; Now trim the arrays and structure to the specified fibers
   if (keyword_set(platefile)) then begin
      if (min(fiblist) LT 1 OR max(fiblist) GT nobj) then $
       message, 'Invalid value for FIBERID: must be between 1 and '+string(nobj)
      objflux = objflux[*,fiblist-1]
      objivar = objivar[*,fiblist-1]
   endif
   zans = zans[zindx]
   synflux = synflux[*,zindx]
   dispflux = dispflux[*,zindx]
   nobj = nfib
   splog, 'Number of fibers = ', nobj

   ;----------
   ; Read line lists and convert to vacuum

   linefile = filepath('emlines.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   yanny_read, linefile, pdata
   linelist = *pdata[0]
   yanny_free, pdata

   vaclambda = linelist.lambda
   airtovac, vaclambda
   linelist.lambda = vaclambda

   ;----------
   ; Create the output data structure and array

   lfitall = create_linestruct()
   lfitall = replicate(lfitall, n_elements(linelist), nobj)

   if (arg_present(yfit) OR keyword_set(outfile)) then $
    yfit = make_array(size=size(objflux), /float)

   ;----------
   ; Generate the additional tags to add to the output structure.
   ; (The values below are copied from the ZANS structure).

   res1 = { plate:    0L, $
            mjd:      0L, $
            fiberid:  0L        }
   res_prepend = make_array(value=res1, dimension=size(lfitall,/dimens))

   ;----------
   ; Loop through each object and do the line-fitting

   fiberlist = 0

   for iobj=0, nobj-1 do begin
      splog, 'Fitting object #', iobj

      ; If for any weird reason PLATE,MJD in the ZANS struct are different
      ; from those values in the header, use the values from the ZANS structure.
      res_prepend[*,iobj].plate = zans[iobj].plate
      res_prepend[*,iobj].mjd = zans[iobj].mjd
      res_prepend[*,iobj].fiberid = zans[iobj].fiberid

      if (strtrim(zans[iobj].class,2) EQ 'QSO') then begin
         thiswidth = 0.030
         npoly = 2
      endif else begin
         thiswidth = 0.005
         npoly = 2
      endelse
      polyback = poly_array(npixobj,npoly)

      t0 = systime(1)
      zguess = zans[iobj].z

      ; The first background term is simply the best-fit dispersion template
      ; for galaxies, or the synthetic spectrum for other types of objects.
      if (strtrim(zans[iobj].class,2) EQ 'GALAXY') then begin
         background = dispflux[*,iobj]
         contrange = [300.,600.]
         sigguess = 1.5d-4
      endif else if (strtrim(zans[iobj].class,2) EQ 'QSO') then begin
         background = 0
         contrange = 0
         sigguess = 0.003d0
      endif else begin
         background = synflux[*,iobj]
         contrange = [300.,600.]
         sigguess = 1.5d-4
      endelse

      if (keyword_set(background)) then fitmask = background NE 0 $
       else fitmask = bytarr(npixobj)
      ipix = where(fitmask)
      if (ipix[0] NE -1) then begin
         minwave = 10^min(objloglam[ipix])
         maxwave = 10^max(objloglam[ipix])
      endif else begin
         minwave = 0
         maxwave = 0
      endelse

      ; Construct the flat background terms around any line that
      ; is not centered within the synthetic or model template.
      nflat = 0
      for iline=0, n_elements(linelist)-1 do begin
         obswave = linelist[iline].lambda * (1 + zguess)
         if (obswave LE minwave OR obswave GE maxwave) then begin
            thismask = objloglam GT alog10(obswave) - 0.5*thiswidth $
             AND objloglam LT alog10(obswave) + 0.5*thiswidth
            junk = where(thismask, nthis)
            if (nthis GT 1) then begin
               if (nflat EQ 0) then flatback = long(thismask) $
                else flatback = [[flatback], [long(thismask)]]
               nflat = nflat + 1
            endif
         endif
      endfor

      ; Make certain that the simple background terms never overlap
      ; the dispersion template
      if (nflat GT 0 AND keyword_set(background)) then begin
         if (nflat EQ 1) then allflat = flatback NE 0 $
          else allflat = total(flatback,2) NE 0
         background = background * (allflat EQ 0)
      endif

      ; Combine the simple background terms that overlap one another
      if (nflat GT 1) then begin
         ; First smush all these terms together
         allflat = long(total(flatback,2) NE 0)

         ; Now split them back up again, and append to the "background".
         ; Normalize the level to the median level in the object flux
         ; in that wavelength region.
         allflat = [0, allflat, 0]
         dflat = allflat[1:npixobj+1] - allflat[0:npixobj]
         istart = where(dflat EQ 1, nflat)
         iend = where(dflat EQ -1, nend) - 1

         for iflat=0, nflat-1 do begin
            medval = djs_median( objflux[istart[iflat]:iend[iflat],iobj] )
            medval = medval + (medval EQ 0) ; Force to non-zero value
            thisterm = fltarr(npixobj,npoly)
            thisterm[istart[iflat]:iend[iflat],*] = medval
            thisterm = thisterm * polyback
            if (NOT keyword_set(background)) then background = thisterm $
             else background = [[background], [thisterm]]
            fitmask[istart[iflat]:iend[iflat],0] = 1
         endfor
      endif

      ; Fill the output structures with default values
      lfitall[*,iobj].linez_err = -1
      lfitall[*,iobj].linesigma_err = -1
      lfitall[*,iobj].linearea_err = -1
      lfitall[*,iobj].linecontlevel_err = -1
      lfitall[*,iobj].linechi2 = -1
      lfitall[*,iobj].linename = linelist.name
      lfitall[*,iobj].linewave = linelist.lambda

      ; Call the line-fitting engine.
      ; Restrict the fit to only those lines within the wavelength region
      ; at the guessed redshift (+/- 6000 km/sec).
      ipix = where(fitmask)
      if (ipix[0] NE -1) then begin
         minwave = 10^min(objloglam[ipix])
         maxwave = 10^max(objloglam[ipix])
      endif else begin
         minwave = 0
         maxwave = 0
      endelse
      iuse = where(linelist.lambda GE minwave/1.02/(1+zguess) $
       AND linelist.lambda LE maxwave*1.02/(1+zguess), nuse)
      if (nuse GT 0) then begin
         lfit1 = linebackfit(linelist[iuse].lambda, objloglam, $
          objflux[*,iobj], invvar=objivar[*,iobj] * fitmask, $
          linename=linelist[iuse].name, $
          zindex=linelist[iuse].zindex, windex=linelist[iuse].windex, $
          findex=linelist[iuse].findex, fvalue=linelist[iuse].fvalue, $
          zguess=zguess, sigguess=sigguess, $
          background=background, contrange=contrange, $
          yfit=yfit1, bfit=bfit, bterms=bterms)

         if (keyword_set(yfit)) then yfit[*,iobj] = yfit1

         ; Fill the output structures
         lfitall[iuse,iobj] = lfit1
      endif

      splog, 'Object #', iobj, ' CPU time for line fitting = ', $
       systime(1)-t0

      if (keyword_set(doplot)) then begin
         igood = where(objivar[*,iobj] GT 0, ngood)
         if (ngood GT 1) then begin
; ???
;stop
;set_plot,'x'
;splot, 10^objloglam[igood], objflux[igood,iobj]
;soplot, 10^objloglam[igood], background[igood,0], color='blue'
;soplot, 10^objloglam[igood], yfit1[igood], color='red'
            plottitle = string(zans[iobj].plate, zans[iobj].mjd, $
             zans[iobj].fiberid, $
             format='("Plate=",i4," MJD=",i5," Fiber=", i3)')
            djs_plot, 10^objloglam[igood], objflux[igood,iobj], $
             xtitle='Wavelength [Ang]', ytitle='Flux', title=plottitle
            djs_oplot, 10^objloglam[igood], yfit1[igood], color='red'

            ; Wait for a keystroke...
            if (keyword_set(debug)) then begin
               print, 'Press any key...'
               cc = strupcase(get_kbrd(1))
            endif
         endif
      endif
   endfor

   ;----------
   ; Concatenate the two output structures into a single structure

   lfitall = struct_addtags(res_prepend, lfitall)

   ;----------
   ; Write the output file

   outhdr = zhdr
   sxaddpar, outhdr, 'DIMS0', n_elements(linelist), ' Number of emission lines'
   sxaddpar, outhdr, 'DIMS1', nobj, ' Number of objects'
   sxaddpar, outhdr, 'VERSLINE', idlspec2d_version(), $
    'Version of idlspec2d for line fitting', after='VERS1D'

   if (keyword_set(outfile)) then begin
      mwrfits, 0, outfile, outhdr, /create ; Retain original header in first HDU
      mwrfits, lfitall, outfile
      mwrfits, yfit, outfile
   endif

   ;----------
   ; Close log file

   if (keyword_set(plotfile)) then dfpsclose

   splog, 'Total time for SPECLINEFIT = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'

   if (keyword_set(platefile)) then $
    splog, 'Successful completion of SPECLINEFIT at ', systime(), /close

   zline = temporary(lfitall)
   return
end
;------------------------------------------------------------------------------
