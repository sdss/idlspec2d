;+
; NAME:
;   speclinefit
;
; PURPOSE:
;   Line-fitting calling script for entire plate(s)
;
; CALLING SEQUENCE:
;   speclinefit, [ platefile, fiberid=, /doplot, /debug ]
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
;                file unless /DEBUG is set.
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
;     ZBESTFILE = 'spZbest-0306-51690.fits'
;     ZLINEFILE = 'spZline-0306-51690.fits'
;
;   The output structure is always dimensioned [NLINE,NOBJ], with blank entries
;   where nothing was measured.
;
; EXAMPLES:
;
; BUGS:
;   At present, only do fits for galaxies, not stars or quasars. ???
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/emlines.par
;
; PROCEDURES CALLED:
;   cpbackup
;   create_linestruct()
;   dfpsclose
;   dfpsplot
;   headfits()
;   linebackfit()
;   mrdfits()
;   mwrfits
;   splog
;   skymask()
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   12-Feb-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro speclinefit, platefile, fiberid=fiberid, doplot=doplot, debug=debug

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
         speclinefit, platefile[i], fiberid=fiberid, doplot=doplot, debug=debug
      endfor
      return
   endelse

   ;----------
   ; Set the half-width of the wavelength region to use for background level
   ; fitting in log10(Ang); default is 10.e-4 = 10 pix. ???

   backwidth = 10.e-4

   ;----------
   ; Determine names of output files

   platemjd = strmid(platefile, 8, 10)

   zbestfile = 'spZbest-' + platemjd + '.fits'
   zlinefile = 'spZline-' + platemjd + '.fits'
   if (NOT keyword_set(logfile)) then $
    logfile = 'spDiagLine-' + platemjd + '.log'

   if (keyword_set(doplot) AND NOT keyword_set(debug)) then begin
      plotfile = 'spDiagLine-' + platemjd + '.ps'
      cpbackup, plotfile
      dfpsplot, plotfile, /color
   endif

   stime0 = systime(1)

   cpbackup, logfile
   splog, filename=logfile
   splog, 'Log file ' + logfile + ' opened ' + systime()
   if (keyword_set(plotfile)) then $
    splog, 'Plot file ' + plotfile
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   if (NOT keyword_set(hdr)) then begin
      splog, 'ABORT: Plate file not valid: ' + platefile, /close
      return
   endif
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
   objloglam = objloglam0 + lindgen(npixobj) * objdloglam

   ;----------
   ; Read the P-1D output file

   zhdr = headfits(zbestfile)
   zans = mrdfits(zbestfile, 1)
   dispflux = mrdfits(zbestfile, 3)

   if (NOT keyword_set(objflux) $
    OR NOT keyword_set(zans) OR NOT keyword_set(dispflux)) then begin
      splog, 'ABORT: Unable to read all input files', /close
      return
   endif

   ;----------
   ; Trim to specified fibers if FIBERID is set

   if (keyword_set(fiberid)) then begin
      if (min(fiberid) LT 0 OR max(fiberid) GT nobj) then $
       message, 'Invalid value for FIBERID: must be between 0 and '+string(nobj)
      objflux = objflux[*,fiberid-1]
      objivar = objivar[*,fiberid-1]
      anyandmask = anyandmask[fiberid-1]
      anyormask = anyormask[fiberid-1]
      plugmap = plugmap[fiberid-1]
      zans = zans[fiberid-1]
      nobj = n_elements(fiberid)
   endif else begin
      fiberid = lindgen(nobj) + 1
   endelse
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
   ; Create the output structure

   lfitall = create_linestruct()
   lfitall = replicate(lfitall, n_elements(linelist), nobj)

   ;----------
   ; Generate the additional tags to add to the output structure.
   ; (The default values below are actually over-written by the values in ZANS).

   res1 = { plate:    long(sxpar(zhdr, 'PLATEID')), $
            tile:     long(sxpar(zhdr, 'TILEID')), $
            mjd:      long(sxpar(zhdr, 'MJD')), $
            fiberid:  0L        }
   res_prepend = make_array(value=res1, dimension=size(lfitall,/dimens))

   ;----------
   ; Loop through each object and do the line-fitting

   fiberlist = 0

   for iobj=0, nobj-1 do begin
      ; If for any weird reason PLATE,TILE,MJD are different in the ZANS structure
      ; than those values in the header, use the values from the ZANS structure.
      res_prepend[*,iobj].plate = zans[iobj].plate
      res_prepend[*,iobj].tile = zans[iobj].tile
      res_prepend[*,iobj].mjd = zans[iobj].mjd
      res_prepend[*,iobj].fiberid = zans[iobj].fiberid

      if (strtrim(zans[iobj].class,2) EQ 'GALAXY') then begin
         splog, 'Fitting object #', iobj+1

         t0 = systime(1)
         zguess = zans[iobj].z

         ; The first background term is simply the best-fit dispersion template
         background = dispflux[*,iobj]
         ipix = where(background NE 0) > 0
         minwave = 10^min(objloglam[ipix])
         maxwave = 10^max(objloglam[ipix])

         ; Construct the flat background terms around any line that
         ; is not centered within the dispersion template.
         nflat = 0
         for iline=0, n_elements(linelist)-1 do begin
            obswave = linelist[iline].lambda * (1 + zguess)
            if (obswave LE minwave OR obswave GE maxwave) then begin
               thismask = objloglam GT alog10(obswave) - backwidth $
                AND objloglam LT alog10(obswave) + backwidth
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
         if (nflat GT 0) then begin
            if (nflat EQ 1) then allflat = flatback NE 0 $
             else allflat = total(flatback,2) NE 0
            background = background * (allflat EQ 0)
         endif

         ; Combine the simple background terms that overlap one another
         if (nflat GT 1) then begin
            ; First smush all these terms together
            allflat = long(total(flatback,2) NE 0)

            ; Now split them back up again, and append to the "background".
            ; Normalize the level to the median level in the object flux.
            allflat = [0, allflat, 0]
            dflat = allflat[1:npixobj+1] - allflat[0:npixobj]
            istart = where(dflat EQ 1, nflat)
            iend = where(dflat EQ -1, nend) - 1

            for iflat=0, nflat-1 do begin
               thismask = lonarr(npixobj)
               thismask[istart[iflat]:iend[iflat]] = $
                djs_median( objflux[istart[iflat]:iend[iflat],iobj] )
               background = [[background], [thismask]]
            endfor
         endif

         ; Call the line-fitting engine
         if (size(background,/n_dimen) EQ 1) then thismask = background NE 0 $
          else thismask = total(background,2) NE 0
         lfit1 = linebackfit(linelist.lambda, objloglam, $
          objflux[*,iobj], invvar=objivar[*,iobj] * thismask, $
          linename=linelist.name, $
          zindex=linelist.zindex, windex=linelist.windex, $
          findex=linelist.findex, fvalue=linelist.fvalue, zguess=zguess, $
          background=background, yfit=yfit1, bfit=bfit, bterms=bterms)

         ; Fill the two output structures
         lfitall[*,iobj] = lfit1

         splog, 'Object #', iobj+1, ' CPU time for line fitting = ', $
          systime(1)-t0

         if (keyword_set(doplot)) then begin
            igood = where(objivar[*,iobj] GT 0, ngood)
            if (ngood GT 1) then begin
               djs_plot, 10^objloglam[igood], objflux[igood,iobj], $
                xtitle='Wavelength [Ang]', ytitle='Flux'
               djs_oplot, 10^objloglam[igood], yfit1[igood], color='red'

               ; Wait for a keystroke...
               if (keyword_set(debug)) then begin
                  print, 'Press any key...'
                  cc = strupcase(get_kbrd(1))
               endif
            endif
         endif
      endif else begin
         splog, 'Skipping object #', iobj+1
      endelse
   endfor

   ;----------
   ; Concatenate the two output structures into a single structure

   lfitall = struct_addtags(res_prepend, lfitall)

   ;----------
   ; Write the output file

   sxaddpar, zhdr, 'VERSLINE', idlspec2d_version(), $
    'Version of idlspec2d for line fitting', after='VERS1D'

   cpbackup, zlinefile
   mwrfits, 0, zlinefile, zhdr, /create ; Retain original header in first HDU
   mwrfits, lfitall, zlinefile

   ;----------
   ; Close log file

   if (keyword_set(plotfile)) then dfpsclose

   splog, 'Total time for SPECLINEFIT = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPECLINEFIT at ', systime()
   splog, /close

   return
end
;------------------------------------------------------------------------------
