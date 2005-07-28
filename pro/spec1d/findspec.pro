;+
; NAME:
;   findspec
;
; PURPOSE:
;   Routine for finding SDSS spectra that match a given RA, DEC.
;
; CALLING SEQUENCE:
;   findspec, [ra, dec, infile=, outfile=, searchrad=, slist=, $
;    /best, /silent ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   ra         - Right ascension; scalar or array in degrees.
;   dec        - Declination; scalar or array in degrees.
;   infile     - Input file with RA, DEC positions, one per line.
;                If set, then this over-rides values passed in RA,DEC.
;   outfile    - If set, then print matches to this file.
;   searchrad  - Search radius in degrees; default to 3./3600 (3 arcsec).
;   best       - If set, then return the best match for each location, where
;                best is defined to be the closest object on the plate with
;                the best S/N.
;                This also forces the return of one structure element in SLIST
;                per position, so that you get exactly a paired list between
;                RA,DEC and SLIST.
;   silent     - If set, then suppress printing outputs to the terminal.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   slist      - Structure containing information for each match.
;
; COMMENTS:
;   The search radius is set to within 1.55 degress of a plate center,
;   then within 3 arcsec of an object.
;
; EXAMPLES:
;   Make a file "file.in" with the following two lines:
;     218.7478    -0.3745007
;     217.7803    -0.8900855
;
;   Then run the command:
;     IDL> findspec,infile='file.in'
;
;   This should print:
;     PLATE   MJD FIBERID            RA            DEC
;     ----- ----- ------- ------------- --------------
;       306 51637     101      218.7478     -0.3745007
;       306 51637     201      217.7803     -0.8900855
;
; BUGS:
;
; PROCEDURES CALLED:
;  djs_readcol
;  djs_diff_angle()
;  platelist
;  readspec
;  struct_print
;
; REVISION HISTORY:
;   15-Feb-2001  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro findspec, ra, dec, infile=infile, outfile=outfile, searchrad=searchrad, $
 slist=slist, best=best, silent=silent

   common com_findspec, plist, nlist, platesn

   if (NOT keyword_set(plist)) then begin
      platelist, plist=plist
      if (NOT keyword_set(plist)) then $
       message, 'Plate list (platelist.fits) not found in $SPECTRO_DATA'
      nlist = n_elements(plist)
      platesn = fltarr(nlist)
      for i=0L, nlist-1L do $
       platesn[i] = min([ plist[i].sn2_g1, plist[i].sn2_g2, $
        plist[i].sn2_i1, plist[i].sn2_i2 ])
   endif

   ;----------
   ; Read an input file if specified

   if (keyword_set(infile)) then begin
      djs_readcol, infile, ra, dec, format='(D,D)'
   endif

   if (NOT keyword_set(searchrad)) then searchrad = 3./3600.

   ;----------
   ; Call this routine recursively if RA,DEC are arrays

   nvec = n_elements(ra)
   if (nvec GT 1) then begin
      slist = 0
      for i=0L, nvec-1L do begin
         findspec, ra[i], dec[i], searchrad=searchrad, $
          slist=slist1, best=best, /silent
         if (keyword_set(slist1)) then $
          slist = keyword_set(slist) ? [slist, slist1] : slist1
      endfor
      if (NOT keyword_set(silent)) then struct_print, slist
      if (keyword_set(outfile)) then struct_print, slist, filename=outfile
      return
   endif

   ;----------
   ; Create output structure

   blanklist = create_struct(name='slist', $
    'plate'   , 0L, $
    'mjd'     , 0L, $
    'fiberid' , 0L, $
    'ra'      , 0.d, $
    'dec'     , 0.d, $
    'matchrad', 0.0 )

   ;----------
   ; Loop through each possible plate, looking for any matches.

   adist = djs_diff_angle(ra, dec, plist.ra, plist.dec)

   slist = 0
   for iplate=0L, nlist-1L do begin
      if (adist[iplate] LT 1.55 + searchrad) then begin
         readspec, plist[iplate].plate, mjd=plist[iplate].mjd, plugmap=plugmap
         if (keyword_set(plugmap)) then begin
            objdist = djs_diff_angle(ra, dec, plugmap.ra, plugmap.dec)
            ikeep = where(objdist LE searchrad, nkeep)
            if (nkeep GT 0) then begin
               slist1 = replicate(blanklist, nkeep)
               slist1.plate = plist[iplate].plate
               slist1.mjd = plist[iplate].mjd
               slist1.fiberid = plugmap[ikeep].fiberid
               slist1.ra = plugmap[ikeep].ra
               slist1.dec = plugmap[ikeep].dec
               slist1.matchrad = objdist[ikeep]
               if (keyword_set(slist)) then begin
                  slist = [slist, slist1]
                  snval = [snval, platesn[iplate]]
               endif else begin
                  slist = slist1
                  snval = platesn[iplate]
               endelse
            endif
         endif
      endif
   endfor

   ;----------
   ; Select the exposure on the plate with the best S/N

   if (keyword_set(best) AND n_elements(slist) GT 1) then begin
      ; First trim to the plate with the best S/N
      junk = max(snval, ibest)
      indx = where(slist.plate EQ slist[ibest].plate AND slist.mjd EQ slist[ibest].mjd)
      slist = slist[indx]

      ; Then trim to the closest object on that plate
      if (n_elements(slist) GT 1) then begin
         junk = min(slist.matchrad, ibest)
         slist = slist[ibest]
      endif
   endif
   if (keyword_set(best) AND NOT keyword_set(slist)) then slist = blanklist

   if (NOT keyword_set(silent)) then struct_print, slist
   if (keyword_set(outfile)) then struct_print, slist, filename=outfile

   return
end
;------------------------------------------------------------------------------
