;+
; NAME:
;   findspec
;
; PURPOSE:
;   Routine for finding SDSS spectra that match a given RA, DEC.
;
; CALLING SEQUENCE:
;   findspec, [ra, dec, infile=, outfile=, searchrad=, slist=, $
;    /best, /print ]
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
;   print      - If set, then print matches to the terminal.
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
 slist=slist, best=best, print=print1

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
   idone = where(strmatch(plist.status1d,'Done*'), ndone)
   if (ndone EQ 0) then $
    message, 'No reduced plates!'

   ;----------
   ; Read an input file if specified

   if (keyword_set(infile)) then begin
      djs_readcol, infile, ra, dec, format='(D,D)'
   endif
   nobj = n_elements(ra)

   if (NOT keyword_set(searchrad)) then searchrad = 3./3600.

   ;----------
   ; Create output structure

   blanklist = create_struct(name='slist', $
    'plate'   , 0L, $
    'mjd'     , 0L, $
    'fiberid' , 0L, $
    'ra'      , 0.d, $
    'dec'     , 0.d, $
    'matchrad', 0.0 )

   ; Set default return values
   if (keyword_set(best)) then slist = replicate(blanklist, nobj) $
    else slist = 0

   ;----------
   ; Match all plates with objects

   nplate = n_elements(plist)
   spherematch, ra, dec, plist[idone].ra, plist[idone].dec, searchrad+1.55, $
    imatch1, itmp, dist12, maxmatch=0
   if (imatch1[0] EQ -1) then return
   imatch2 = idone[itmp]

   ;----------
   ; Read all relevant plates

   iplate = imatch2[uniq(imatch2,sort(imatch2))]
   for i=0L, n_elements(iplate)-1L do begin
      readspec, plist[iplate[i]].plate, mjd=plist[iplate[i]].mjd, $
       plugmap=plugmap1, /silent
      if (i EQ 0) then plugmap = replicate(plugmap1[0], 640L*n_elements(iplate))
      copy_struct_inx, plugmap1, plugmap, index_to=i*640L+lindgen(640)
   endfor
   spherematch, ra, dec, plugmap.ra, plugmap.dec, searchrad, $
    i1, i2, d12, maxmatch=0
   if (i1[0] EQ -1) then return
   i2plate = iplate[long(i2/640)] ; index within the PLIST structure

   if (NOT keyword_set(best)) then begin

      ;------------------------------------------------------------------------
      ; RETURN ALL MATCHES
      ;------------------------------------------------------------------------

      slist = replicate(blanklist, n_elements(i1))
      slist.plate = plist[i2plate].plate
      slist.mjd = plist[i2plate].mjd
      slist.fiberid = plugmap[i2].fiberid
      slist.ra = plugmap[i2].ra
      slist.dec = plugmap[i2].dec
      slist.matchrad = d12

   endif else begin
      ;------------------------------------------------------------------------
      ; RETURN ONLY BEST MATCH PER OBJECT
      ;------------------------------------------------------------------------

      ; Read the median S/N for each spectrum
      readspec, plist[i2plate].plate, mjd=plist[i2plate].mjd, $
       plugmap[i2].fiberid, zans=zans, /silent

      ; We have all the possible matches.  Now sort those in the order
      ; of each input object, but with the last entry being the best
      ; according to the SN_MEDIAN value.
      isort = sort(i1 + (zans.sn_median>0)/max(zans.sn_median+1.))
      i1 = i1[isort]
      i2 = i2[isort]
      i2plate = i2plate[isort]
      d12 = d12[isort]

      iuniq = uniq(i1)
      slist[i1[iuniq]].plate = plist[i2plate[iuniq]].plate
      slist[i1[iuniq]].mjd = plist[i2plate[iuniq]].mjd
      slist[i1[iuniq]].fiberid = plugmap[i2[iuniq]].fiberid
      slist[i1[iuniq]].ra = plugmap[i2[iuniq]].ra
      slist[i1[iuniq]].dec = plugmap[i2[iuniq]].dec
      slist[i1[iuniq]].matchrad = d12[iuniq]
   endelse

   ;----------
   ; Print to terminal and/or output file

   if (keyword_set(print1)) then struct_print, slist
   if (keyword_set(outfile)) then struct_print, slist, filename=outfile

   return
end
;------------------------------------------------------------------------------
