;+
; NAME:
;   findspec
;
; PURPOSE:
;   Routine for finding SDSS spectra that match a given RA, DEC.
;
; CALLING SEQUENCE:
;   findspec, [ra, dec, infile=, outfile=, slist=, /silent ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   ra         - Right ascension; scalar or array in degrees.
;   dec        - Declination; scalar or array in degrees.
;   infile     - Input file with RA, DEC positions, one per line.
;                If set, then this over-rides values passed in RA,DEC.
;   outfile    - If set, then print matches to this file.
;   silent     - If set, then suppress printing outputs to the terminal.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   slist      - Structure containing information for each match.
;
; COMMENTS:
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
pro findspec, ra, dec, infile=infile, outfile=outfile, slist=slist,  $
 silent=silent

   common com_findspec, plist, nlist, minsn

   if (NOT keyword_set(plist)) then begin
      platelist, plist=plist
      nlist = n_elements(plist)
      minsn = fltarr(nlist)
      for i=0, nlist-1 do $
      minsn[i] = min([ plist[i].sn2_g1, plist[i].sn2_g2, $
       plist[i].sn2_i1, plist[i].sn2_i2 ])
   endif

   ;----------
   ; Read an input file if specified

   if (keyword_set(infile)) then begin
      djs_readcol, infile, ra, dec, format='(D,D)'
   endif

   ;----------
   ; Call this routine recursively if RA,DEC are arrays

   nvec = n_elements(ra)
   if (nvec GT 1) then begin
      for i=0, nvec-1 do begin
         findspec, ra[i], dec[i], slist=slist1, /silent
         if (i EQ 0) then slist = slist1 $
          else slist = [slist, slist1]
      endfor
      if (NOT keyword_set(silent)) then struct_print, slist
      if (keyword_set(outfile)) then struct_print, slist, filename=outfile
      return
   endif

   ;----------
   ; Create output structure

   slist = create_struct(name='slist', $
    'plate'   , 0L, $
    'mjd'     , 0L, $
    'fiberid' , 0L, $
    'ra'      , 0.d, $
    'dec'     , 0.d )
   slist = replicate(slist, nlist)

   ;----------
   ; Loop through each possible plate, looking for nearest object

   adist = djs_diff_angle(ra, dec, plist.ra, plist.dec)

   for iplate=0, nlist-1 do begin
      if (adist[iplate] LT 1.55) then begin
         readspec, plist[iplate].plate, mjd=plist[iplate].mjd, plugmap=plugmap
         objdist = djs_diff_angle(ra, dec, plugmap.ra, plugmap.dec)
         mindist = min(objdist, imin)
         if (mindist LT 3.0/3600.) then begin
            slist[iplate].plate = plist[iplate].plate
            slist[iplate].mjd = plist[iplate].mjd
            slist[iplate].fiberid = plugmap[imin].fiberid
            slist[iplate].ra = plugmap[imin].ra
            slist[iplate].dec = plugmap[imin].dec
         endif
      endif
   endfor

   ;----------
   ; Select the exposure on the plate with the best S/N

   indx = where(slist.plate NE 0)
   if (indx[0] EQ -1) then begin
      slist = slist[0]
      slist.ra = ra
      slist.dec = dec
      indx = 0
   endif

   junk = max(minsn[indx], ibest)
   slist = slist[indx[ibest]]

   if (NOT keyword_set(silent)) then struct_print, slist
   if (keyword_set(outfile)) then struct_print, slist, filename=outfile

   return
end
;------------------------------------------------------------------------------
