;+
; NAME:
;   findspec
;
; PURPOSE:
;   Routine for finding SDSS spectra that match a given RA, DEC.
;
; CALLING SEQUENCE:
;   findspec, [ra, dec, infile=, outfile=, searchrad=, slist=, $
;    topdir=, run2d=, run1d=, /best, /print, /sdss ]
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
;   topdir     - Optional override value for the environment
;                variable $BOSS_SPECTRO_REDUX.
;   run2d      - Override environment variable RUN2D, searching only for
;                those RUN2D reductions
;   run1d      - Override environment variable RUN1D, searching only for
;                those RUN1D reductions
;   best       - If set, then return the best match for each location, where
;                best is defined to be the closest object on the plate with
;                the best S/N.
;                This also forces the return of one structure element in SLIST
;                per position, so that you get exactly a paired list between
;                RA,DEC and SLIST.
;   print      - If set, then print matches to the terminal.
;   sdss       - This is a shortcut option that is exactly equivalent to
;                topdir=GETENV('SPECTRO_REDUX'), run2d='26'
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
;   findspec.pro is currently setup to use BOSS data, but it still works
;   with SDSS-I,II data with the following procedure:
;
;   findspec, ra, dec, topdir=GETENV('SPECTRO_REDUX'), run2d='26'
;
;   or
;
;   findspec, ra, dec, /sdss
;
; EXAMPLES:
;   Make a file "file.in" with the following two lines:
;     218.7478    -0.3745007
;     217.7803    -0.8900855
;
;   Then run the command:
;     IDL> findspec,infile='file.in',/sdss,/print
;
;   This should print:
;
;PLATE   MJD FIBERID            RA            DEC    MATCHRAD
;----- ----- ------- ------------- -------------- -----------
;  306 51637     101      218.7478     -0.3745007      0.0000
;  306 51690     117      218.7478     -0.3745007      0.0000
;  306 51637     201      217.7803     -0.8900855  2.9996e-05
;  306 51690     217      217.7803     -0.8900855  2.9996e-05
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
 slist=slist, topdir=topdir1, run2d=run2d1, run1d=run1d1, $
 best=best, print=print1, sdss=sdss

   common com_findspec, plist, nlist, lasttopdir

   if (keyword_set(sdss)) then begin
      topdir1 = getenv('SPECTRO_REDUX')
      run2d1 = '26'
   endif
   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d
   if (keyword_set(run1d1)) then run1d = strtrim(run1d1,2) $
    else run1d = getenv('RUN1D')
   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')

   if (keyword_set(lasttopdir) EQ 0) then lasttopdir = ''
   if (keyword_set(plist) EQ 0 OR topdir NE lasttopdir) then begin
      lasttopdir = topdir
      ; platelist, plist=plist, topdir=topdir
      platelist_files = file_search(topdir+'/platelist.fits')
      plates_files = file_search(topdir+'/plates-*.fits')
      if (strlen(platelist_files[0]) GT 0) then $
         plist = mrdfits(platelist_files[0],1,/silent)
      if (strlen(plates_files[0]) GT 0) then $
         plist = mrdfits(plates_files[0],1,/silent)
      if (NOT keyword_set(plist)) then $
       message, 'Plate list (platelist.fits or plates-*.fits) not found in ' + topdir
      nlist = n_elements(plist)
   endif

   qdone = strmatch(plist.status1d,'Done*')
   if (tag_exist(plist, 'RUN2D')) then begin
      qdone2d = bytarr(nlist)
      for i=0, n_elements(run2d)-1 do $
       qdone2d = qdone2d OR strmatch(strtrim(plist.run2d,2), run2d[i])
      qdone = qdone AND qdone2d
   endif
   if (tag_exist(plist, 'RUN1D')) then begin
      qdone1d = bytarr(nlist)
      for i=0, n_elements(run1d)-1 do $
       qdone1d = qdone1d OR strmatch(strtrim(plist.run1d,2), run1d[i])
      qdone = qdone AND qdone1d
   endif
   idone = where(qdone, ndone)
   if (ndone EQ 0) then begin
      splog, 'No reduced plates!'
      slist = 0
      return
   endif

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
   spherematch, ra, dec, plist[idone].racen, plist[idone].deccen, $
    searchrad+1.55, imatch1, itmp, dist12, maxmatch=0
   if (imatch1[0] EQ -1) then return
   imatch2 = idone[itmp]

   ;----------
   ; Read all relevant plates

   ; Read the number of fibers per plate, or assume 640 if old SDSS-I platelist
   if (tag_exist(plist,'N_TOTAL')) then n_total = plist.n_total $
    else n_total = replicate(640L, n_elements(plist))

   iplate = imatch2[uniq(imatch2,sort(imatch2))]
   i0 = 0L
   for i=0L, n_elements(iplate)-1L do begin
      readspec, plist[iplate[i]].plate, mjd=plist[iplate[i]].mjd, $
       plugmap=plugmap1, topdir=topdir, run2d=run2d, /silent
      if (n_elements(plugmap1) EQ 1) then $
       message, 'Failed reading plate ='+strtrim(plist[iplate[i]].plate,2) $
        + ' MJD='+strtrim(plist[iplate[i]].mjd,2)
      if (i EQ 0) then plugmap = replicate(create_struct('PLATE',0L,'MJD',0L, $
       'FIBERID',0L,'RA',0d,'DEC',0d), total(n_total[iplate]))
      index_to=i0+lindgen(n_total[iplate[i]])
      plugmap[index_to].plate = plist[iplate[i]].plate
      plugmap[index_to].mjd = plist[iplate[i]].mjd
      plugmap[index_to].fiberid = plugmap1.fiberid
      plugmap[index_to].ra = plugmap1.ra
      plugmap[index_to].dec = plugmap1.dec
;      copy_struct_inx, plugmap1, plugmap, index_to=index_to

      i0 += n_total[iplate[i]]
   endfor
   spherematch, ra, dec, plugmap.ra, plugmap.dec, searchrad, $
    i1, i2, d12, maxmatch=0
   if (i1[0] EQ -1) then return

   if (NOT keyword_set(best)) then begin

      ;------------------------------------------------------------------------
      ; RETURN ALL MATCHES
      ;------------------------------------------------------------------------

      slist = replicate(blanklist, n_elements(i1))
      slist.plate = plugmap[i2].plate
      slist.mjd = plugmap[i2].mjd
      slist.fiberid = plugmap[i2].fiberid
      slist.ra = plugmap[i2].ra
      slist.dec = plugmap[i2].dec
      slist.matchrad = d12

   endif else begin
      ;------------------------------------------------------------------------
      ; RETURN ONLY BEST MATCH PER OBJECT
      ;------------------------------------------------------------------------

      ; Read the median S/N for each spectrum
      readspec, plugmap[i2].plate, mjd=plugmap[i2].mjd, $
       plugmap[i2].fiberid, zans=zans, topdir=topdir, $
       run2d=run2d, run1d=run1d, /silent

      ; We have all the possible matches.  Now sort those in the order
      ; of each input object, but with the last entry being the best
      ; according to the SN_MEDIAN value.
      isort = sort(i1 + (zans.sn_median>0)/max(zans.sn_median+1.))
      i1 = i1[isort]
      i2 = i2[isort]
      d12 = d12[isort]

      iuniq = uniq(i1)
      slist[i1[iuniq]].plate = plugmap[i2[iuniq]].plate
      slist[i1[iuniq]].mjd = plugmap[i2[iuniq]].mjd
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
