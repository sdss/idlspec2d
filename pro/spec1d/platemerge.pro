;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with tsObj files.
;
; CALLING SEQUENCE:
;   platemerge, [zfile, outroot=, /public]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   zfile       - Redshift file(s) from spectro-1D; default to all files
;                 specified by the PLATELIST routine.
;   outroot     - Root name for output files; default to '$SPECTRO_DATA/spAll';
;                 the files are then 'spAll.fits' and 'spAll.dat'.
;                 If /PUBLIC is set, then add '-public' to the root name.
;   public      - If set, then limit to plates in platelist with PUBLIC != ''
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
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   headfits()
;   mrdfits()
;   mwrfits
;   plug2tsobj()
;   platelist
;   struct_print
;   sxpar()
;
; REVISION HISTORY:
;   30-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platemerge, zfile, outroot=outroot, public=public

   if (NOT keyword_set(outroot)) then begin
      outroot = 'spAll'
      if (keyword_set(public)) then outroot = outroot + '-public'
      outroot = djs_filepath(outroot, root_dir=getenv('SPECTRO_DATA'))
   endif

   ;----------
   ; Find the list of spZ files.

   if (NOT keyword_set(zfile)) then begin
      platelist, plist=plist
      if (NOT keyword_set(plist)) then return

      indx = where(strtrim(plist.status1d,2) EQ 'Done' AND $
       (strtrim(plist.platequality,2) EQ 'good' $
       OR strtrim(plist.platequality,2) EQ 'marginal'))
      if (indx[0] EQ -1) then return
      if (keyword_set(public)) then $
       indx = indx[ where(strtrim(plist[indx].public)) ]
      if (indx[0] EQ -1) then return
      plist = plist[indx]

      nfile = n_elements(plist)
      fullzfile = strarr(nfile)
      fullzfile = 'spZbest-' + string(plist.plate, format='(i4.4)') $
       + '-' + string(plist.mjd, format='(i5.5)') + '.fits'
      zsubdir = string(plist.plate, format='(i4.4)')
      for i=0L, nfile-1 do $
       fullzfile[i] = djs_filepath(fullzfile[i], $
        root_dir=getenv('SPECTRO_DATA'), subdirectory=zsubdir[i])
   endif else begin
      fullzfile = findfile(zfile, count=nfile)
   endelse

   print, 'Found ', nfile, ' files'
   if (nfile EQ 0) then return
   fullzfile = fullzfile[ sort(fullzfile) ]

   nout = nfile * 640L
   print, 'Total number of objects = ', nout

   ;----------
   ; Find the corresponding spPlate files (needed only for PRIMTARGET+SECTARGET
   ; flags, which are incorrect in the tsObj files.

   fullplatefile = repstr(fullzfile, 'spZbest', 'spPlate')

   ;----------
   ; Find the first tsObj file that exists for use in constructing the
   ; output structure.

   ifile = 0
   while (NOT keyword_set(tsobj0)) do begin
      tsobj0 = plug2tsobj(plist[ifile].plate, 0, 0)
      ifile = ifile + 1
      if (ifile EQ nfile) then $
       message, 'No tsObj files found!'
   endwhile

   ;----------
   ; Loop through each file

   for ifile=0L, nfile-1 do begin
      print,'File ',ifile+1, ' of ', nfile,': '+fullzfile[ifile]

      hdr = headfits(fullzfile[ifile])
      plate = sxpar(hdr, 'PLATEID')
      zans = mrdfits(fullzfile[ifile], 1)
      tsobj = plug2tsobj(plate, zans.plug_ra, zans.plug_dec)
      if (NOT keyword_set(tsobj)) then $
       splog, 'WARNING: No tsObj file found for plate ', plate

      if (NOT keyword_set(outdat)) then begin
         pstuff = create_struct( $
          'progname'    , ' ', $
          'chunkname'   , ' ', $
          'platequality', ' ', $
          'platesn2'    , 0.0, $
          'specprimary' ,  0L )
         outdat = replicate(create_struct(pstuff, zans[0], tsobj0), nout)
         struct_assign, {junk:0}, outdat ; Zero-out all elements
      endif

      indx = lindgen(640)+640L*ifile
      copy_struct_inx, zans, outdat, index_to=indx
      if (keyword_set(tsobj)) then $
       copy_struct_inx, tsobj, outdat, index_to=indx

      ; Fill in the first columns of this output structure
      outdat[indx].progname = plist[ifile].progname
      outdat[indx].chunkname = plist[ifile].chunkname
      outdat[indx].platequality = plist[ifile].platequality
      outdat[indx].platesn2 = plist[ifile].platesn2

      ; Over-write PRIMTARGET+SECTARGET with those values from spPlate file.
      plugmap = mrdfits(fullplatefile[ifile], 5)
      outdat[indx].primtarget = plugmap.primtarget
      outdat[indx].sectarget = plugmap.sectarget

      ; Over-write the MJD with that from the plate file name ???
      ; Early versions of 2D (such as v4_3_1) could have an inconsistent value.
      thismjd = long( strmid(fileandpath(fullplatefile[ifile]), 13, 5) )
      outdat[indx].mjd = thismjd
   endfor

   ;----------
   ; Write the output file

   mwrfits, outdat, outroot+'.fits', /create

   ;----------
   ; Create the structure for ASCII output

   adat = create_struct( $
    'plate'      ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'zwarning'   ,  0L, $
    'rchi2'      , 0.0, $
    'ra'         , 0.0d, $
    'dec'        , 0.0d, $
    'plate_sn2'  ,  0.0, $
    'counts_model', fltarr(5), $
    'objc_type'  ,  '', $
    'primtarget' ,  0L, $
    'sectarget'  ,  0L, $
    'objtype'    ,  '' )
   adat = replicate(adat, nout)
   struct_assign, outdat, adat

   adat.plate_sn2 = (outdat.spec1_g + outdat.spec1_i + outdat.spec2_g $
    + outdat.spec2_i) / 4.0

   ii = where(strtrim(adat.class,2) EQ '')
   if (ii[0] NE -1) then adat[ii].class = '""'

   ii = where(strtrim(adat.subclass,2) EQ '')
   if (ii[0] NE -1) then adat[ii].subclass = '""'

   objtypes = ['UNKNOWN', 'CR', 'DEFECT', 'GALAXY', 'GHOST', 'KNOWNOBJ', $
    'STAR', 'TRAIL', 'SKY']
   adat.objc_type = objtypes[outdat.objc_type]
outdat = 0 ; Free memory

   struct_print, adat, filename=outroot+'.dat'

   return
end
;------------------------------------------------------------------------------
