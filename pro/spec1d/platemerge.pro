;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with tsObj files.
;
; CALLING SEQUENCE:
;   platemerge, [zfile, outfile=, ascfile=, /qsurvey]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   zfile       - Redshift file(s) from spectro-1D; default to all files
;                 specified by the PLATELIST routine.
;   outfile     - Output FITS file name; default to '$SPECTRO_DATA/spAll.fits'
;   ascfile     - Output ASCII file name; default to '$SPECTR_DATA/spAll.dat'
;   qsurvey     - If set, then limit to plates in platelist with QSURVEY=1.
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
pro platemerge, zfile, outfile=outfile, ascfile=ascfile, qsurvey=qsurvey

   if (NOT keyword_set(outfile)) then $
    outfile = djs_filepath('spAll.fits', root_dir=getenv('SPECTRO_DATA'))
   if (NOT keyword_set(ascfile)) then $
    ascfile = djs_filepath('spAll.dat', root_dir=getenv('SPECTRO_DATA'))

   ;----------
   ; Find the list of spZ files.

   if (NOT keyword_set(zfile)) then begin
      platelist, plist=plist
      if (NOT keyword_set(plist)) then return
      if (keyword_set(qsurvey)) then begin
         indx = where(plist.qsurvey AND strtrim(plist.status1d,2) EQ 'Done')
      endif else begin
         indx = where(strtrim(plist.status1d,2) EQ 'Done')
      endelse
      if (indx[0] EQ -1) then return
      plist = plist[indx]
      nfile = n_elements(plist)
      fullzfile = strarr(nfile)
      fullzfile = 'spZbest-' + string(plist.plate, format='(i4.4)') $
       + '-' + string(plist.mjd, format='(i5.5)') + '.fits'
      zsubdir = string(plist.plate, format='(i4.4)')
      for i=0, nfile-1 do $
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
   ; Loop through each file

   for ifile=0, nfile-1 do begin
      print,'File ',ifile, ' of ', nfile,': '+fullzfile[ifile]

      hdr = headfits(fullzfile[ifile])
      plate = sxpar(hdr, 'PLATEID')
      zans = mrdfits(fullzfile[ifile], 1)
      tsobj = plug2tsobj(plate, zans.plug_ra, zans.plug_dec)
      if (NOT keyword_set(tsobj)) then $
       splog, 'WARNING: No tsObj file found for plate ', plate

      if (ifile EQ 0) then begin
         outdat = replicate( create_struct(zans[0], tsobj[0]), nout)
         struct_assign, {junk:0}, outdat ; Zero-out all elements
      endif

      indx = lindgen(640)+640L*ifile
      copy_struct_inx, zans, outdat, index_to=indx
      if (keyword_set(tsobj)) then $
       copy_struct_inx, tsobj, outdat, index_to=indx

      ; Over-write PRIMTARGET+SECTARGET with those values from spPlate file.
      plugmap = mrdfits(fullplatefile[ifile], 5)
      outdat[indx].primtarget = plugmap.primtarget
      outdat[indx].sectarget = plugmap.sectarget
   endfor

   ;----------
   ; Write the output file

   mwrfits, outdat, outfile, /create

   ;----------
   ; Create the structure for ASCII output

   adat = create_struct( $
    'plate'      ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'rchi2'      , 0.0, $
    'dof'        ,  0L, $
    'ra'         , 0.0d, $
    'dec'        , 0.0d, $
    'plate_sn2'  ,  0.0, $
    'counts_model', fltarr(5), $
    'objc_type'  ,  '', $
    'primtarget' ,  0L, $
    'sectarget'  ,  0L )
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

   struct_print, adat, filename=ascfile

   return
end
;------------------------------------------------------------------------------
