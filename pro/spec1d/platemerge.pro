;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with tsObj files.
;
; CALLING SEQUENCE:
;   platemerge, [zfile, outfile=, ascfile=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   zfile       - Redshift file(s) from spectro-1D; default to all files
;                   matching '*/spZbest*.fits'
;   outfile     - Output FITS file name; default to 'spAll.fits'
;   ascfile     - Output ASCII file name; default to 'spAll.dat'
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
;   struct_print
;   sxpar()
;
; REVISION HISTORY:
;   30-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platemerge, zfile, outfile=outfile, ascfile=ascfile

   if (NOT keyword_set(zfile)) then zfile = '*/spZbest*.fits'

   fullzfile = findfile(zfile, count=nfile)
   print, 'Found ', nfile, ' files'
   if (nfile EQ 0) then return
   fullzfile = fullzfile[ sort(fullzfile) ]

   if (NOT keyword_set(outfile)) then outfile = 'spAll.fits'
   if (NOT keyword_set(ascfile)) then ascfile = 'spAll.dat'

   nout = nfile * 640L
   print, 'Total number of objects = ', nout

   ;----------
   ; Loop through each file

   for ifile=0, nfile-1 do begin
      print,'File ',ifile, ' of ', nfile,': '+fullzfile[ifile]

      hdr = headfits(fullzfile[ifile])
      plate = sxpar(hdr, 'PLATEID')
      zans = mrdfits(fullzfile[ifile], 1)
      tsobj = plug2tsobj(plate, zans.plug_ra, zans.plug_dec)

      if (ifile EQ 0) then begin
         outdat = replicate( create_struct(zans[0], tsobj[0]), nout)
         struct_assign, {junk:0}, outdat ; Zero-out all elements
      endif

      copy_struct_inx, zans, outdat, index_to=lindgen(640)+640L*ifile
      if (keyword_set(tsobj)) then $
       copy_struct_inx, tsobj, outdat, index_to=lindgen(640)+640L*ifile

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
    'modelcounts', fltarr(5), $
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

   struct_print, adat, filename=ascfile

   return
end
;------------------------------------------------------------------------------
