;+
; NAME:
;   platelist
;
; PURPOSE:
;   Make list of reduced plates
;
; CALLING SEQUENCE:
;   platelist, [fullplatefile, outfile=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   infile      - Either a list of spPlancomb-*.par files, or a list of
;                 spPlate*.fits files; default to all files matching
;                 '$SPECTRO_DATA/*/spPlancomb-*.par'.
;   outfile     - If set, then write output to this file.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   If INFILE is a list of plan files, i.e.
;     spPlancomb-0306-51690.par
;   then look for the following files:
;     spPlate-0306-51690.fits (as specified by 'combinefile' in the plan file)
;     spZbest-0306-51690.fits
;
;   Otherwise, assume that INFILE is a list of FITS files (spPlate-*.fits),
;   and look for the redshift files (spZbest-*.fits).
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   fileandpath()
;   headfits()
;   mrdfits()
;   sxpar()
;   yanny_par()
;   yanny_read
;
; REVISION HISTORY:
;   29-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platelist, infile, outfile=outfile, plist=plist

   if (NOT keyword_set(infile)) then $
    infile = djs_filepath('spPlancomb-*.par', $
     root_dir=getenv('SPECTRO_DATA'), subdirectory='*')
   fullfile = findfile(infile, count=nfile)
   if (nfile EQ 0) then return

   ;----------
   ; Sort these files

   fullfile = fullfile[sort(fullfile)]

   ;----------
   ; Create output structure

   plist = create_struct( $
    'plateid' , 0L, $
    'mjd'     , 0L, $
    'ra'      , 0.0, $
    'dec'     , 0.0, $
    'snvec'   , fltarr(4), $
    'nums'    , lonarr(5) )
   plist = replicate(plist, nfile)

   ;----------
   ; Open output file

   if (keyword_set(outfile)) then $
    openw, olun, outfile, /get_lun $
   else $
    olun = -1L

   ;----------
   ; Loop through all files

   printf, olun, 'PLATE  MJD   RA    DEC   SN2_G1 SN2_I1 SN2_G2 SN2_I2 ' $
    + 'Ngal Nqso Nsta Nunk Nsky'
   printf, olun, '-----  ----- ----- ----- ------ ------ ------ ------ ' $
    + '---- ---- ---- ---- ----'

   for ifile=0, nfile-1 do begin

      junk = fileandpath(fullfile[ifile], path=path)  ; Determine PATH

      ;----------
      ; Test if INFILE specifies Yanny param files for spPlancomb.

      if (strmid(fullfile[ifile],strlen(fullfile[ifile])-4) EQ '.par') $
       then begin
         yanny_read, fullfile[ifile], hdr=hdrp
         platefile = $
          djs_filepath(yanny_par(hdrp, 'combinefile'), root_dir=path)
      endif else begin
         platefile = fullfile[ifile]
      endelse

      ;----------
      ; Determine names of associated files

      platemjd = strmid(fileandpath(platefile), 8, 10)
      zbestfile = 'spZbest-' + platemjd + '.fits'

      fullzfile = djs_filepath(zbestfile, root_dir=path)

      hdr1 = headfits(platefile)
      hdr2 = headfits(fullzfile)

      if (size(hdr1, /tname) EQ 'STRING') then begin
;         plateid = sxpar(hdr1, 'PLATEID')
;         mjd = sxpar(hdr1, 'MJD')
         ra = sxpar(hdr1, 'RA')
         dec = sxpar(hdr1, 'DEC')
         snvec = [ sxpar(hdr1, 'SPEC1_G'), $
                   sxpar(hdr1, 'SPEC1_I'), $
                   sxpar(hdr1, 'SPEC2_G'), $
                   sxpar(hdr1, 'SPEC2_I') ]
      endif else begin
         ra = 0
         dec = 0
         snvec = [0,0,0,0]
      endelse

      ; Get the following from the file names, since sometimes they
      ; are wrong in the file headers!!
      plateid = long( strmid(fileandpath(platefile), 8, 4) )
      mjd = long( strmid(fileandpath(platefile), 13, 5) )

      ;----------
      ; The RA,DEC in the header is sometimes wrong, so try to derive
      ; the field center from the plug-map information.  Choose the
      ; coordinates of the object closest to the center of the plate
      ; as defined by XFOCAL=YFOCAL=0.

      plug = mrdfits(platefile, 5)
      if (keyword_set(plug)) then begin
         iobj = where(strtrim(plug.holetype,2) EQ 'OBJECT')
         junk = min( plug[iobj].xfocal^2 + plug[iobj].yfocal^2, imin)
         ra = plug[iobj[imin]].ra
         dec = plug[iobj[imin]].dec
      endif

      if (size(hdr2, /tname) EQ 'STRING') then begin
         zans = mrdfits(fullzfile, 1, /silent)
         class = strtrim(zans.class,2)
         nums = [ total(class EQ 'GALAXY'), $
                  total(class EQ 'QSO'), $
                  total(class EQ 'STAR'), $
                  total(class EQ 'UNKNOWN'), $
                  total(class EQ 'SKY') ]
      endif else begin
         nums = lonarr(5)
      endelse

      printf, olun, plateid, mjd, ra, dec, snvec, nums, $
       format='(i5,i7,f6.1,f6.1,4f7.1,5i5)'
      flush, olun

      plist[ifile].plateid = plateid
      plist[ifile].mjd = mjd
      plist[ifile].ra = ra
      plist[ifile].dec = dec
      plist[ifile].snvec = snvec
      plist[ifile].nums = nums
   endfor

   ;----------
   ; Close output file

   if (keyword_set(outfile)) then begin
      close, olun
      free_lun, olun
   endif

end
;------------------------------------------------------------------------------
