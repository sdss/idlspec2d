;+
; NAME:
;   make_regress1d
;
; PURPOSE:
;   Generate a regression table from spectro 1-D from Chicago-1D outputs files.
;
; CALLING SEQUENCE:
;   make_regress1d, [listfile, indir=, outfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   listfile   - List with plate number and MJD (one such pair per line)
;                for plates to use.  Default to the file
;                $IDLSPEC2D_DIR/etc/regress1d_all.plates, which lists plates
;                300 to 309.
;   indir      - Input directory for the Chicago-1D output files, of the
;                format 'spDiag1d-'+mjd+'-'+plate+'.par, i.e.
;                'spDiag1d-51666-0300.par'; default to '.'
;   outfile    - Name of output file; default to
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   An ASCII file with minimal object information and redshifts is tabulated
;   from the Chicago-1D outputs.
;
;   The Chicago outputs are changed in some special cases that are assumed
;   to be missing or wrong.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_readcol
;
; INTERNAL SUPPORT ROUTINES:
;   create_zregress()
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/regress1d_all.plates
;
; REVISION HISTORY:
;   27-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function create_zregress, num

   ftemp = create_struct( name='ZREGRESS', $
    'PLATE', 0L, $
    'MJD', 0L, $
    'FIBERID', 0L, $
    'Z', 0.D, $
    'CLASS', '', $
    'PRIMTARGET', 0L, $
    'SECTARGET', 0L, $
    'COMMENTS', '' )

   zregress = replicate(ftemp, num)

   return, zregress
end

;------------------------------------------------------------------------------
pro make_regress1d, listfile, indir=indir, outfile=outfile

   if (NOT keyword_set(listfile)) then $
    listfile = filepath('regress1d_all.plates', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(outfile)) then outfile = 'regress1d_all.dat'

   djs_readcol, listfile, platevec, mjdvec

   nplate = n_elements(platevec)

   zregress = create_zregress(nplate * 640)

   ii = 0
   for iplate=0, nplate-1 do begin
      ; Read the Yanny param file for Chicago-1D outputs for this plate
      mjdstr = string(mjdvec[iplate],format='(i5.5)')
      platestr = string(platevec[iplate],format='(i4.4)')
      filename = filepath('spDiag1d-' + mjdstr + '-' + platestr + '.par', $
       root_dir=indir)
      yanny_read, filename, tt, /quick
      cdiag = *tt[0]
      yanny_free, tt

      ; Change the CLASS for 'SPEC_HIZ_QSO' to 'SPEC_QSO'
      kk = where(cdiag.class EQ 'SPEC_HIZ_QSO')
      if (kk[0] NE -1) then cdiag[kk].class = 'SPEC_QSO'

      for ifiber=0, 639 do begin
         zregress[ii].plate = platevec[iplate]
         zregress[ii].mjd = mjdvec[iplate]
         zregress[ii].fiberid = ifiber + 1
         jj = (where(cdiag.fiberid EQ ifiber+1))[0]
         if (jj NE -1) then begin
            zregress[ii].primtarget = cdiag[jj].primtarget
            zregress[ii].sectarget = cdiag[jj].sectarget
            zregress[ii].z = cdiag[jj].zfinal
            zregress[ii].class = '  ' + strmid(cdiag[jj].class,5) + '          '
         endif else begin
            zregress[ii].class = '  MISSING          '
         endelse

         ; Find out if this is a sky fiber, and change CLASS to 'SKY'
         if (zregress[ii].sectarget AND 2^4) then begin
            zregress[ii].z = 0
            zregress[ii].class = '  SKY'
         endif

         ; Declare as 'UNKNOWN' the following cases:
         ; Any object with redshift less than -1000
         ; A QSO with a redshift less than 1000
         if (zregress[ii].z LT -1000./3e5 OR $
          (zregress[ii].z LT 1000./3d5 AND $
           strtrim(zregress[ii].class,2) EQ 'QSO')) then begin
            zregress[ii].z = 0
            zregress[ii].class = '  UNKNOWN'
         endif

         ii = ii + 1
      endfor
   endfor

   get_lun, olun
   openw, olun, outfile
   printf, olun, $
    '#PLATE MJD FIBER  Z        CLASS     PRIMTARGET   SECTARGET  COMMENTS'
   for jj=0, n_elements(zregress)-1 do begin
      printf, olun, zregress[jj], format='(i5,i6,i4,f10.5,a10,i12,i12,a)'
   endfor
   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------
