;+
; NAME:
;   spmanual
;
; PURPOSE:
;   Read manual inspection files (spInspect) for redshifts
;
; CALLING SEQUENCE:
;   zmanual = spmanual(plateid, mjd= )
;
; INPUTS:
;   plateid    - Plate number (scalar)
;   mjd        - MJD of observation
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   zmanual    - Structure with manual inspection data, with one entry
;                per fiber even if the manual file is incomplete
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Read manual inspection files in:
;    $SPINSPECT_DIR/data/*/spInspect-$PLATE-$MJD.par
;
; EXAMPLES:
;
; BUGS:
;   Currently only reads the French Participation Group (fpg) files.
;   Output structure will contain whatever the input files contain.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   09-Feb-2010  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
function spmanual, plateid, mjd=mjd

   common com_spmanual, plates

   if (n_elements(plateid) NE 1 OR n_elements(mjd) NE 1) then $
    message, 'PLATEID and MJD must be scalars'

   root_dir = getenv('SPINSPECT_DIR')

   platestr = string(plateid, format='(i4.4)')
   mjdstr = string(mjd, format='(i5.5)')
   files = findfile(djs_filepath('spInspect-'+platestr+'-'+mjdstr+'.par', $
    root_dir=root_dir, subdir=['data','fpg']), count=nfile)
   if (nfile EQ 0) then return, 0

   ;----------
   ; Read the manual inspection file

   bossobject = yanny_readone(files, /anon)
   nthis = n_elements(bossobject)
   if (size(bossobject,/tname) NE 'STRUCT') then return, 0

   ;----------
   ; If exactly 1000 entries, assume 1000 fibers; otherwise read the platePlans
   ; file to determine the number of fibers

   if (nthis EQ 1000) then begin
      nfiber = 1000
   endif else begin
      if (NOT keyword_set(plates)) then begin
         plates = yanny_readone(djs_filepath('platePlans.par', $
          root_dir=getenv('PLATELIST_DIR')))
      endif
      if (NOT keyword_set(plates)) then $
       message, 'Unable to read platePlans.par file'
      i = where(plates.plateid EQ plateid, ct)
      if (ct EQ 0) then begin
         nfiber = 1000
      endif else begin
         case plates[i[0]].survey of
         'boss': nfiber = 1000
         'marvels_pre': nfiber = 1000
         'sdss': nfiber = 640
         'segue1': nfiber = 640
         'segue2': nfiber = 640
         else: message, 'Unsupported survey type '+plates[i[0]].survey
         endcase
      endelse
   endelse

   ;----------
   ; Copy into the output structure

   blankdat = bossobject[0]
   struct_assign, {junk:0}, blankdat
   if (tag_exist(blankdat,'PLATE') EQ 0) then $
    blankdat = create_struct(blankdat, {plate:0L})
   if (tag_exist(blankdat,'MJD') EQ 0) then $
    blankdat = create_struct(blankdat, {mjd:0L})
   if (tag_exist(blankdat,'FIBERID') EQ 0) then $
    blankdat = create_struct(blankdat, {fiberid:0L})
   zmanual = replicate(blankdat, nfiber)
   zmanual[bossobject.fiberid-1] = bossobject
   zmanual.plate = plateid
   zmanual.mjd = mjd
   zmanual.fiberid = lindgen(nfiber) + 1

   return, zmanual
end
;------------------------------------------------------------------------------
