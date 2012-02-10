;+
; NAME:
;   inspectgen
;
; PURPOSE:
;   Create empty spInspect file(s) for manual inspection of classifications
;
; CALLING SEQUENCE:
;   inspectgen, plate, mjd=, [ specinspect=, $
;    spectext=, specblend=, hdr=, structs= ]
;
; INPUTS:
;   plate      - Plate number
;   MJD        - MJD number
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   specinspect- SPECINSPECT structure (640 elements)
;   spectext   - SPECTEXT structure, with one empty element with fiberid=0
;   specblend  - SPECBLEND structure, with one empty element with fiberid=0
;   hdr        - String array with Yanny header
;   structs    - String array with Yanny structure definitions
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   readspec
;   sdss_flagval()
;
; REVISION HISTORY:
;   20-May-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro inspectgen, plate, mjd=mjd, specinspect=specinspect, $
 spectext=spectext, specblend=specblend, hdr=hdr, structs=structs

   ;----------
   ; Read the redshifts

   readspec, plate, mjd=mjd, zans=zans, zhdr=zhdr
   if (NOT keyword_set(zans)) then return
   nfiber = n_elements(zans)

   ;----------
   ; Create the output structures

   specinspect = create_struct( $
    name = 'SPECINSPECT', $
    'plate'          ,  long(plate), $
    'mjd'            ,  long(mjd), $
    'fiberid'        ,  0L, $
    'z_pipeline'     , 0.0, $
    'z_person'       , 0.0, $
    'z_conf_person'  ,  string(thisname), $
    'comments'       ,  '' )

   spectext = create_struct( $
    name = 'SPECTEXT', $
    'plate'          ,  long(plate), $
    'mjd'            ,  long(mjd), $
    'fiberid'        ,  0L, $
    'inspector'      ,  string(thisname), $
    'manual_text'    ,  '' )

   specblend = create_struct( $
    name = 'SPECBLEND', $
    'plate'          ,  long(plate), $
    'mjd'            ,  long(mjd), $
    'fiberid'        ,  0L, $
    'inspector'      ,  string(thisname), $
    'manual_class'   ,  '', $
    'manual_subclass',  '', $
    'manual_z',        0.0, $
    'manual_comments',  '' )

   ;----------
   ; Fill in the SPECINSPECT structure

   specinspect = replicate(specinspect, nfiber)
   specinspect.fiberid = zans.fiberid
   specinspect.class = strtrim(zans.class,2)
   specinspect.subclass = strtrim(zans.subclass,2)
   specinspect.z = zans.z

   ;----------
   ; Mark objects with CLASS = 'SKY' or 'UNKNOWN'

   qsky = (zans.zwarning AND sdss_flagval('ZWARNING','SKY')) NE 0
   isky = where(qsky, nsky)
   if (nsky GT 0) then specinspect[isky].class = 'SKY'

   ibad = where(zans.zwarning NE 0 AND qsky EQ 0, nbad)
   if (nbad GT 0) then specinspect[ibad].class = 'UNKNOWN'

   ;----------
   ; Set the Yanny header

   hdr = ['plate ' + strtrim(string(plate)), $
          'mjd ' + strtrim(string(mjd)), $
          'vers2d ' + sxpar(zhdr,'VERS2D'), $
          'verscomb ' + sxpar(zhdr,'VERSCOMB'), $
          'vers1d ' + sxpar(zhdr,'VERS1D'), $
          'inspector ' + thisname ]

   ;----------
   ; Explicitly set the Yanny structure description in order to
   ; set the string lengths, whose lengths must be explicitly set
   ; for the benefit of the fTCL Yanny reader.

   structs = [ $
    'typedef struct {', $
    ' int PLATE;', $
    ' int MJD;', $
    ' int FIBERID;', $
    ' char CLASS[20];', $
    ' char SUBCLASS[40];', $
    ' float Z;', $
    ' char INSPECTOR[20];', $
    ' char MANUAL_CLASS[20];', $
    ' char MANUAL_SUBCLASS[40];', $
    ' float MANUAL_Z;', $
    ' char MANUAL_COMMENTS[80];', $
    '} SPECINSPECT;', $
    ' ', $
    'typedef struct {', $
    ' int PLATE;', $
    ' int MJD;', $
    ' int FIBERID;', $
    ' char INSPECTOR[20];', $
    ' char MANUAL_TEXT[800];', $
    '} SPECTEXT;', $
    ' ', $
    'typedef struct {', $
    ' int PLATE;', $
    ' int MJD;', $
    ' int FIBERID;', $
    ' char INSPECTOR[20];', $
    ' char MANUAL_CLASS[20];', $
    ' char MANUAL_SUBCLASS[40];', $
    ' float MANUAL_Z;', $
    ' char MANUAL_COMMENTS[80];', $
    '} SPECBLEND;']

   return
end
;------------------------------------------------------------------------------
