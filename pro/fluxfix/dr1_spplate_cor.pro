;-------------------------------------------------------------------------------
;+
; NAME:
;   dr1_spplate_cor
;
; PURPOSE:
;   Modify a list of spPlate files to improve the spectrophotometry 
;
; CALLING SEQUENCE:
;   dr1_spplate_cor, dr1_platelist, logfile =, spectro_data_dir = 
;   
; INPUTS:
;   dr1_platelist    - Name of yanny parameter file listing the plate and mjd
;                      numbers to modify
; 
; OPTIONAL KEYWORDS:
;   logfile          - Name of log file
;   spectro_data_dir - Directory path to original spPlate files.  If not set
;                      then it is assumed that the environment variable
;                      'SPECTRO_DATA' is set.
;
; OUTPUTS:
;   A modified spPlate-$PLATE-$MJD.fits file is written to the current directory
;
; PROCEDURES CALLED:
;   splog
;   spplate_correct
;   yanny_read
;
; REVISION HISTORY:
;   Written by C. Tremonti, 1 Nov 2002, JHU
;-
;--------------------------------------------------------------------------------
pro dr1_spplate_cor, platelist_par, logfile = logfile, $
    spectro_data_dir = spectro_data_dir

  splog, filename = logfile

  ;-----------------------------------------------------------------------------
  ; Read back yanny parameter file to get list of plates and mjds
  ;-----------------------------------------------------------------------------

  yanny_read, platelist_par, pdata

  plate = (*pdata).plate
  mjd = (*pdata).mjd

  ;-----------------------------------------------------------------------------
  ; Determine path name to original spPlate files  
  ;-----------------------------------------------------------------------------

  platestr = string(plate, format='(I4.4)')
  mjdstr = string(mjd, format='(I5)')

  if (NOT keyword_set(spectro_data_dir)) then $
  spectro_data_dir = getenv('SPECTRO_DATA')
  if (NOT keyword_set(spectro_data_dir)) then $
      message, 'Environment variable SPECTRO_DATA must be set!'

  platename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'

  ;-----------------------------------------------------------------------------
  ; Loop through plates
  ;-----------------------------------------------------------------------------

  for indx = 0, n_elements(plate) - 1 do begin
 
    orig_spplate = filepath(platename[indx], root_dir=spectro_data_dir, $
                 subdirectory=platestr[indx])

    ;------------
    ; If input file does not exist then skip it

    if file_test(orig_spplate) ne 1 then begin
       splog, 'Skipping ' + platename[indx]
       splog, 'No File found in ' + spectro_data_dir + platestr[indx]
       continue
    endif
 
    ;------------
    ; If file already exists in the output directory then skip it

    if file_test(platename[indx]) eq 1 then begin
       splog, 'Output file already exists. Skipping ' + platename[indx]
       continue
    endif else splog, 'Modifying ' + platename[indx]
    
    ;------------
    ; Modify the plate

    spplate_correct, plate[indx], mjd[indx], spectro_data_dir = $
      spectro_data_dir, /remove_smear, /zeropoint_cor

  endfor

end
