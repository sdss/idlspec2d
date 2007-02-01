;+
; NAME:
;   elodie_filelist
;
; PURPOSE:
;   Get a listing of all valid Elodie filenames.
;
; CALLING SEQUENCE:
;   res = elodie_filelist(minwave=)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   res      - list of full filenames for all unmasked Elodie spectra
;
; OPTIONAL OUTPUTS:
;   minwave  - if the filter.par file exists and contains a minwave
;              keyword, its value. Else returns 1.0
;
; COMMENTS:
;   If the filter.par file does not exist, all files in the
;   $ELODIE_DIR/LL_ELODIE/ directory are returned, and minwave is set
;   to 0
;
; EXAMPLES:
;
; BUGS:
;   The masking logic should be inverted: only listed files should be matched.
;
; PROCEDURES CALLED:
;   yanny_read()
;-
;------------------------------------------------------------------------------
function elodie_filelist, minwave=minwave

                                ; List all Elodie path- and filenames
   elodie_path = getenv('ELODIE_DIR')
   allfiles = findfile(filepath('0*', root_dir=elodie_path, $
    subdir='LL_ELODIE'), count=nstar)
   lastslash = strpos(allfiles[0], '/', /reverse_search)
   fparts = strmid(allfiles, lastslash+1)

                                ; --------
                                ; Newer versions of the Elodie product
                                ; contain a .par file filtering the
                                ; catalog spectra. For older versions,
                                ; just return a listing of all the files.
   minwave = 1.0
   filtfile = findfile(filepath('filter.par', root_dir=elodie_path, subdir='etc'))
   if filtfile eq '' then return, allfiles

   yanny_read, filtfile, filt, hdr=hdr, stnames=stnames, /anonymous
   especlist = *filt
   yanny_free, pdata
   emasked = where(especlist.state eq 'MASKED')

   minwave = yanny_par(hdr, 'minwave')
   if minwave eq '' then minwave = 1.0
   minwave = float(minwave)
                                ; I'm sure there's an IDL idiom for
                                ; this...
   maskedspecs = intarr(n_elements(fparts))
   for i=0, n_elements(emasked)-1 do begin
      maskedname = especlist[i].filename
      match = where(strpos(fparts, maskedname) eq 0, cnt)
      if cnt EQ 0 then $
         splog, "WARNING: masked elodie file " + maskedname + " not found" $
      else $
         maskedspecs[match] = 1
   end

   return, allfiles[where(maskedspecs eq 0)]
end
