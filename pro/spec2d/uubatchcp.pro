;+
; NAME:
;   uubatchcp
;
; PURPOSE:
;   To copy run2d and run1d reductions from scratch to boss_spectro_redux
;   For an individual CPU cycle called by uubatchpbs
;   To function on clusters without node sharing (e.g. University of Utah)
;
; CALLING SEQUENCE:
;   uubatchcp, [ planfile, topdir=, run2d=, run1d=, scratchdir=, /zcode, /galaxy, /skip2d]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of pipeline plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   topdir     - Optional override value for the environment
;                variable $BOSS_SPECTRO_REDUX.
;   run2d      - Optional override value for the environment variable $RUN2D
;   run1d      - Optional override value for the environment variable $RUN1D
;   scratchdir   - If set, then treat this as topdir until the computation is complete
;   zcode      - If set, cp reductions produced by Zcode in auto mode.
;   galaxy     - If set, cp reductions produced by Galaxy (Portsmouth, PCA) Suite of Products.
;   skip2d     - If set, then skip the Spectro-2D reductions.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;   01-May-2011  Adapted for uubatchpbs by Joel R. Brownstein, University of Utah, 
;                to allow reductions to be first written to $BOSS_SCRATCH_DIR
;                primarily to prevent uubatchpbs from crashing in case 
;                $BOSS_SPECTRO_REDUX is full, and also to select which subset of the 
;                pipeline products are moved into BOSS_SPECTRO_REDUX for retention.
;-
;------------------------------------------------------------------------------
pro uubatchcp_init, plan

    common uubatchcp_data, has_platemjd, platemjd, plate, mjd
    
    ;----------
    ; Find the SPEXP structure

    allseq = yanny_readone(plan, 'SPEXP', hdr=hdr, /anon)
    has_platemjd = (N_elements(allseq) GT 0)
    if (has_platemjd) then begin

        ;----------
        ; Find keywords from the header

        mjd = string(yanny_par(hdr, 'MJD'), format='(i05.5)')
        plate = string(yanny_par(hdr,'plateid'), format='(i04.4)')
        platemjd = '-' + plate + '-' + mjd
        
    endif
      
end
;------------------------------------------------------------------------------
pro uubatchcp, planfile, topdir=topdir, run2d=run2d, run1d=run1d, scratchdir=scratchdir, $
 zcode=zcode, galaxy=galaxy, skip2d=skip2d

    common uubatchcp_data

   ;----------
   ; Determine the top-level of the output directory tree

   if (not keyword_set(topdir)) then topdir = getenv('BOSS_SPECTRO_REDUX')
   if strpos(topdir,'/',strlen(topdir)-1) lt 0 then topdir+='/'

   if (not keyword_set(scratchdir)) then scratchdir = getenv('BOSS_SCRATCH_DIR')
   if strpos(scratchdir,'/',strlen(scratchdir)-1) lt 0 then scratchdir+='/'

   if (keyword_set(scratchdir)) then begin
       print, "UUBATCHCP:  Please set scratch_dir keyword or $BOSS_SCRATCH_DIR"
       return
   endif
   
   if (topdir eq scratchdir) then begin
       print, "UUBATCHCP:  Nothing to do because topdir=scratchdir"
       return
   endif
   
   if (keyword_set(run2d)) then run2d = strtrim(run2d,2) $
    else run2d = getenv('RUN2D')
   if (keyword_set(run1d)) then run1d = strtrim(run1d,2) $
    else run1d = getenv('RUN1D')
    
    if (NOT keyword_set(planfile)) then begin
        planfiles = djs_filepath('spPlan2d-*.par', root_dir=topdir, subdir=run2d+'/*')
        planfile = findfile(planfiles)
    endif
        
    ;----------
    ; If multiple plan files exist, then call this script recursively
    ; for each such plan file.

    if (N_elements(planfile) GT 1) then begin
        for i=0, N_elements(planfile)-1 do $
        uubatchcp, planfile[i], topdir=topdir, run2d=run2d, run1d=run1d, scratchdir=scratchdir, $
 zcode=zcode, galaxy=galaxy, skip2d=skip2d
        return
    endif

    uubatchcp_init, planfile[0]
    
    if (has_platemjd) then begin
    
        plate2d = plate + '/' + run2d

        topdir2d = djs_filepath('', root_dir=topdir, subdir=plate2d)
        scratchdir2d = djs_filepath('', root_dir=scratchdir, subdir=plate2d)

        if (not file_test(scratchdir2d)) then begin
           print, "UUBATCHCP: Nothing to copy. "+scratchdir2d+" does not exist."
           return
        endif
       
        if (not file_test(topdir2d)) then file_mkdir, topdir2d
        
        cameras = ['r1','r2','b1','b2']
       
        if (not keyword_set(skip2d)) then begin
            scratchdir2d_files = djs_filepath('*'+platemjd+'.*', root_dir=scratchdir2d)
            file_copy, scratchdir2d_files, topdir2d, /overwrite, /require_directory
            for c=0,n_elements(cameras)-1 do begin
                camera_files = djs_filepath('*'+cameras[c]+'*.*', root_dir=scratchdir2d)
                wh = where(file_test(djs_filepath(file_basename(camera_files), root_dir=topdir2d)) eq 0,nc)
                if (nc gt 0) then file_copy, camera_files[wh], topdir2d, /require_directory
            endfor
        endif
        
        topdir1d = djs_filepath('', root_dir=topdir2d, subdir=run1d)
        scratchdir1d = djs_filepath('', root_dir=scratchdir2d, subdir=run1d)

        if (not file_test(scratchdir1d)) then begin
           print, "UUBATCHCP: Nothing to copy. "+scratchdir1d+" does not exist."
           return
        endif

        if (not file_test(topdir1d)) then file_mkdir, topdir1d

        scratchdir1d_files = djs_filepath('*'+platemjd+'.*', root_dir=scratchdir1d)
        file_copy, scratchdir1d_files, topdir1d, /overwrite, /require_directory

       
    endif
   
end