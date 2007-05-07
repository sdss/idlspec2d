;+
; NAME:
;   trial_lrg_cat
; PURPOSE:
;   make trial LRG catalog for APO-LSS testing purposes
; CALLING SEQUENCE:
;   trial_lrg_cat
; COMMENTS:
;   uses model mags, not cmodel
; REVISION HISTORY:
;   2006-Dec-18  Blanton, NYU
;-
;------------------------------------------------------------------------------
pro trial_lrg_cat, rerun=rerun

if(NOT keyword_set(rerun)) then rerun=['137', '161']
outfile=getenv('DATA')+'/sdss/trial_lrg_2/trial_lrg_cat.fits'

runs=sdss_runlist(rerun=rerun)
primary=sdss_flagval('RESOLVE_STATUS', 'SURVEY_PRIMARY')

first=1
for irun=0L, n_elements(runs)-1L do begin
    for camcol=1L, 6L do begin
        splog, strtrim(string(runs[irun].run),2)+'/'+strtrim(string(camcol),2)
        objs=sweep_readobj(runs[irun].run, camcol, rerun=runs[irun].rerun, $
                           type='gal')
        if(n_tags(objs) gt 0) then begin
            fibermag = 22.5 - 2.5*alog10(objs.fiberflux>0.01)
            ilist = boss_galaxy_select(objs)
            ilrg = where(ilist GT 0 AND fibermag[3,*] GT 16 AND $
                         (objs.resolve_status AND primary) gt 0, nlrg)
            if(nlrg gt 0) then begin
                if(first eq 1) then begin
                    mwrfits, objs[ilrg], outfile, /create
                    first=0
                endif else begin
                    mwrfits_chunks, objs[ilrg], outfile, /append
                endelse
                
            endif
        endif
    endfor
endfor
    
return
end

