;+
; NAME:
;   plate_spec_image
; PURPOSE:
;   Create spectroscopic images for a plate
; CALLING SEQUENCE:
;   plate_spec_image, plate, mjd=, run2d=, run1d=, topdir=
; INPUTS
;  plate - plate #
;  mjd - MJD of observation
; OPTIONAL INPUTS:
;  run1d - 1d rerun number ('' by default)
;  run2d - 2d rerun number ('' by default)
;  topdir - directory to use instead of $BOSS_SPECTRO_REDUX
;  xra - [2] wavelength range to plot (Ang)
; OPTIONAL KEYWORDS:
;  /noclobber - do not clobber existing image PNG files
;               (will clobber HTML)
; COMMENTS:
;  Creates files in:
;     $BOSS_SPECTRO_REDUX/images/[run2d]/[run1d]/[plate4]-[mjd]
;  with the naming convention:
;     spec-image-[plate4]-[mjd]-[fiber4].png
;     spec-image-[plate4]-[mjd]-[fiber4].thumb.png
;  Also creates an index.html file which is just a table of
;  the thumbnails.
; REVISION HISTORY:
;   15-May-2009  by M. Blanton, NYU
; VERSION:
;   $Id: plate_spec_image.pro 64936 2015-09-21 18:50:32Z julianbautista $
;------------------------------------------------------------------------------
pro plate_spec_image, plate, mjd=mjd, run2d=run2d, run1d=run1d, $
                      topdir=topdir,  xra=xra, silent=silent, $
                      noclobber=noclobber, legacy=legacy, plates=plates

 RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR
 CPU, TPOOL_NTHREADS = 1

if(NOT keyword_set(run2d)) then run2d=''
if(NOT keyword_set(run1d)) then run1d=''
if(NOT keyword_set(topdir)) then topdir= getenv('BOSS_SPECTRO_REDUX')


get_field_type, fieldid=plate, mjd=mjd, legacy=legacy, plates=plates, fps=fps

pmjd=field_to_string(plate)+'-'+string(mjd, f='(i5.5)')
outdir=topdir+'/images/'+run2d+'/'+run1d+'/'+pmjd
file_mkdir, outdir
readspec, plate, mjd=mjd, zans=zans, run2d=run2d, run1d=run1d, $
  topdir=topdir,legacy=legacy,plates=plates, plug=plug

if(n_tags(zans) eq 0) then begin
   splog, 'No 1d reductions for this plate.'
   return
endif

if file_test(outdir+'/index.html') then rmfile, outdir+'/index.html'
openw, unit, outdir+'/tmp-index.html', /get_lun
printf, unit, '<?xml version="1.0" encoding="UTF-8"?>'
printf, unit, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
; printf, unit, '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">'
printf, unit, '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'
printf, unit, '<head>'
printf, unit, '<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'
printf, unit, '<title>'+pmjd+'</title>'
printf, unit, '<style type="text/css">'
printf, unit, 'body { background: #111; }'
printf, unit, 'td { color: gray; }'
printf, unit, '</style>'
; printf, unit, '<link href="layout.css" rel="stylesheet" type="text/css"></link>'
printf, unit, '</head>'
printf, unit, '<body>'
printf, unit, '<table border="0" cellspacing="5">'
nper=5L
for i=0L, n_elements(zans)-1L do begin
    if((i mod nper) eq 0) then $
      printf, unit, '<tr>'
    fiber=zans[i].target_index
    ra=zans[i].fiber_ra
    dec=zans[i].fiber_dec
    if keyword_set(legacy) then begin
       pmjdf= pmjd+'-'+string(f='(i4.4)', fiber)
    endif else begin
       if keyword_set(plates) then begin
         catalogid=plug[i].catalogid
         if catalogid eq 0 then begin
            pmjdf= pmjd+'-'+strtrim(fiber,2)
         endif else begin
             pmjdf= pmjd+'-'+strtrim(string(catalogid),2)
         endelse
       endif else begin
         catalogid=plug[i].catalogid
         pmjdf=pmjd+'-'+strtrim(string(catalogid),2)
       endelse
    endelse
    stamp='http://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg?TaskName=Skyserver.Chart.Image&ra='+strtrim(string(f='(f40.5)', ra),2)+'&dec='+strtrim(string(f='(f40.5)', dec),2)+'&scale=0.1&width=512&height=512&opt=G&query=&Grid=on'
    currbase='spec-image-'+pmjdf
    outbase=outdir+'/'+currbase
    msg='plotting target_indx:'+strtrim(fiber,2)
    if not keyword_set(legacy) then msg+=' catalogid:'+strtrim(catalogid,2)
    msg+= ' ('+strtrim(i+1,2)+'/'+strtrim(n_elements(zans),2)+')'
    splog, msg
    sdss_spec_image, outbase, plate, fiber, mjd=mjd, $
      run2d=run2d, run1d=run1d, topdir=topdir, xra=xra, silent=silent, $
      noclobber=noclobber, legacy=legacy, plates=plates
    printf, unit, '<td><a href="'+stamp+'">'+pmjdf+ $
      '</a><br /><a href="'+currbase+'.png">'
    printf, unit, '<img src="'+currbase+'.thumb.png" alt="'+pmjdf+'" /></a>'
    printf, unit, '</td>'
    if(((i mod nper) eq nper-1L) OR (i eq n_elements(zans)-1L)) then $
      printf, unit, '</tr>'
endfor
printf, unit, '</table>'
printf, unit, '</body>'
printf, unit, '</html>'
close, unit
free_lun, unit
file_move, outdir+'/tmp-index.html', outdir+'/index.html', /overwrite, /verbose
return
end
;------------------------------------------------------------------------------
