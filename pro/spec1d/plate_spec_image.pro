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
;   $Id$
;------------------------------------------------------------------------------
pro plate_spec_image, plate, mjd=mjd, run2d=run2d, run1d=run1d, topdir=topdir, $
                      xra=xra, silent=silent

if(NOT keyword_set(run2d)) then run2d=''
if(NOT keyword_set(run1d)) then run1d=''
if(NOT keyword_set(topdir)) then topdir= getenv('BOSS_SPECTRO_REDUX')

pmjd=string(plate, f='(i4.4)')+'-'+string(mjd, f='(i5.5)')
outdir=topdir+'/images/'+run2d+'/'+run1d+'/'+pmjd
; spawn, 'mkdir -p '+outdir
file_mkdir, outdir
readspec, plate, mjd=mjd, zans=zans, run2d=run2d, run1d=run1d, $
  topdir=topdir

if(n_tags(zans) eq 0) then begin
   splog, 'No 1d reductions for this plate.'
   return
endif

rmfile, outdir+'/index.html'
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
    fiber=zans[i].fiberid
    pmjdf= pmjd+'-'+string(f='(i4.4)', fiber)
    currbase='spec-image-'+pmjdf
    outbase=outdir+'/'+currbase
    sdss_spec_image, outbase, plate, fiber, mjd=mjd, $
      run2d=run2d, run1d=run1d, topdir=topdir, xra=xra, silent=silent
    printf, unit, '<td>'+pmjdf+ $
      '<br /><a href="'+currbase+'.png">'
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
; spawn, 'mv '+outdir+'/tmp-index.html '+outdir+'/index.html'
file_move, outdir+'/tmp-index.html', outdir+'/index.html', /overwrite, /verbose
return
end
;------------------------------------------------------------------------------
