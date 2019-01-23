; obtain the flux scaling factor from prepspec
; obsolete, use rm_get_pt instead

function get_prepspec_pt, rmid, datadir=datadir,pterr=pterr

if ~keyword_set(datadir) then datadir='/data3/quasar/yshen/work/lags/prepspec/'

; find the file
rmtag='rm'+string(rmid,format='(i3.3)')
ptfile=datadir+rmtag+'/'+rmtag+'_p0_t.dat'
if file_test(ptfile) eq 0 then download_prepspec_fits,suffix='p0_t.dat',subdir='/ACF/',rmid=rmid

; now check ptfile again
;print, ptfile
if file_test(ptfile) eq 1 then begin
  readcol,ptfile,format='x,d,d', pt, pterr,/silent
endif else begin
  splog, 'no ptfile available, use original spectrum:'
  pt=replicate(1., 32)
endelse

return, pt
end
