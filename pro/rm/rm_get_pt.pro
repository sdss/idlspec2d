; Get the PrepSpec p0(t) corrections
; note that p0(t) is supposed to be divided by the raw spectrum

function rm_get_pt, rmid, topdir=topdir, err=err, nep=nep, mjd_all=mjd_all

;if ~keyword_set(nep) then nep=32L
nep = n_elements(mjd_all)

if ~keyword_set(topdir) then $
  topdir='/data3/yshen/ftp/sdssrm/collab/prepspec/ACBFJ/'
rmtag = 'rm'+string(rmid,format='(i3.3)')

ptfile = topdir + rmtag + '/' + rmtag+'_p0_t.dat'
pt=replicate(1., nep) & err=0.*pt
if file_test(ptfile) eq 1 then begin
  readcol,ptfile,format='d,d,d',mjd,lnpt, lnpt_err, /silent
  if n_elements(mjd) gt 0 then begin
    mjd = floor(mjd)
    nep1=n_elements(mjd)
    for i=0L, nep1 - 1 do begin
       ind=where( abs(floor(mjd_all) - (mjd[i]+50000.) ) le 1d-2)
       if ind[0] ne -1 then begin
         pt[ind]=exp(lnpt[i]) & err[ind]=pt[ind]*lnpt_err[i]
       endif
    endfor
  endif
endif

return, pt

end
