; compile the rms profile (from prepspec) and avg profile (from qsofit)
; for a given broad line
; use the results to study the rms/avg profile of broad lines

pro make_rmsavg_line, rmid, line=line, result=result

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_final.fits'
tt = mrdfits(file,1)

if ~keyword_set(line) then line = 'Hb'
if n_elements(rmid) eq 0 then rmid = 0
zsys = tt[rmid].zsys

if line eq 'Hb' then begin
  lam0=4862.68D
  maxv = 1d4 ; km/s
  line_qsofit = 'Hbeta_br'
  line_prepspec = 'hb'
endif
if line eq 'MgII' then begin
  lam0=2798.75D
  maxv = 1d4 ; km/s
  line_qsofit = 'MgII_br'
  line_prepspec = 'mg2'
endif
if line eq 'CIV' then begin
  lam0=1549.06D
  maxv = 1d4 ; km/s
  line_qsofit = 'CIV_br'
  line_prepspec = 'c4'
endif


nobj = n_elements(rmid)

; +- maxv km/s from the center
npix = maxv/3d5/(1d-4*alog(10.))
loglam = dindgen(npix*2 + 1)*1d-4 + alog10(lam0) - npix*1d-4
nwave = n_elements(loglam)
wave = 10.D^loglam

result = {rmid:-1L,zsys:-1.D,line:line,lam0:lam0, wave:wave,loglam:loglam,rms:dblarr(nwave),rms_err:dblarr(nwave), $
   avg:dblarr(nwave), rms_avg:dblarr(nwave)-1.D }
result = replicate(result, nobj)
result.rmid=rmid & result.zsys=zsys

; currently using 2014a prepspec results, but should switch to 2014b later
prepspec_dir = '/data3/yshen/ftp/sdssrm/collab/prepspec/2014b/'
qsofit_dir = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/civ_3gauss/fits/'
for i=0, nobj - 1 do begin

  ; rms line-only profile from prepspec
  tag = 'rm' + string(rmid[i],format='(i3.3)')
  prepspec_file = prepspec_dir + tag +  '/' + tag + '_' + line_prepspec + '_w.dat'

  if file_test(prepspec_file) eq 1 then begin
    readcol, prepspec_file,format='d,d,d',lam,flux,err, /silent
    if n_elements(lam) gt 10 then begin 
      result[i].rms = interpol(flux, alog10(lam),loglam )
      result[i].rms_err = interpol(err, alog10(lam),loglam )

      ; avg line-only profile from multi-Gaussian model from qsofit
      fits_file = qsofit_dir + '0000-56837-' + string(rmid[i]+1,format='(i4.4)')+'.fits'
      fits = mrdfits(fits_file,1,/silent)
      linefit = fits.line_fit & linename = strtrim(fits.linename)
      ind = where(linename eq line_qsofit)
      lnlam1 = alog(wave*(1. + zsys[i])/(1. + fits.z))
      result[i].avg = manygauss(lnlam1, linefit[ind])

      result[i].rms_avg = result[i].rms/result[i].avg
    endif
  endif
endfor

end
