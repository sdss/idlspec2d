; compare the synthetic mag from spectrum and
; mags from imaging

pro comp_syn_photo_mags, rmid, prepspec=prepspec, calibdir=calibdir, $
  photodir=photodir,outdir=outdir

; output figfile
if ~keyword_set(outdir) then outdir='/data3/quasar/yshen/work/lags/prepspec/output/' 
if ~keyword_Set(calibdir) then calibdir='wh_skysub/' ; specdata
if ~keyword_Set(photodir) then photodir='/data3/quasar/yshen/work/lags/prelim_photo/' ; imaging data
if n_elements(prepspec) then prepspec=1 ; default is to fix the flux scaling using prepspec results

; find the plate-fiber-mjd for each epoch
common rm_fibermap, target_info, epoch_info
if n_elements(target_info) eq 0 then begin
  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target_info=mrdfits(file,1)
  epoch_info=mrdfits(file,2)
endif
plate=target_info[rmid].plate
fiber=target_info[rmid].fiberid
mjd=target_info[rmid].mjd
mean_mjd=epoch_info.mean_mjd
nep=n_elements(plate)

; first read in the spectrum
synmag=dblarr(5, nep) & synmag_err=dblarr(5,nep)
rm_readspec,plate,fiber,mjd=mjd,calibdir=calibdir,wave=wave,$
   flux=flux,invvar=ivar,/silent
if keyword_set(prepspec) then begin ; applying the flux correction from prepspec
    ; this is a constant scaling applied to the entire spectrum
  pt=get_prepspec_pt(rmid)
  for i=0L, nep-1 do begin
     flux[*,i] = flux[*,i]/pt[i]
     ivar[*,i] = ivar[*,i]*(pt[i])^2
  endfor
endif

; now compute synthetic mag and mag_err
for i=0,nep-1 do begin
  mag=rm_spec2mag(wave[*,i],flux[*,i],ivar[*,i],synmag_err=magerr)
  synmag[*,i]=mag & synmag_err[*,i]=magerr
endfor

; now get photometric LC
rmtag='rm'+string(rmid,format='(i3.3)')
lcfile=photodir+'rm'+string(rmid,format='(i3.3)')+'photlcg.txt'
readcol, lcfile, format='d,d,d,a',mjd_g,mag_g,err_g,tag_g
lcfile=photodir+'rm'+string(rmid,format='(i3.3)')+'photlci.txt'
readcol, lcfile, format='d,d,d,a',mjd_i,mag_i,err_i,tag_i

; now make a plot
outfile=outdir+'rm'+string(rmid,format='(i3.3)')+'_lc_check.ps'
begplot, name=outfile, /color,/landscape
ploterror, mean_mjd, synmag[1,*],synmag_err[1,*], title=rmtag+',g',xrange=[56650.,56850.], $
   xtickformat='(i0)', psym=symcat(6)
oploterror,mjd_g,mag_g,err_g,psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red')
  
  ; i band
ploterror, mean_mjd, synmag[3,*],synmag_err[3,*], title=rmtag+',i',xrange=[56650.,56850.], $
   xtickformat='(i0)', psym=symcat(6)
oploterror,mjd_i,mag_i,err_i,psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red')

endplot
cgfixps, outfile

end
