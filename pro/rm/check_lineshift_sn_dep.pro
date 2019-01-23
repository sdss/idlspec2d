; check the SN dependence of lineshift measurements


pro check_lineshift_sn_dep, dwave_arr=dwave_arr, sn=sn, linename=linename

file='/data3/yshen/work/lineshifts/lineshift.fits'
result=mrdfits(file,1)
tag0=tag_names(result)

file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target=mrdfits(file,1)

plate=target[0:848].plate
fiber=target[0:848].fiberid
mjd=target[0:848].mjd
med_sn=target[0:848].med_sn

; now for each object, get the difference between single-epoch and coadded epoch
if ~keyword_set(linename) then linename='CIV'
topdir='/data3/yshen/spectro/bossredux/v5_7_1/'
for i=0, 31L do begin

  plate1=plate[i,0] & mjd1=mjd[i,0]
  pstr=string(plate1, format='(i4.4)')
  mstr=string(mjd1,format='(i5.5)')

  fitsfile=topdir + pstr+'/qsofit/qso_prop-'+pstr+'-'+mstr+'.fits'
  fit=mrdfits(fitsfile, 1)

  tag1=tag_names(fit)
  ind0=where(tag0 eq linename or tag0 eq linename+'_BR')
  ind_err0=where(tag0 eq linename+'_ERR' or tag0 eq linename+'_BR_ERR')

  ind1=where(tag1 eq linename or tag1 eq linename+'_BR')
  ind_err1=where(tag1 eq linename+'_ERR' or tag1 eq linename+'_BR_ERR')

  ind=where( (result.(ind_err0))[0,*] gt 0 and (fit.(ind_err1))[0,*] gt 0, nnn )

  if nnn gt 0 then begin

     diff_wav= (fit.(ind1))[0,ind] - (result.(ind0))[0, ind]
     err_wav = sqrt( ((result.(ind_err0))[0,ind])^2 + ((fit.(ind_err1))[0,ind])^2  )
     if n_elements(dwave_arr) eq 0 then begin
         dwave_arr=reform(diff_wav/err_wav) 
         sn=reform(med_sn[i, ind])
     endif else begin
         dwave_arr=[dwave_arr, reform(diff_wav/err_wav)]
         sn=[sn, reform(med_sn[i, ind])]
     endelse

  endif 

endfor


end
