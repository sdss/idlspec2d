; get the multi-gauss parameters for a set of objects

function get_gauss_para, plate,fiber,mjd,fitsdir=fitsdir,errdir=errdir,fname=fname,$
    line=line

if ~keyword_set(fitsdir) then $
fitsdir='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/fits/'

if ~keyword_set(errdir) then $
errdir='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/err/'

if ~keyword_set(fname) then fname=string(plate, format='(i4.4)')+'-'+string(mjd,format='(i5.5)')+ $
  '-'+string(fiber,format='(i4.4)')+'.fits'

fitsfile=fitsdir+fname
para=mrdfits(fitsfile,1,/silent)
errfile=errdir+fname
err=mrdfits(errfile,1,/silent)

linename = strtrim(para.linename)
fitflag = para.fitflag & line_fit = para.line_fit
conti_redchi2=para.conti_redchi2
ind=where( strmatch(strupcase(linename), strupcase(line+'*') ) eq 1)
pp = line_fit[ind]
line_redchi2=(para.line_redchi2)[ind[0]]

; now get error
pp_err=(err.line_fit)[ind,*]

; fix the missing MC trials
n_mc=n_elements(pp_err[0,*])
if n_mc lt 50 then pp_err=[ [pp_err], [pp_err[*,0:(50 - n_mc - 1)]] ]

output={plate:plate,fiber:fiber,mjd:mjd,zfit:para.z,line:line,conti_redchi2:conti_redchi2, line_redchi2:line_redchi2, $
  gauss_pp:pp, gauss_pp_MC:pp_err, fitflag:fitflag[ind[0]]}

return, output

end


; routine to compile multi-gauss parameters for all RM epochs
pro compile_gauss_para, platelist, mjdlist,line=line

if ~keyword_set(line) then line='CIV_BR'


if n_elements(platelist) eq 0 then begin
  file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  platelist=(target[0].plate)[0:31]
  mjdlist=(target[0].mjd)[0:31]
  fiberlist=(target.fiberid)[0:31, 0:848]
endif

if platelist[0] eq 0 then topdir='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/' $
else topdir='/data3/yshen/spectro/bossredux/v5_7_1/'

nnn=n_elements(platelist)

for i=0L, nnn - 1 do begin

   plate=platelist[i] & mjd=mjdlist[i]
   pstr=string(plate, format='(i4.4)') & mstr=string(mjd, format='(i5.5)')
   tt=temporary(result) ; clear this array
   for j=0L, 848 do begin

      if plate ne 0 then begin
         fiber=fiberlist[i, j]
         fitsdir=topdir+pstr+'/qsofit/fits/'
         errdir=topdir+pstr+'/qsofit/err/'
      endif else begin ; for the coadded plate
         fiber=j+1
         fitsdir=topdir+'fits/'
         errdir=topdir+'err/'
      endelse
      fname=string(plate, format='(i4.4)')+'-'+string(mjd,format='(i5.5)')+ $
         '-'+string(fiber,format='(i4.4)')+'.fits'

      ;print, errdir+fname
      ;message, 'stop'

      if file_test(errdir+fname) eq 1 and file_test(fitsdir+fname) eq 1 then begin
         struct=get_gauss_para(plate,fiber,mjd,fitsdir=fitsdir, errdir=errdir,line=line)
         struct=struct_addtags({RMID:j}, struct)
         if n_elements(result) eq 0 then result=struct else result=[result, struct]
      endif

      splog, 'Finished Plate RMID: ', plate, " ", j
   endfor
   if plate ne 0 then outfile=topdir+string(plate, format='(i4.4)')+'/qsofit/'+line+'_par_'+pstr+'-'+mstr+'.fits' $
     else outfile=topdir+line+'_par_'+pstr+'-'+mstr+'.fits'
   mwrfits, result, outfile, /create

endfor

end
