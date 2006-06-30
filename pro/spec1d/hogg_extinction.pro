;+
; COMMENTS:
;  - Check out Burki et al, A&ASS, 112, 383 (1995).
; BUGS:
;  - Many things hard-coded.
;  - No comment header.
;-
pro hogg_extinction
prefix= 'hogg_extinction'
camname= ['b1','b2','r1','r2']
ncam= 4

outfilename= prefix+'.fits'
if (NOT file_test(outfilename)) then begin
    if (NOT keyword_set(efficiency)) then begin
        savefile= 'plot_thru.ss'
        if (NOT file_test(savefile)) then plot_thru
        splog, 'restoring '+savefile
        restore, savefile
    endif
    k0= dblarr(n_elements(loglam))
    k0_invvar= k0
    camamp= dblarr(n_elements(loglam),ncam)
    camamp_invvar= camamp

    for ii=0L,n_elements(loglam)-1L do begin
        if ((ii MOD 100) EQ 0) then splog, 1D1^(loglam[ii])
        vgood= where((efficiency[ii,*] GT 0.0) AND $
                     (airmass GT 1.0) AND $
                     (airmass LT 1.5),nvgood)
        if (nvgood GT 100) then begin
            lne= reform(alog(efficiency[ii,vgood]),nvgood)
            aa= [[airmass[vgood]]]
            thiscamlist= [-1]
            for cc=0,ncam-1 do begin
                thiscam= double(strmid(explist[vgood],0,2) EQ camname[cc])
                if (total(thiscam) GT 0.0) then begin
                    thiscamlist= [thiscamlist,cc]
                    aa= [[aa],[thiscam]]
                endif
            endfor
            thiscamlist= thiscamlist[1:n_elements(thiscamlist)-1]
            ww= replicate(1d0,nvgood)
            hogg_iter_linfit, transpose(aa),lne,ww,xx,covar=covar
;                 plot, airmass[good[vgood]],lne,psym=1
;                 oplot, !X.CRANGE,transpose([[1,1],[!X.CRANGE]]##xx),psym=0
            k0[ii]= xx[0]
            k0_invvar[ii]= 1.0/covar[0,0]
            for cc=0,n_elements(thiscamlist)-1 do begin
                camamp[ii,thiscamlist[cc]]= xx[1+cc]
                camamp_invvar[ii,thiscamlist[cc]]= 1.0/covar[1+cc,1+cc]
            endfor
        endif
    endfor
    splog, 'writing file '+outfilename
    mwrfits, [[loglam],[k0],[k0_invvar],[camamp],[camamp_invvar]], $
      outfilename,/create
endif

set_plot, 'ps'
device, filename= prefix+'.ps'
hogg_plot_defaults
readcol, '~/Longslit/calib/extinction/atm_trans_am1.0.dat', $
  longwave,longthru
longwave= longwave*1D4
readcol, '~/primus/data/atmosphere.dat', $
  primwave,primthru
filename= prefix+'.fits'
foo= mrdfits(filename)
good= where((foo[*,2] GT 0.0))
hoggwave= 1D1^foo[good,0]
; hoggthru= 0.10+0.58*foo[good,1]
plot, hoggwave,hoggthru,psym=10,/nodata, $
  xrange=[3500,9500],xtitle= 'wavelength  (A)', $
  yrange=[-1.5,0.1],ytitle= 'd ln(throughput) / d airmass', $
  title='all cameras'
oplot, longwave,alog(longthru),color=djs_icolor('grey'),psym=10
oplot, primwave,-primthru,color=djs_icolor('grey'),psym=10
oplot, hoggwave,hoggthru,psym=10
; oplot, 1D1^foo[good,0],foo[good,1]+2.0/sqrt(foo[good,2]),psym=10
; oplot, 1D1^foo[good,0],foo[good,1]-2.0/sqrt(foo[good,2]),psym=10
for cc=0,ncam-1 do begin
    plot, 1D1^foo[good,0],foo[good,cc+3],psym=10, $
      xrange=[3500,9500],xtitle= 'wavelength  (A)', $
      yrange=[-5,0],ytitle= 'ln(throughput)', $
      title=camname[cc]
;    oplot, 1D1^foo[good,0],foo[good,3]+2.0/sqrt(foo[good,4]),psym=10
;    oplot, 1D1^foo[good,0],foo[good,3]-2.0/sqrt(foo[good,4]),psym=10
endfor
device,/close

return
end
