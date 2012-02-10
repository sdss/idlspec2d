;+
;
; NAME:
;   bolton_sn_fit
;
; PURPOSE:
;   robust fit to SN^2 versus fibermag
;
; USAGE:
;   sn2ref = bolton_sn_fit(fibermag, sn2vector [, $
;    lims=lims, refmag=refmag, quantile=quantile, doplot=doplot])
;
; ARGUMENTS:
;   fibermag: x coordinate for fit
;   sn2vector: y coordinate, log10'ed before fitting
;
; OPTIONAL ARGUMENTS:
;   lims: 2-element vector giving low-to-high (i.e.,
;    bright-to-faint) fibermag limit over which to fit.
;    (default = [20., 22])
;   refmag: reference fibermag at which to eval final fit
;    (default = 21.)
;   quantile: reject points beyond this quantile in terms
;    of absolute deviation from first-pass fit.
;    (default = 0.8 -- i.e., reject most deviant 20%)
;   doplot: set keyword for simple splot diagnostic.
;
; COMMENTS:
;   Does LADFIT, rejects outliers, does LADFIT again!
;
;  ** CURRENT DEFAULTS ARE FOR BOSS R-CHANNEL!! **
;  ** FOR B-CHANNEL, OVERRIDE DEFAULTS!! **
;
;  So far, only spot-tested, but seems fairly robust.
;
; WRITTEN:
;   Bolton @ Utah 2010 May, with pair-coding by Olmstead.
;
;-

function bolton_sn_fit, fibermag, sn2vector, lims=lims, refmag=refmag, $
 quantile=quantile, doplot=doplot,camera=camera,sn2=sn2


if (not keyword_set(refmag)) then refmag = 21.
if (not keyword_set(quantile)) then quantile = 0.8

if (not keyword_set(camera)) then camera =' '
if (keyword_set(camera)) then begin
    if camera eq 'b1' or camera eq 'b2' then begin 
        lims =[21,23]
        refmag=22
    endif
endif
if (not keyword_set(lims)) then lims = [20., 22]

wh = where((fibermag gt lims[0]) and (fibermag lt lims[1]), nwh)
xloc = fibermag[wh]
yloc = alog10(sn2vector[wh])
lfit = ladfit(xloc, yloc)
yfit1 = poly(xloc, lfit)
resid = abs(yloc - yfit1)
resid_srt = resid[sort(resid)]
cutval = resid_srt[long(quantile * nwh)]
wsub = where(resid lt cutval)
lfit2 = ladfit(xloc[wsub], yloc[wsub])
yfit2 = poly(xloc, lfit2)
sn2_fid = 10.^poly(refmag, lfit2)

if keyword_set(doplot) then begin
    splot, xloc, yloc, ps=6, color=2,ytitle='S/N^2   '+sn2_fid,xtitle='fibermag  '+camera
    soplot, xloc[wsub], yloc[wsub], ps=4, color=1
    soplot, xloc, yfit1, ps=1, color=3
    soplot, xloc, yfit2, ps=1, color=6
endif
;return,camera
if keyword_set(camera) then print,camera
if keyword_set(sn2) then print,sn2
return, sn2_fid 
;if (keyword_set(camera)) then  begin 
;    return, camera
;endif

end

