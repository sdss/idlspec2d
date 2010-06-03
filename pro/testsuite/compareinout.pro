;Compares the input and output of spreduce2d/spcombine_v5 to test
;whether the output is "good enough"

function compareinout, sdss2=sdss2, sdss3=sdss3, epsilon=epsilon

if not(keyword_set(epsilon)) then epsilon = 0.05

if keyword_set(sdss3) then begin
    outputflux =  mrdfits('$IDLSPEC2D_DIR/pro/testdata/sdss3/0001/spPlate-0001-55025.fits', 0, h)
stop
    m = sxpar(h, 'CD1_1')
    b = sxpar(h, 'CRVAL1')
    so = size(outputflux)
    outputlambda = 10.^(m*findgen(so[1]) + b)

    inputflux =   mrdfits('$IDLSPEC2D_DIR/pro/testdata/sdss3/0001/spPlate-0001-55025-input.fits', 0, h)
    inputlambda = mrdfits('$IDLSPEC2D_DIR/pro/testdata/sdss3/0001/spPlate-0001-55025-input.fits', 1, h)
    objtype = mrdfits('$IDLSPEC2D_DIR/pro/testdata/sdss3/0001/spPlate-0001-55025-input.fits', 2, h)
    sky = where(objtype eq 'SKY')
    inputflux[*, sky] = 0.
    resampledinputflux = fltarr(so[1], so[2])
endif

if keyword_set(sdss2) then begin
    outputflux = mrdfits('$IDLSPEC2D_DIR/pro/testdata/sdss2/0400/spPlate-0400-51820.fits', 0, h)
    m = sxpar(h, 'CD1_1')
    b = sxpar(h, 'CRVAL1')
    so = size(outputflux)
    outputlambda = 10.^(m*findgen(so[1]) + b)

    inputflux = mrdfits('$IDLSPEC2D_DIR/pro/testdata/sdss2/0400/spPlate-0400-51820-input.fits', 0, h)
    m = sxpar(h, 'CD1_1')
    b = sxpar(h, 'CRVAL1')
    si = size(inputflux)
    inputlambda = 10.^(m*findgen(si[1]) + b)
endif

quant = 2
minoutputlambda = min(outputlambda) 
maxoutputlambda = max(outputlambda)


for i = 0, so[2]-1 do begin
    nowlambda      = where(inputlambda[*, i] ge minoutputlambda and inputlambda[*, i] le maxoutputlambda)
    nowinputlambda = [inputlambda[min(nowlambda)-1, i], inputlambda[nowlambda, i], $
                      inputlambda[max(nowlambda)+1, i]]
    nowinput       = [inputflux[min(nowlambda)-1, i], inputflux[nowlambda, i], inputflux[max(nowlambda)+1, i]]
    minw = 0
    maxw = n_elements(nowinputlambda)-1 
    for j = 0, so[1]-1 do begin 
        rr = abs(lambda[j] - nowinputlambda)
        wplc = where(rr eq min(rr))
        wplc = wplc[0]
            if wplc-quant le 0 then begin
                wnow = nowinputlambda[minw:wplc+quant]
                qenow = nowinput[minw:wplc+quant]
            endif
            if wplc+quant ge n_elements(nowinputlambda) then begin
                wnow = nowinputlambda[wplc-quant:maxw]
                qenow = nowinput[wplc-quant:maxw]
            endif
            if wplc-quant gt 0 and wplc+quant lt n_elements(nowinputlambda) then begin
                wnow  = nowinputlambda[wplc-quant:wplc + quant]
                qenow = nowinput[wplc-quant:wplc + quant]
            endif
            xy2traceset, wnow, qenow, set, ncoeff=2.*quant-1, yfit=yfit, /silent
            traceset2xy, set, lambda[j], result
            resampledinput[j, i] = result
    endfor
endfor

si = size(resampledinput)
chisq = fltarr(so[2])

for i = 0, so[2]-1 do chisq[i] = max((resampledinput - outputflux)/resampledinput)

bad = where(chisq gt epsilon, countbad)
good = where(chisq le epsilon, countgood)

if countbad le 4 then begin
    !p.multi = [2, 2]
    for i = 0, countbad- 1 do begin
        plot, outputlambda, outputflux[*, bad[i]], $
          title = 'Failed Reduction', xtitle = 'Wavelength [A]', ytitle = 'Flux electrons'
        oplot, inputlambda[*, bad[i]], inputflux[*, bad[i]], color=100
    endfor
    if countgood gt 0 then begin
        plot, outputlambda, outputflux[*, good[0]], $
          title = 'Successful Reduction', xtitle = 'Wavelength [A]', ytitle = 'Flux electrons'
        oplot, inputlambda[*, good[0]], inputflux[*, good[0]], color=100
    endif
endif


if countbad gt threshold then return, 0 else return, 1
stop

end
