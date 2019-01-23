;+
; Make QA plots to access the quality of spectral measurements
; from rm_qsofit

pro plot_fits_qa, plate, mjd, range=range, plot_conti=plot_conti, comment=comment,$
    psfile=psfile


    platestr=string(plate,format='(i4.4)')
    mjdstr = string(mjd,format='(i5.5)')

    outfile = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
       + '/' + platestr + '/qsofit/qso_prop-' + platestr+'-'+mjdstr+'.fits'

    result=mrdfits(outfile,1)
    ; Only use a subset of the objects
    if keyword_set(range) then result = result[range[0]:range[1]]

    ; Simulated QSOs used in RM proposal
    file='/data2/quasar/yshen/work/RM/TDSS/level2/sim_qso.fits'
    simqso = mrdfits(file,1)

    target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
    fibermap = mrdfits(target_file,1,/silent)
    imag = reform((fibermap.psfmag)[3,*])

    if keyword_set(psfile) then begplot, name=psfile,/encap,/color,xsize=7,ysize=6,/cmyk

    csize=1.5 & ssize=0.5 & cthick=5

    if keyword_set(plot_conti) then !P.multi=[0,1,2]
    title = platestr + '-' + mjdstr
    if keyword_set(comment) then title = title + ', ' + comment
    plot, simqso.imag, 1d-16/simqso.f_line, psym=3, xtitle=textoidl('i_{PSF}'), ytitle='Line Flux Error / Line Flux' $
    , xrange=[18, 22], yrange = [1d-3, 2d0],/ylog,charsize=csize,title=title,charthick=1.5,$
    ytickname=textoidl(['10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^1']),pos=[0.13,0.13,0.96,0.92]

    ; plot the observed line measurement quality
    linelist = ['HALPHA', 'HBETA', 'MGII', 'CIII_ALL', 'CIV']
    ;linelist = reverse(linelist)
    color=cgcolor(['red','magenta','cyan','green','blue'])
    sym = [4,17,6, 9,16]
    tagname = tag_names(result)

    nline = n_elements(linelist)
    for i=2,nline-1 do begin

       itag = where(tagname eq strupcase(linelist[i]) )
       ind = where( (result.(itag))[2,*] gt 0 )

       imag_line = imag[ind] 
       itag1 = where(tagname eq strupcase(linelist[i] + '_ERR') )
       err_line = ( (result.(itag1))[2,ind] )
    
       oplot, imag_line, err_line*alog(10.D), psym=symcat(sym[i],thick=cthick),color=color[i], $
         symsize=ssize

    endfor

    for i=0,1 do begin

       itag = where(tagname eq strupcase(linelist[i]) )
       ind = where( (result.(itag))[2,*] gt 0 )

       imag_line = imag[ind]
       itag1 = where(tagname eq strupcase(linelist[i] + '_ERR') )
       err_line = ( (result.(itag1))[2,ind] )

       oplot, imag_line, err_line*alog(10.D), psym=symcat(sym[i],thick=cthick+1),color=color[i], $
         symsize=ssize+0.5

    endfor


    oplot, [18,22],[0.1,0.1],thick=cthick,line=2, color=cgcolor('red')

    legend, linelist, box=0,pos=[0.15, 0.89],/norm, color=color,textcolor=color,charsize=csize
    legend, 'Sim_QSO', box=0,pos=[0.32,0.89],/norm, textcolor=cgcolor('opposite'),charsize=csize

    if keyword_Set(plot_conti) then begin
    plot, simqso.imag, 1d-15/simqso.f_conti, psym=3, xtitle='PSF i Mag', ytitle='Conti Error/Flux' $
        , xrange=[18, 22], yrange = [1d-3, 2d0],/ylog,charsize=csize

    contiwave = [1350.,3000.,5100.]
    conti_lum_str = 'logL'+string(contiwave,format='(i0)')
    nconti = n_elements(conti_lum_str)
    for i=0L, nconti-1 do begin
       itag = where(tagname eq strupcase(conti_lum_str[i]) )
       ind = where( (result.(itag)) gt 0)
       imag_line = imag[ind]
       itag1 = where(tagname eq strupcase(conti_lum_str[i]+'_ERR') )
       err_line = ( (result.(itag1))[ind])

       oplot, imag_line, err_line*alog(10.D), psym=symcat(sym[i]),color=color[i],symsize=ssize
    endfor

    oplot, [18,22],[0.1,0.1],thick=cthitk,line=2, color=cgcolor('red')

    legend, conti_lum_str, box=0,pos=[0.15,0.45],/norm, color=color,textcolor=color,charsize=csize
    legend, 'Sim_QSO', box=0,pos=[0.28,0.45],/norm, textcolor=cgcolor('opposite'),charsize=csize

    !P.multi = 0
    endif

    if keyword_set(psfile) then endplot

end
