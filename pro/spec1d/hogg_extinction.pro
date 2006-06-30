;+
; BUGS:
;  - Many things hard-coded.
;  - No comment header.
;-
pro hogg_extinction
prefix= 'hogg_extinction'

camname= ['b1','b2','r1','r2']
for cc= 0,3 do begin
    outfilename= prefix+'_'+camname[cc]+'.fits'

    if (NOT file_test(outfilename)) then begin
        if (NOT keyword_set(efficiency)) then begin
            savefile= 'plot_thru.ss'
            if (NOT file_test(savefile)) then plot_thru
            splog, 'restoring '+savefile
            restore, savefile
        endif
        lne0= dblarr(n_elements(loglam))
        lne0_invvar= lne0
        k0= lne0
        k0_invvar= lne0

        splog, 'starting work on camera '+camname[cc]
        good= where(strmid(explist,0,2) EQ camname[cc],ngood)
        for ii=0L,n_elements(loglam)-1L do begin
            vgood= where((efficiency[ii,good] GT 0.0) AND $
                         (airmass[good] GT 1.0) AND $
                         (airmass[good] LT 2.0),nvgood)
            if (nvgood GT 100) then $
              if ((max(airmass[good[vgood]])-min(airmass[good[vgood]])) GT 0.2) $
              then begin
                lne= reform(alog(efficiency[ii,good[vgood]]),nvgood)
                aa= transpose([[replicate(1d0,nvgood)],[airmass[good[vgood]]]])
                ww= replicate(1d0,nvgood)
                hogg_iter_linfit, aa,lne,ww,xx,covar=covar
;            plot, airmass[good[vgood]],lne,psym=1
;            oplot, !X.CRANGE,transpose([[1,1],[!X.CRANGE]]##xx),psym=0
                lne0[ii]= xx[0]
                lne0_invvar[ii]= 1.0/covar[0,0]
                k0[ii]= xx[1]
                k0_invvar[ii]= 1.0/covar[1,1]
            endif
        endfor
        splog, 'writing file '+outfilename
        mwrfits, [[loglam],[k0],[k0_invvar],[lne0],[lne0_invvar]], $
          outfilename,/create
    endif
endfor

set_plot, 'ps'
device, filename= prefix+'.ps'
hogg_plot_defaults
for cc=0,3 do begin
    filename= prefix+'_'+camname[cc]+'.fits'
    foo= mrdfits(filename)
    good= where(foo[*,2] GT 0.0)
    plot, 1D1^foo[good,0],foo[good,1],psym=10, $
      xrange=[3500,9500],xtitle= 'wavelength  (A)', $
      yrange=[-1.5,0.1],ytitle= 'd ln(throughput) / d airmass', $
      title=camname[cc]
endfor
device,/close

return
end
