; NAME:
;   spcalib_qa
;
; PURPOSE:
;   Compare photometric accuracy of standards
;
; CALLING SEQUENCE:
;   SpCalib_QA, [run2d=, fieldid=, mjd=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   field       - field to include
;   mjd         - MJD to include
;   run2d       - RUN2D version of reduction
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Depends on the spAll files (either full run2d version or field-mjd version)
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; Function Called:
;   mpfitfun
;   field_to_string
;   djs_filepath
;   mrdfits
;   sdss_flagval
;
; External PROCEDURES CALLED:
;   plot
;   XYOUTS
;   cpbackup
;
; Internal PROCEDURES CALLED:
;   std_hist
;
; REVISION HISTORY:
;   21-June-2022  Written by S. Morrison (UIUC)
;------------------------------------------------------------------------------

pro std_hist, fratio, bs=bs, xmin=xmin, xmax=xmax, filt=filt

    if not keyword_set(bs) then bs=0.015
    if not keyword_set(xmin) then xmin=-.3
    if not keyword_set(xmax) then xmax=-.2
    if not keyword_set(filt) then filt='g'

    logf = alog10(fratio)
    yhist=histogram(logf,locations=xhist,binsize=bs,omin=hmin,omax=hmax,min=xmin, max=xmax)
    p=mpfitfun('gauss1', xhist, yhist, yhist,weights=1.D, [0,0.2,1],/QUIET)
    xhist=[hmin-bs, xhist,hmax+bs/2.0]
    yhist=[0, yhist, 0]
    plot, xhist, yhist, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[0,max(yhist)+2], xticks=5, xminor=10, FONT=0, thick=2.75, psym=10, xtitle='log (synflux/calibflux) ['+filt+']', ytitle='Nstd'
    XYOUTS, xmin+0.02,(max(yhist)*.9), 'mean, sigma!C'+string(p[0],format='(F5.2)')+','+string(p[1],format='(F5.3)'),  CHARSIZE=1.2,  FONT=0
    xfit = findgen((xmax-xmin)/0.001, increment=0.001,start=xmin)
    oplot, xfit, GAUSS1(xfit, p), psym=0, linestyle=1
    return
end

pro SpCalib_QA, run2d=run2d, fieldid=fieldid, mjd=mjd
    RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR

    if not keyword_set(run2d) then run2d = getenv('RUN2D')
    outname='spCalib_QA-'+run2d
    if keyword_set(fieldid) then begin
        spallfile='spAll-'+field_to_string(fieldid)+'-'+strtrim(mjd,2)+'.fits'
        spallfile = djs_filepath(spallfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/spectra/full/'+field_to_string(fieldid)+'/'+strtrim(mjd,2))
        outname = outname+'-'+field_to_string(fieldid)+'-'+strtrim(mjd,2)
        outname=djs_filepath(outname+'.ps', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/'+field_to_string(fieldid))
        cpbackup, outname
    endif else begin
        spallfile = 'spAll-'+run2d+'.fits'
        spallfile = djs_filepath(spallfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d)
        outname = djs_filepath(outname+'.ps', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d)
    endelse

    spall=mrdfits(spallfile,1)
    ind = where(strmatch(spall.objtype, 'SPECTROPHOTO_STD',/fold_case) $
            and ((spall.zwarning AND sdss_flagval('ZWARNING', 'UNPLUGGED')) eq 0)) ; std after removing unplugged ones
    f_ratio = spall[ind].SPECTROSYNFLUX/spall[ind].calibflux
    ind_good = where(f_ratio[1,*] gt 0)

    mydevice = !D.NAME
    SET_PLOT, 'ps'
    DEVICE, FILENAME=outname,/LANDSCAPE, /times

    !p.multi=[0,1,3]
    !psym=10
    !y.margin=[4,4]
    !x.ticklen=0.075
    !X.thick=2.0
    !Y.thick=2.0
    !X.charsize=2
    !Y.charsize=2


    std_hist, f_ratio[1, ind_good], bs=0.015, xmin=-.3, xmax=.2, filt='g'
    std_hist, f_ratio[2, ind_good], bs=0.015, xmin=-.3, xmax=.2, filt='r'
    std_hist, f_ratio[3, ind_good], bs=0.015, xmin=-.3, xmax=.2, filt='i'



    if keyword_set(fieldid) then XYOUTS, 0.5,0.98, 'Field='+field_to_string(fieldid)+' MJD='+strtrim(mjd,2),  CHARSIZE=1.2,  FONT=0, alignment=0.5,/NORMAL

    DEVICE, /CLOSE
; Return plotting to the original device:
    SET_PLOT, mydevice
end
