
pro write_spflat, flatinfoname, iflat, flatstruct, flathdr, $
                  arcname, nbright, ymodel, scatter, $
                  nowrite=nowrite, writeflatmodel=writeflatmodel, outdir=outdir
                  

    sxaddpar, flathdr, 'NBRIGHT', nbright, $
          'Number of bright pixels (>10^5) in extracted flat-field'

    sxaddpar, flathdr, 'FLATNAME', flatstruct[iflat].name, $
            'Name of flat'
    sxaddpar, flathdr, 'ARCNAME', arcname[flatstruct[iflat].iarc], $
            'Name of associated arc'
    sxaddpar, flathdr, 'PROFTYPE', flatstruct[iflat].PROFTYPE, $
            'Extract Profile (1:guass 2:ExpCubic 3:DoubleGauss 4:2+3)'
            
    quad = ['LL','LR','UL','UR']
    foreach mw,flatstruct[iflat].medwidth, idx do $
            sxaddpar, flathdr, 'MEDWIDT'+strtrim(idx,2), mw, $
                    ' Median spatial dispersion widths in '+quad[idx]+' quadrant'



    if not keyword_set(nowrite) then begin
        if keyword_set(outdir) then begin
            file_mkdir, outdir
            flatinfofile = filepath(string(format='(a,i8.8,a)',flatinfoname, $
                                               sxpar(flathdr, 'EXPOSURE'),'.fits'),$
                                root_dir=outdir)
        endif else begin
            flatinfofile = string(format='(a,i8.8,a)',flatinfoname, $
                sxpar(flathdr, 'EXPOSURE'), '.fits')
        endelse
        mwrfits_named, *flatstruct[iflat].fflat, flatinfofile, hdr = flathdr, name='FFLAT', /create
        mwrfits_named, *flatstruct[iflat].tset, flatinfofile, name='TSET'
        mwrfits_named, *flatstruct[iflat].fibermask, flatinfofile, name='FIBERMASK'
        mwrfits_named, *flatstruct[iflat].widthset, flatinfofile, name='WIDTHSET'
        mwrfits_named, *flatstruct[iflat].superflatset, flatinfofile, name='SUPERFLATSET'
        mwrfits_named, *flatstruct[iflat].xsol, flatinfofile, name='XSOL'


        spawn, ['gzip', '-f', flatinfofile], /noshell
        
        ; ASB: write flat image model info if requested:
        if keyword_set(writeflatmodel) then begin
            flatmodelfile = string(format='(a,i8.8,a)',flatinfoname + $
                'MODELIMG-', sxpar(flathdr, 'EXPOSURE'), '.fits')
            mwrfits_named, flatimg, flatmodelfile, name = 'FLATIMG'/create
            mwrfits_named, flativar, flatmodelfile, name = 'FLATIVAR'
            mwrfits_named, ymodel + scatter,flatmodelfile, name = 'ymodel + scatter'
            spawn, ['gzip', '-f', flatmodelfile], /noshell
        endif
    endif
    
    flatstruct[iflat].hdr = ptr_new(flathdr)

    return
end
