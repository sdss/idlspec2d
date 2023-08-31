
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
        sxaddpar, flathdr, 'EXTNAME', 'fflat'
        mwrfits, *flatstruct[iflat].fflat, flatinfofile, flathdr, /create
            
        sxaddpar, hdr, 'EXTNAME', 'tset'
        mwrfits, *flatstruct[iflat].tset, flatinfofile, hdr
            
        sxaddpar, hdr, 'EXTNAME', 'fibermask'
        mwrfits, *flatstruct[iflat].fibermask, flatinfofile, hdr
        
        sxaddpar, hdr, 'EXTNAME', 'widthset'
        mwrfits, *flatstruct[iflat].widthset, flatinfofile, hdr
            
        sxaddpar, hdr, 'EXTNAME', 'superflatset'
        mwrfits, *flatstruct[iflat].superflatset, flatinfofile, hdr
            
 ;       sxaddpar, hdr, 'EXTNAME', 'qbad'
 ;       mwrfits, flatstruct[iflat].qbad, flatinfofile, hdr
            
        sxaddpar, hdr, 'EXTNAME', 'xsol'
        mwrfits, *flatstruct[iflat].xsol, flatinfofile, hdr
            
 ;       sxaddpar, hdr, 'EXTNAME', 'name'
 ;       mwrfits, flatstruct[iflat].name, flatinfofile, hdr


        spawn, ['gzip', '-f', flatinfofile], /noshell
        
        ; ASB: write flat image model info if requested:
        if keyword_set(writeflatmodel) then begin
            flatmodelfile = string(format='(a,i8.8,a)',flatinfoname + $
                'MODELIMG-', sxpar(flathdr, 'EXPOSURE'), '.fits')
            mwrfits, flatimg, flatmodelfile, /create
            mwrfits, flativar, flatmodelfile
            mwrfits, ymodel + scatter, flatmodelfile
            spawn, ['gzip', '-f', flatmodelfile], /noshell
        endif
    endif
    
    flatstruct[iflat].hdr = ptr_new(flathdr)

    return
end
