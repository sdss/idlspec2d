
pro aporeduce, exposure, indir, outdir=outdir, plugmapfile=plugmapfile

     if (NOT keyword_set(outdir)) then outdir=indir

     expstr = string(exposure,format='(i8.8)')
     files = findfile(filepath('sdR*'+expstr $
                        +'.fit',root_dir=indir), count=nfiles)

     if (nfiles NE 4) then return
     files = files[sort(files)]

     camnames = ['b1','b2','r1','r2']

     ; find flavor, plate and MJD

     hdr = headfits(files[0])

     if (size(hdr,/tname) NE 'STRING') then return

     flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
     if (flavor EQ 'target') then flavor = 'science'

     plate = sxpar(hdr, 'PLATEID')
     platestr = string(plate, format='(i4.4)')
     mjd = sxpar(hdr, 'MJD')
    
     if (NOT keyword_set(plugmapfile)) then $
         plugmapfile = '/data/spectro/plugmap/plPlugMapM-'+$
             sxpar(hdr,'NAME')+'.par'

     for icam=0, 3 do begin

        flatname = filepath('tset-'+platestr+'-'+camnames[icam]+'.fits', $
                             root_dir=outdir)
        flatexist = findfile(flatname) NE ''
        arcname = filepath('wset-'+platestr+'-'+camnames[icam]+'.fits', $
                             root_dir=outdir)
        arcexist = findfile(arcname) NE ''

        case flavor of
          'flat' : begin 
                if (NOT flatexist) then $
                   quicktrace, files[icam], flatname
                end

          'arc'  : begin
                if (flatexist AND (NOT arcexist)) then $
                   quickwave, files[icam], arcname, flatname
                end

          'science': begin
	        exptime = sxpar(hdr, 'EXPTIME')
                outname = filepath('sci-'+platestr+'-'+camnames[icam]+$
                             '-'+expstr+'.fits', root_dir=outdir)

                if (flatexist AND arcexist AND exptime GT 300) then begin
                    quickextract, flatname, arcname, files[i], $
                     outname, plugmapfile
                endif
                end

           else : return
        endcase

      endfor

end
