
pro aporeduce, filename, indir, outdir=outdir, plugmapfile=plugmapfile, $
     plugdir=plugdir

     if (NOT keyword_set(indir)) then indir='./'
     if (NOT keyword_set(plugdir)) then plugdir='./'
     if (NOT keyword_set(outdir)) then outdir=indir

     if (size(filename,/tname) NE 'STRING') then begin
        message, 'filename is not a string', /cont
        return
     endif

     filer = strmid(filename,0,3)  ; root 'sdR'
     filec = strmid(filename,4,2)  ; camera number
     filee = strmid(filename,7,8)  ; exposure number

     camnames = ['b1','b2','r1','r2']

     cam = (where(filec EQ camnames))[0]

     if (filer NE 'sdR' OR cam EQ -1) then begin
        message, 'Cannot parse filename '+filename, /cont
        return
     endif

     
     file = findfile(filepath(filename,root_dir=indir), count=nfiles)

     if (nfiles NE 1) then begin
        message, 'Found '+string(nfiles)+' instead of 1', /cont
        return
     endif

     file = file[0]

     ; find flavor, plate and MJD

     hdr = headfits(file)

     if (size(hdr,/tname) NE 'STRING') then return

     flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
     if (flavor EQ 'target') then flavor = 'science'

     plate = sxpar(hdr, 'PLATEID')
     platestr = string(plate, format='(i4.4)')
     mjd = sxpar(hdr, 'MJD')
    
     if (NOT keyword_set(plugmapfile)) then $
         plugmapfile = 'plPlugMapM-'+sxpar(hdr,'NAME')+'.par'
     fullplugmapfile = filepath(plugmapfile, root_dir=plugdir)

     plugmapexist = findfile(fullplugmapfile) NE ''

     flatname = filepath('tset-'+platestr+'-'+filec+'.fits', $
                             root_dir=outdir)
     flatexist = findfile(flatname) NE ''
     arcname = filepath('wset-'+platestr+'-'+filec+'.fits', $
                             root_dir=outdir)
     arcexist = findfile(arcname) NE ''

     case flavor of
       'flat' : if (NOT flatexist[0] AND plugmapexist[0]) then $
            quicktrace, file, flatname, fullplugmapfile

       'arc'  : if (flatexist[0] AND (NOT arcexist[0])) then $
                   quickwave, file, arcname, flatname

       'science': begin
	  exptime = sxpar(hdr, 'EXPTIME')
          outname = filepath('sci-'+platestr+'-'+filec+'-'+filee+'.fits',$
                root_dir=outdir)

          if (flatexist[0] AND arcexist[0] AND exptime GT 61) then $
             quickextract, flatname, arcname, file, outname

          ; Now we should check to see if all 4 images have been reduced
         
          all4 = findfile(filepath('sci-'+platestr+'-??-'+filee+'.fits',$
                root_dir=outdir))
          if (n_elements(all4) EQ 4) then message, 'should check S/N', /cont
             
        end

        else : return
     endcase

end
