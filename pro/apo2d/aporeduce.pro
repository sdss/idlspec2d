pro aporeduce, filename, indir=indir, outdir=outdir, $
 plugmapfile=plugmapfile, plugdir=plugdir, minexptime=minexptime

   if (n_params() LT 1) then begin
      print, 'Syntax - aporeduce, filename, [ indir=, outdir=, $'
      print, ' plugmapfile=, plugdir=, minexptime= ]'
      return
   endif

   if (NOT keyword_set(indir)) then indir='./'
   if (NOT keyword_set(plugdir)) then plugdir='./'
   if (NOT keyword_set(outdir)) then outdir=indir
   if (NOT keyword_set(minexptime)) then minexptime = 61

   if (size(filename, /tname) NE 'STRING') then begin
      message, 'FILENAME is not a string', /cont
      return
   endif

   filer = strmid(filename,0,3)  ; root 'sdR'
   filec = strmid(filename,4,2)  ; camera number
   filee = strmid(filename,7,8)  ; exposure number

   camnames = ['b1','b2','r1','r2']

   cam = (where(filec EQ camnames))[0]

   if (filer NE 'sdR' OR cam EQ -1) then begin
      message, 'Cannot parse FILENAME '+filename, /cont
      return
   endif

   fullname = (findfile(filepath(filename, root_dir=indir), count=ct))[0]
   if (ct NE 1) then begin
      message, 'Found '+string(nfiles)+' instead of 1', /cont
      return
   endif

   ;----------
   ; Find flavor, plate and MJD

   hdr = headfits(fullname)

   if (size(hdr,/tname) NE 'STRING') then return

   flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
   if (flavor EQ 'target') then flavor = 'science'

   plate = sxpar(hdr, 'PLATEID')
   platestr = string(plate, format='(i4.4)')
   mjd = sxpar(hdr, 'MJD')

   logfile = filepath('logfile-' + strtrim(string(mjd),2) + '.fits', $
    root_dir=outdir)
    
   if (NOT keyword_set(plugmapfile)) then $
       plugmapfile = 'plPlugMapM-'+sxpar(hdr,'NAME')+'.par'
   fullplugmapfile = filepath(plugmapfile, root_dir=plugdir)

   plugexist = keyword_set( findfile(fullplugmapfile) )

   outflat = filepath('tset-'+platestr+'-'+filec+'.fits', $
                           root_dir=outdir)
   flatexist = keyword_set( findfile(outflat) )
   outarc = filepath('wset-'+platestr+'-'+filec+'.fits', $
                           root_dir=outdir)
   arcexist = keyword_set( findfile(outarc) )

   ;----------
   ; Reduce file depending on its flavor: flat, arc, or science

   rstruct = 0
   case flavor of
      'flat' : begin
         if (NOT flatexist AND plugexist) then $
          rstruct = quicktrace(fullname, outflat, fullplugmapfile)
      end
      'arc' : begin
         if (flatexist AND (NOT arcexist)) then $
          rstruct = quickwave(fullname, outflat, outarc)
      end
      'science': begin
          exptime = sxpar(hdr, 'EXPTIME')
          outsci = filepath('sci-'+platestr+'-'+filec+'-'+filee+'.fits',$
                 root_dir=outdir)

          if (flatexist AND arcexist AND exptime GT minexptime) then $
           rstruct = quickextract(outflat, outarc, fullname, outsci)

          ; Now we should check to see if all 4 images have been reduced
;          all4 = findfile(filepath('sci-'+platestr+'-??-'+filee+'.fits',$
;                root_dir=outdir))
;          if (n_elements(all4) EQ 4) then message, 'should check S/N', /cont
       end
       else : begin
          splog, 'Unknown flavor: ', flavor
       end
   endcase

   if (keyword_set(rstruct)) then begin
      ; Append to binary FITS log file a structure with info for this frame
      rstruct = create_struct('MJD', mjd, $
                              'PLATE', plate, $
                              'EXPNUM', filee, $
                              rstruct )
      mwrfits, rstruct, logfile
   endif

end
