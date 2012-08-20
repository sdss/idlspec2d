pro bosslog

   ;----------
   ; Find all the spFrame files, and sort by exposure number

   spawn, 'ls -d ' + getenv('BOSS_SPECTRO_DATA') + '/?????', mjdlist
   
   for imjd=0L, n_elements(mjdlist)-1 do begin
      splog, 'Searching for files in ' + mjdlist[imjd]
      flist = findfile( $
       filepath('sdR-??-????????.fit*', root_dir=mjdlist[imjd]), count=ct)
      if (ct GT 0) then begin
         if (NOT keyword_set(files)) then files = flist $
          else files = [files, flist]
      endif
   endfor

   files = files[ sort(fileandpath(files)) ]
   nfile = n_elements(files)

   ;----------
   ; Construct a structure with all the header info
   ; Use a specific file as a template

   file_template = files[(where(strmatch(files,'*sdR-b1-00061798.fit*')))[0]]
   hdr = headfits(file_template)
   struct1 = create_struct('SIMPLE', 'T')
   hname = 'SIMPLE'
 
 ;Tag names to extract from the template file
 cards = ['decdeg', 'exposure', 'exptime', 'ff', 'ffs', 'flavor', $
          'hgcd', 'mjd', 'ne', 'plateid', 'ra', 'radeg', 'cameras', $
          'filename', 'tai-beg', 'dec', 'cartid', 'name', 'mapid', 'hartmann']
 cards = strupcase(cards)
          
   for j=1L, n_elements(hdr)-1L do begin
   hname1 = strtrim(strmid(hdr[j],0,8),2)
   junk = where(cards eq hname1, count)
   if (count ne 0) then begin   
     thisname = repstr(hname1,'-','_')
     if (thisname EQ 'NE') then thisname = 'NEON'
     if (hname1 NE ' ' AND hname1 NE 'END' ) then begin
       value = (sxpar(hdr,thisname))[0]
       struct1 = create_struct(struct1, thisname, value)
       hname = [hname, hname1] ; Header names
     endif
   endif
   endfor

   structs = replicate(struct1, nfile)
   tags = tag_names(struct1) ; Tag names

   for ifile=0, nfile-1L do begin
       splog, 'Reading file ', ifile, ' of ', nfile, $
         ': ', fileandpath(files[ifile])
         
        catch, theError
        if theError ne 0 then begin
          print, 'Something wrong happened'
          continue
        endif

        hdr = headfits(files[ifile]) 
        
        ; Get the mjd, camera, and exposure from the filename itself
        junk = strsplit(files[ifile], '/', length=l)
        n = n_elements(junk)
        mjd = strmid(files[ifile], junk[n-2], l[n-2])
        filename = strmid(files[ifile], junk[n-1], l[n-1])
        junk = strsplit(filename, '-', length=l)
        camera = strmid(filename, junk[1], l[1])
        exposure = strmid(filename, junk[2], 8)
                        
        for j=0L, n_elements(hname)-1L do begin
          value = sxpar(hdr, hname[j])
          if (hname[j] eq 'MJD') then value = mjd
          if (hname[j] eq 'CAMERAS') then value = camera
          if (hname[j] eq 'EXPOSURE') then value = exposure
          structs[ifile].(j) = value 
        endfor
   endfor
   
   catch, /cancel
   
 mwrfits, structs, 'bosstest.fits' ,/create

 return 
 end
   
   
   
   
   
   