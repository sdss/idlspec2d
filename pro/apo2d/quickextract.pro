pro quickextract, flatname, arcname, file, outname

  if (N_PARAMS() LT 4) then begin
    print, 'Syntax - quickextract, flatname, arcname, file, outname, plugmapfile'
  endif

  if (n_elements(flatname) NE 1) then return 
  if (n_elements(arcname) NE 1) then return 

  nobj  = n_elements(file)
  nout  = n_elements(outname)

  if (nout NE nobj) then begin
       print, "Number of files OUT is different then number of files IN"
       return
  endif

  for i=0,nout - 1 do begin

     sdssproc, file[i], image, hdr=hdr, color=color

     camname = strtrim(sxpar(hdr, 'CAMERAS'),2)

     tset = mrdfits(flatname,1)
     plugsort = mrdfits(flatname,2)
     fibermask = mrdfits(flatname,3)
     wset = mrdfits(arcname,1)

     traceset2xy, tset, ycen, xcen

     ; boxcar extract
     flux = extract_boxcar(image, xcen, radius = 3.0)

     ; estimate fluxivar
     fluxivar = 1.0/(abs(flux) + 10.0)

     mwrfits, flux, outname[i], /create
     mwrfits, fluxivar, outname[i]

   endfor
end

     

