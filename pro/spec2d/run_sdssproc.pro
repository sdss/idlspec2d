;+
;
; NAME:
;  run_sdssproc
;
; PURPOSE:
;  run SDSSPROC over many raw sdR files and write image and invvar files.
;
; USAGE:
;  run_sdssproc [, indir=indir, outdir=outdir, mjdlist=mjdlist, $
;   clobber=clobber, gzip=gzip, minflat=minflat, maxflat=maxflat]
;
; ARGUMENTS:
;  indir: raw data directory, defaulting to $RAWDATA_DIR
;  outdir: top directory to write processed files, defaulting
;          to working dir ('.')
;  mjdlist: list of MJDs to process, defaulting to all MJDs
;           found under indir
;  clobber: set this to foce clobbering of existing image and invvar files.
;  gzip: set this to cause output files to be gzipped (default is no gzip)
;  minflat: sdssproc "minflat" keyword, default to 0.8
;  maxflat: sdssproc "maxflat" keyword, default to 1.2
;
; WRITTEN:
;  bolton@utah 2011jan
;
;-

pro run_sdssproc, indir=indir, outdir=outdir, mjdlist=mjdlist, clobber=clobber, $
 gzip=gzip, minflat=minflat, maxflat=maxflat

; Defaults and prep:
if (n_elements(minflat) ne 1) then minflat = 0.8
if (n_elements(maxflat) ne 1) then maxflat = 1.2
if (not keyword_set(indir)) then indir = getenv('RAWDATA_DIR')
if (not keyword_set(outdir)) then outdir = '.'
if (not keyword_set(mjdlist)) then mjdlist = fileandpath(file_search(indir+'/?????'))
   mjdlist = fileandpath(mjdlist)

s_mjdlist = strtrim(string(mjdlist), 2)
nmjd = n_elements(s_mjdlist)

if (nmjd lt 1) then begin
   print, 'No MJDs specified or found.'
   return
endif

; Loop over MJDs:
for i = 0L, nmjd - 1 do begin
   print, 'MJD = ' + s_mjdlist[i] + ':'
   out_this = outdir + '/' + s_mjdlist[i]
   if (file_test(out_this) eq 0) then spawn, 'mkdir ' + out_this
   sdr_full = file_search(indir + '/' + s_mjdlist[i] + '/' + 'sdR-[b,r][1,2]-????????.fit*')
   nsdr = n_elements(sdr_full)
   print, ' Found ' + strtrim(string(nsdr),2) + ' sdR files.'
   for j = 0L, nsdr - 1 do begin
      sdr_file = fileandpath(sdr_full[j])
      sdr_base = (strsplit(sdr_file, '.', /extract))[0]
      print, ' File ' + strtrim(string(j),2) + ' (' + sdr_base + ')'
      outfile = out_this + '/' + sdr_base + '-IMAGE.fits'
      varfile = out_this + '/' + sdr_base + '-INVVAR.fits'
      if (keyword_set(clobber) or (file_test(outfile+'*') eq 0) or (file_test(varfile+'*') eq 0)) then begin
         spawn, 'rm -f ' + outfile+'*'
         spawn, 'rm -f ' + varfile+'*'
         sdssproc, sdr_full[j], outfile=outfile, varfile=varfile, /applycrosstalk, $
          /applypixflat, /applybias, minflat=minflat, maxflat=maxflat, /silent
         if keyword_set(gzip) then begin
            spawn, 'gzip ' + outfile
            spawn, 'gzip ' + varfile
         endif
      endif
   endfor
endfor

return
end
