; Script to run SOS as observer@sos3.apo
; The plPlugMap files must be in your current directory.
; Log files are written in your current directory to indicate which
;   files have already been reduced.
; One can use the 4 CPUs by starting one instance of this proc
;   for each CAMNAME; default to CAMNAME='??' to do all cameras with 1 CPU.
; If MJD is not specified, then use the most recent.
; Loop forever, and roll over to new MJDs (if MJD not specified).
;------------------------------------------------------------------------------
pro apod1, mjd=mjd, camname=camname

   if (keyword_set(mjd)) then begin
      mjdstr = strtrim(string(mjd[0]),2)
   endif else begin
      spawn, '\ls -d /data/spectro/5????', alldir
      alldir = fileandpath(alldir)
      foo = max(long(alldir), ilast)
      mjdstr = alldir[ilast]
   endelse
   if (NOT keyword_set(camname)) then camname = '??'

   indir='/data/spectro/'+mjdstr
   outdir='/data/boss/sos/'+mjdstr
   plugdir='./'
   copydir='/data/boss/sos/combined'
;indir='/Users/schlegel/rawdata/'+mjdstr
;outdir='./'
;plugdir='/Users/schlegel/rawdata/plugmaps'
;copydir=''
   filename = fileandpath(findfile(indir+'/sdR-'+camname+'*.fit.gz',count=nfile))
   filename = filename[sort(filename)]
   for i=0L,nfile-1L do begin
      logfile=filename[i]+'.log'
      thisfile = findfile(logfile)
      if (keyword_set(thisfile) EQ 0) then begin
         aporeduce, filename[i], indir=indir, outdir=outdir, $
          plugdir=plugdir, copydir=copydir
         splog, filename=logfile, 'Done '+filename[i], /close
      endif
   endfor

   return
end
;------------------------------------------------------------------------------
pro apod, mjd=mjd, camname=camname

   while (1) do begin
      apod1, mjd=mjd, camname=camname
      wait, 10
   endwhile

   return
end
;------------------------------------------------------------------------------
