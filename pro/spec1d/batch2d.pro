; This doesn't yet check to see if things have already been run for some plates
pro batch2d, topindir, _EXTRA=EXTRA

   if (NOT keyword_set(topindir)) then topindir='2d_test'

   ;----------
   ; Create list of plate numbers

   platelist = get_mjd_dir(topindir, _EXTRA=EXTRA)
   fullplatelist = concat_dir(topindir, platelist)
   nplate = n_elements(platelist)

   ;----------
   ; Generate the IDL script files

   fq = "'"
   fullscriptfile = strarr(nplate)
   for iplate=0, nplate-1 do begin
      i = rstrpos(platelist[iplate], '.')
      if (i EQ -1) then i = strlen(platelist[iplate])
      fullscriptfile[iplate] = djs_filepath(platelist[iplate]+'.batch', $
       root_dir=fullplatelist[iplate])
      openw, olun, fullscriptfile[iplate], /get_lun
      printf, olun, 'cd, ' + fq+fullplatelist[iplate]+fq
      printf, olun, 'spreduce2d'
      printf, olun, 'spcombine'
      printf, olun, 'exit'
      close, olun
      free_lun, olun
   endfor

   ;----------
   ; Create list of input files

   infile = transpose( [ [fullscriptfile] ] )

   ;----------
   ; Create list of expected output files

   outfile = infile ; ???

   ;----------
   ; Prioritize to do the most recent plates first

   priority = lonarr(nplate)
   isort = sort(platelist)
   priority[isort] = lindgen(nplate) + 1

   ;----------
   ; Determine which computers to use for these reductions

   hostfile = filepath('batch2d.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   yanny_read, hostfile, pp
   hostconfig = *pp[0]
   yanny_free, pp

   ;----------
   ; Begin the batch jobs

hostconfig[*].remotedir = '/scr0/data'
   command = 'nice +19 idl ' + fullscriptfile
   batch, infile, outfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   return
end

