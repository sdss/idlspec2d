pro batch1d, fullplatefile

   ;----------
   ; Create list of plate files

   if (NOT keyword_set(fullplatefile)) then $
    fullplatefile = findfile('*/spPlate-*.fits')
   platefile = fileandpath(fullplatefile, path=localpath)
   nplate = n_elements(platefile)

   ;----------
   ; Find which programs are already done by looking at the log files

   platemjd = strmid(platefile, 8, 10)
   diaglog = localpath + '/spDiag1d-' + platemjd + '.log'

   qdone = lonarr(nplate)
   endstring = 'Successful completion of SPREDUCE1D'
   for iplate=0, nplate-1 do begin
      spawn, 'tail -1 ' + diaglog[iplate], tailstring
      qdone[iplate] = strpos(tailstring[0], endstring) NE -1
      if (qdone[iplate]) then $
       print, 'File ' + platefile[iplate] + ' already reduced'
   endfor

   indx = where(qdone EQ 0, nplate)
   if (nplate EQ 0) then begin
      print, 'All plates have been reduced already'
      return
   endif

   fullplatefile = fullplatefile[indx]
   platefile = platefile[indx]
   localpath = localpath[indx]
   platemjd = platemjd[indx]
   diaglog = diaglog[indx]

   ;----------
   ; Generate the IDL script files

   fq = "'"
   fullscriptfile = strarr(nplate)
   for iplate=0, nplate-1 do begin
      i = rstrpos(platefile[iplate], '.')
      if (i EQ -1) then i = strlen(platefile[iplate])
      fullscriptfile[iplate] = $
       djs_filepath(strmid(platefile[iplate],0,i)+'.batch', $
       root_dir=localpath[iplate])
      openw, olun, fullscriptfile[iplate], /get_lun
      printf, olun, 'cd, ' + fq+localpath[iplate]+fq
      printf, olun, 'spreduce1d, ' + fq+platefile[iplate]+fq
      printf, olun, 'exit'
      close, olun
      free_lun, olun
   endfor

   ;----------
   ; Create list of input files

   infile = transpose( [ [fullscriptfile], [fullplatefile] ] )

   ;----------
   ; Create list of expected output files

   zallfile = localpath + '/spZall-' + platemjd + '.fits'
   zbestfile = localpath + '/spZbest-' + platemjd + '.fits'
   diagps = localpath + '/spDiag1d-' + platemjd + '.ps'
   outfile = transpose( [ [diaglog], [diagps], [zallfile], [zbestfile] ] )

   ;----------
   ; Prioritize to do the most recent plates first

   priority = lonarr(nplate)
   isort = sort(platemjd)
   priority[isort] = lindgen(nplate) + 1

   ;----------
   ; Determine which computers to use for these reductions

   hostfile = filepath('batch1d.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   yanny_read, hostfile, pp
   hostconfig = *pp[0]
   yanny_free, pp

   ;----------
   ; Begin the batch jobs

   command = 'nice +19 idl ' + fullscriptfile
   batch, infile, outfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   return
end

