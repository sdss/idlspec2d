; If using other machines in Peyton, set
;   topdir='/peyton/scr/spectro0/data/2d_test'

;------------------------------------------------------------------------------
function batch2d_rawfiles, planfile, outfile=outfile

   nplan = n_elements(planfile)
   if (nplan GT 1) then begin
      infiles = batch2d_rawfiles(planfile[0], outfile=outfile)
      for i=1, nplan-1 do begin
         infiles = [infiles, $
          batch2d_rawfiles(planfile[i], outfile=tmpout)]
         outfile = [outfile, tmpout]
      endfor
      return, infiles
   endif

   yanny_read, planfile, pp, hdr=hdr
   extractdir = yanny_par(hdr, 'extractdir')
   thismjd = yanny_par(hdr, 'MJD')
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')
   for ii=0, n_elements(pp)-1 do begin
      sname = tag_names(*pp[ii], /structure_name)
      if (sname EQ 'SPEXP') then begin
         expfiles = ((*pp[ii]).name)[*]
      endif
   endfor

   mjdstr = string(thismjd, format='(i05.5)')
   infiles = djs_filepath(expfiles, root_dir=mjdstr)

   ;----------
   ; Replace the prefix 'sdR' with 'spFrame' in the science frames
   ; and the suffix '.fit' with '.fits'

   newnames = expfiles
   for i=0, n_elements(newnames)-1 do begin
      jj = strpos(newnames[i], '-')
      kk = rstrpos(newnames[i], '.')
      if (jj NE -1 AND kk NE -1) then $
       newnames[i] = 'spFrame' + strmid(newnames[i], jj, kk-jj) $
        + '.fits'
   endfor

   outfile = [ djs_filepath(logfile, root_dir=extractdir), $
               djs_filepath(plotfile, root_dir=extractdir), $
               newnames ]

   yanny_free, pp

   return, infiles
end

;------------------------------------------------------------------------------
function batch2d_combfiles, planfile, outfile=outfile

   nplan = n_elements(planfilecomb)
   if (nplan GT 1) then begin
      rawfiles = batch2d_combfiles(planfile[0], outfile=outfile)
      for i=1, nplan-1 do begin
         infiles = [infiles, $
          batch2d_rawfiles(planfile[i], outfile=tmpout)]
         outfile = [outfile, tmpout]
      endfor
      return, infiles
   endif

   yanny_read, planfile, pp, hdr=hdr
   extractdir = yanny_par(hdr, 'extractdir')
   combinedir = yanny_par(hdr, 'combinedir')
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')
   plotsnfile = yanny_par(hdr, 'plotsnfile')
   fcalibprefix = yanny_par(hdr, 'fcalibprefix')
   combinefile = yanny_par(hdr, 'combinefile')
   for ii=0, n_elements(pp)-1 do begin
      sname = tag_names(*pp[ii], /structure_name)
      if (sname EQ 'SPEXP') then begin
         expfiles = ((*pp[ii]).name)[*]
         expnums = strmid((*pp[ii]).name[0], 11, 8)
      endif
   endfor

   fcorrprefix = 'spFluxcorr-' + expnums

   infiles = [ djs_filepath(expfiles, root_dir=mjdstr) ]

   junk = fileandpath(planfile, path=thisdir)
   outfile = [ djs_filepath(logfile, root_dir=thisdir), $
               djs_filepath(plotfile, root_dir=thisdir), $
               djs_filepath(plotsnfile, root_dir=thisdir), $
               djs_filepath(fcalibprefix+'-*.fits', root_dir=thisdir), $
               djs_filepath(fcorrprefix+'-*.fits', root_dir=thisdir), $
               djs_filepath(combinefile, root_dir=thisdir) ]

   yanny_free, pp

   return, infiles
end

;------------------------------------------------------------------------------
; This doesn't yet check to see if things have already been run for some plates
pro batch2d, platenums, topdir=topdir, mjd=mjd, mjstart=mjstart, mjend=mjend

   if (NOT keyword_set(platenums)) then platenums = '*'
   if (NOT keyword_set(topdir)) then begin
      cd, current=topdir
      cd, topdir
   endif

   ;----------
   ; Create symoblic link from current directory to raw data directory

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable RAWDATA_DIR'
   junk = findfile('rawdata', count=ct)
   if (ct EQ 0) then $
    spawn, 'ln -s ' + rawdata_dir + ' rawdata'

   ;----------
   ; Create list of plate directories

   spawn, 'ls -d ' + string(platenums+' ', $
    format='(99(a," "))'), platedirs
   if (NOT keyword_set(platedirs[0])) then return
   ndir = n_elements(platedirs)

   ;----------
   ; In each plate directory, find all 'spPlancomb*.par' files

   for idir=0, ndir-1 do begin
      planfile = findfile( $
       djs_filepath('spPlancomb*.par', root_dir=platedirs[idir]), count=nfile)

      for ifile=0, nfile-1 do begin
         yanny_read, planfile[ifile], hdr=hdr
         thismjd = yanny_par(hdr, 'MJD')
; GET ONLY THE MJD's THAT WE WANT ???
         if (keyword_set(platelist)) then begin
            platelist = [platelist, platedirs[idir]]
            planlist = [planlist, planfile[ifile]]
         endif else begin
            platelist = platedirs[idir]
            planlist = planfile[ifile]
         endelse
      endfor
   endfor

   nplate = n_elements(planlist)

   ;----------
   ; Create pointers for input and output files

   pinfile = ptrarr(nplate)
   poutfile = ptrarr(nplate)

   ;----------
   ; For each combine plan file, generate the IDL script files

   fq = "'"
   fullscriptfile = strarr(nplate)
   for iplate=0, nplate-1 do begin
      ; Find all relevant 2D plan files
      yanny_read, planlist[iplate], hdr=hdr
      planfile2d = yanny_par(hdr, 'planfile2d')

      ; Split the combine plan file name into a directory and file name
      planfilecomb = fileandpath(planlist[iplate], path=pathcomb)

      ; Construct the name of the batch file
      i = rstrpos(planfilecomb, '.')
      if (i EQ -1) then i = strlen(planfilecomb)
      fullscriptfile[iplate] = djs_filepath(strmid(planfilecomb,0,i)+'.batch', $
       root_dir=pathcomb)

      ; Write the batch file
      openw, olun, fullscriptfile[iplate], /get_lun
      printf, olun, '; Auto-generated batch file '+systime()
      printf, olun, 'cd, ' + fq+pathcomb+fq
      for i=0, n_elements(planfile2d)-1 do $
       printf, olun, 'spreduce2d, ' + fq+planfile2d[i]+fq
      printf, olun, 'spcombine, ' + fq+planfilecomb+fq
      printf, olun, 'exit'
      close, olun
      free_lun, olun

      ; List of input files
      planfile2d = filepath(planfile2d, root_dir=pathcomb)
      rawfiles = batch2d_rawfiles(planfile2d, outfile=outfile2d)
      rawfiles = djs_filepath(rawfiles, root_dir='rawdata')
      outfile2d = djs_filepath(outfile2d, root_dir=pathcomb)
      junk = batch2d_combfiles(planlist[iplate], outfile=outfilecomb)

      pinfile[iplate] = ptr_new([ fullscriptfile[iplate], $
       planlist[iplate], planfile2d, rawfiles ])

      ; List of output files
      poutfile[iplate] = ptr_new([ outfile2d, outfilecomb ])

   endfor

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

   command = 'nice +19 idl ' + fullscriptfile
   batch, topdir, pinfile, poutfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   ;----------
   ; Remove symbolic link to raw data

   spawn, 'rm -f rawdata'

   return
end
;------------------------------------------------------------------------------
