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
;   clobber=clobber, gzip=gzip, minflat=minflat, maxflat=maxflat, $
;   pbsnodes=pbsnodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, $
;   pbs_walltime=pbs_walltime, pbs_scratch=pbs_scratch, ember=ember]

; ARGUMENTS:
;  indir: raw data directory, defaulting to $BOSS_SPECTRO_DATA
;  outdir: top directory to write processed files, defaulting
;          to $PROCDATA_DIR if set, otherwise working dir ('.')
;  mjdlist: list of MJDs to process, defaulting to all MJDs
;           found under indir
;  clobber: set this to foce clobbering of existing image and invvar files.
;  gzip: set this to cause output files to be gzipped (default is no gzip)
;  minflat: sdssproc "minflat" keyword, default to 0.8
;  maxflat: sdssproc "maxflat" keyword, default to 1.2
;  pbs_nodes  - If set, collect the pbs qsub commands into pbs_nodes script files
;                in order to run on clusters without node sharing (ie Utah).
;                default to node sharing, and keep the pbs qsub commands in the
;                individual plate-mjd SCRIPT files.
;   pbs_ppn    - If set, use #PBS -l nodes=1:ppn=pbs_ppn, otherwise
;                default to #PBS -l nodes=1
;   pbs_a      - If set, use #PBS -A pbs_a, otherwise
;                default to none
;   pbs_walltime   - If set, use #PBS -l walltime=pbs_walltime, otherwise
;                default to none
;   pbs_scratch- If set, then write to the pbs_scratch  before rsync'ing with
;                the plate-mjd output directories in $BOSS_SPECTRO_REDUX,
;                default to none
;   ember      - If set, then setup the defaults for the University of Utah:
;                pbs_nodes = 8 (for 8 nodes, without node sharing) 
;                pbs_ppn = 12 (12 processors per node)
;                pbs_a = 'bolton-em' (Bolton's account, limited to 8 nodes: ember253 - ember260)
;                pbs_walltime='48:00:00'
;                pbs_scratch = $PBS_SCRATCH/userID/procdata/'  if $PBS_SCRATCH is set, otherwise /scratch/local/userID/procdata/
;
;
; WRITTEN:
;  bolton@utah 2011jan
;   10-Jan-2011  modified  by Joel R. Brownstein, University of Utah, 
;                to generalize to cluster computers by writing PBS commands to bundled pbs batch files,
;                which source job files containing the calls to sdssproc,
;                via the keywords pbs_nodes, pbs_ppn, pbs_a
;                and with University of Utah defaults preset via the
;                keyword ember.
;   14-jan-2011 bolton fix of non-pbs bug.
;   14-Jan-2011 Brownstein added pbs_scratch
;-

pro run_sdssproc, indir=indir, outdir=outdir, mjdlist=mjdlist, clobber=clobber, $
 gzip=gzip, minflat=minflat, maxflat=maxflat, $
 pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, $
 pbs_walltime=pbs_walltime, pbs_scratch=pbs_scratch, ember=ember

  ; Defaults and prep:
  if (n_elements(minflat) ne 1) then minflat = 0.8
  if (n_elements(maxflat) ne 1) then maxflat = 1.2
  if (not keyword_set(indir)) then indir = getenv('BOSS_SPECTRO_DATA')
  if (not keyword_set(outdir)) then if getenv('PROCDATA_DIR') ne '' then outdir=getenv('PROCDATA_DIR') else outdir = '.'
  if (not keyword_set(mjdlist)) then mjdlist = fileandpath(file_search(indir+'/?????'))
  mjdlist = fileandpath(mjdlist)
  
  if keyword_set(pbs_nodes) or keyword_set(ember) then begin ; initialize for a pbs/ember user 
    home = getenv('HOME')
    if (home ne '') then begin
      pos0 = strpos(home,'/',/reverse_search)+1
      pos1 = strlen(home)
      userID = strmid(home,pos0,pos1-pos0)
    endif else userID = 'user'
    print, 'RUN_SDSSPROC: Starting for user:  ',userID
  endif

  if keyword_set(ember) then begin
     if not keyword_set(pbs_nodes) then pbs_nodes=8
     if not keyword_set(pbs_ppn) then pbs_ppn=12
     if not keyword_set(pbs_walltime) then pbs_walltime='48:00:00'
     if not keyword_set(pbs_a) then pbs_a = 'bolton-em'
     if not keyword_set(pbs_scratch) then begin
       scratch_dir = getenv('PBS_SCRATCH')
       if (scratch_dir eq '') then scratch_dir = '/scratch/local/'
       pbs_scratch = djs_filepath(userID+'/procdata',root_dir=scratch_dir)
     endif
   endif
  
  s_mjdlist = strtrim(string(mjdlist), 2)
  nmjd = n_elements(s_mjdlist)
  
  if (nmjd lt 1) then begin
     print, 'No MJDs specified or found.'
     return
  endif
  
  if keyword_set(pbs_nodes) then begin ; initialize the pbs node/processor hierarchy      
    pbs_root_dir = getenv('BOSS_PBS')
    if (pbs_root_dir eq '') then pbs_root_dir = djs_filepath('pbs',root_dir=outdir) $
    else pbs_root_dir = djs_filepath('procdata',root_dir=pbs_root_dir)
    pbs_dir = djs_filepath(userID,root_dir=pbs_root_dir) + '/'
    if file_test(pbs_dir) then begin
       shift_pbs_dir = djs_filepath(userID+'.*',root_dir=pbs_root_dir) + '/'
       shift_pbs = file_search(shift_pbs_dir, count=nshift_pbs)
       max_shift = -1L
       for i=0,nshift_pbs-1 do begin
         pos0 = strpos(shift_pbs[i],'/'+userID+'.',/reverse_search)+strlen('/'+userID+'.')
         pos1 = strlen(shift_pbs[i])
         next_shift = fix(strmid(shift_pbs[i],pos0,pos1-pos0))
         max_shift = (next_shift gt max_shift) ? next_shift : max_shift
       endfor
       shift_pbs_dir = djs_filepath(userID+'.'+strtrim(max_shift+1,2),root_dir=pbs_root_dir) + '/'
       print, 'RUN_SDSSPROC: Renaming previous PBS directory to: '+shift_pbs_dir
       file_move, pbs_dir, shift_pbs_dir
       file_mkdir, pbs_dir
    endif else file_mkdir, pbs_dir
  
    source_cmd = 'source '
    if keyword_set(gzip) then tar_cmd = 'tar -czf ' else tar_cmd = 'tar -cf '
    linefeed = String(10B)
    
    pbs_batch_dir = djs_filepath('batch',root_dir=pbs_dir)
    pbs_job_dir = djs_filepath('job/' + s_mjdlist,root_dir=pbs_dir)  
    if not file_test(pbs_batch_dir) then file_mkdir, pbs_batch_dir
    pbs_node_index = 'node'+ strtrim(indgen(pbs_nodes),2)
    pbs_node_file = djs_filepath(pbs_node_index+'.pbs',root_dir=pbs_batch_dir)
    if keyword_set(pbs_ppn) then begin
      pbs_ppn_index  = '_proc'+ strtrim(indgen(pbs_ppn),2)
      pbs_batch_file = strarr(pbs_nodes,pbs_ppn)
      pbs_batch_script = strarr(pbs_nodes,pbs_ppn)
      for pbs_node = 0,pbs_nodes-1 do pbs_batch_file[pbs_node,*] = djs_filepath(pbs_node_index[pbs_node]+pbs_ppn_index+'.pbs',root_dir=pbs_batch_dir)
      pbs_proc = 0
    endif    
    pbs_node = 0
  endif

  ; Loop over MJDs:
  for i = 0L, nmjd - 1 do begin
     out_this = outdir + '/' + s_mjdlist[i]

     if keyword_set(pbs_nodes) then begin
       if (file_test(pbs_job_dir[i]) eq 0) then file_mkdir, pbs_job_dir[i]
       if keyword_set(pbs_scratch) then begin
         scratch_this = djs_filepath(s_mjdlist[i],root_dir=pbs_scratch)
         if (file_test(scratch_this) eq 0) then file_mkdir, scratch_this
       endif
     endif

     if (file_test(out_this) eq 0) then file_mkdir, out_this
       
     sdr_full = file_search(indir + '/' + s_mjdlist[i] + '/' + 'sdR-[b,r][1,2]-????????.fit*')
     nsdr = n_elements(sdr_full)
 
     print, 'MJD = ' + s_mjdlist[i] + ': Found ' + strtrim(string(nsdr),2) + ' sdR files.'
  
     if keyword_set(pbs_nodes) then pbs_job_file = djs_filepath('sdssproc'+strtrim(indgen(nsdr),2)+'.job',root_dir=pbs_job_dir[i])
                 
     for j = 0L, nsdr - 1 do begin
        sdr_file = fileandpath(sdr_full[j])
        sdr_base = (strsplit(sdr_file, '.', /extract))[0]
        outfile = out_this + '/' + sdr_base + '-IMAGE.fits'
        varfile = out_this + '/' + sdr_base + '-INVVAR.fits'
        if keyword_set(pbs_scratch) then begin
           scratch_outfile = djs_filepath(s_mjdlist[i],root_dir=pbs_scratch)
           if (file_test(scratch_this) eq 0) then file_mkdir, scratch_this
           scratch_outfile = scratch_this + '/' + sdr_base + '-IMAGE.fits'
           scratch_varfile = scratch_this + '/' + sdr_base + '-INVVAR.fits'
        endif
        
        if (keyword_set(clobber) or (file_test(outfile+'*') eq 0) or (file_test(varfile+'*') eq 0)) then begin
           spawn, 'rm -f ' + outfile+'*'
           spawn, 'rm -f ' + varfile+'*'
           
           if keyword_set(pbs_nodes) then begin ;batch mode
             openw, pbs_job, pbs_job_file[j] ,/get_lun
             if keyword_set(pbs_scratch) then begin
               s_outfile = "'" + scratch_outfile + "'"
               s_varfile = "'" + scratch_varfile + "'"
             endif else begin
               s_outfile = "'" + outfile + "'"
               s_varfile = "'" + varfile + "'"
             endelse
             s_sdr_full = "'" + sdr_full[j] + "'"
             printf, pbs_job, 'idl -e "sdssproc, '+s_sdr_full+', outfile='+s_outfile + $
             ', varfile='+s_varfile+', /applycrosstalk, /applypixflat, /applybias, minflat='+strtrim(minflat,2) + $
             ', maxflat='+strtrim(maxflat,2)+', /silent"'
             if keyword_set(gzip) then begin
               if keyword_set(pbs_scratch) then begin
                 printf, pbs_job, 'gzip ' + scratch_outfile
                 printf, pbs_job, 'gzip ' + scratch_varfile
               endif else begin
                 printf, pbs_job, 'gzip ' + outfile
                 printf, pbs_job, 'gzip ' + varfile
               endelse
             endif
             
             if keyword_set(pbs_scratch) then begin ; finish up: tar to a pipe to copy scratch into place
               printf, pbs_job, 'cd '+scratch_this
               printf, pbs_job, 'tar cf - '+sdr_base + '*.* | (cd '+ out_this + ' ; tar xfp - )'
               printf, pbs_job, 'rm -f ' + sdr_base + '*.*
             endif
             
             close, pbs_job
             free_lun, pbs_job
 
           endif else begin                          ; interactive mode
             print, ' File ' + strtrim(string(j),2) + ' (' + sdr_base + ')'
             sdssproc, sdr_full[j], outfile=outfile, varfile=varfile, /applycrosstalk, $
              /applypixflat, /applybias, minflat=minflat, maxflat=maxflat, /silent
              if keyword_set(gzip) then begin
                spawn, 'gzip ' + outfile
                spawn, 'gzip ' + varfile
             endif
           endelse
           
        endif

        if keyword_set(pbs_nodes) then begin ; finalize the pbs_job's by delegating to pbs_node's and pbs_batch_script's
          if keyword_set(pbs_ppn) then begin
            pbs_batch_script[pbs_node,pbs_proc] += source_cmd + pbs_job_file[j] + linefeed
            pbs_proc += 1
            if pbs_proc ge pbs_ppn then begin
              pbs_node += 1
              pbs_proc = 0
              if pbs_node ge pbs_nodes then pbs_node = 0
            endif
          endif
        endif
    
     endfor
          
  endfor
  
  pbs_batch_complete = 0L
  if keyword_set(pbs_nodes) then begin ; write the pbs_batch_script's to pbs_batch_file's
    for pbs_node = 0, pbs_nodes-1 do begin
      openw, pbs_node_batch, pbs_node_file[pbs_node] ,/get_lun
      printf, pbs_node_batch, '# Auto-generated batch file '+systime()
      if keyword_set(pbs_a) then printf, pbs_node_batch, '#PBS -A '+pbs_a
      if keyword_set(pbs_walltime) then printf, pbs_node_batch, '#PBS -l walltime='+pbs_walltime
      printf, pbs_node_batch, '#PBS -W umask=0022'
      printf, pbs_node_batch, '#PBS -V'
      printf, pbs_node_batch, '#PBS -j oe'
      if keyword_set(pbs_ppn) then begin
        printf, pbs_node_batch, '#PBS -l nodes=1:ppn='+strtrim(pbs_ppn,2)
        for pbs_proc = 0, pbs_ppn-1 do begin
          pbs_batch_complete = (pbs_batch_script[pbs_node,pbs_proc] eq '')
          if not pbs_batch_complete then begin
            openw, pbs_batch, pbs_batch_file[pbs_node,pbs_proc], append=pbs_ppn_append, /get_lun
            printf, pbs_batch, pbs_batch_script[pbs_node,pbs_proc]
            close, pbs_batch
            free_lun, pbs_batch
            printf, pbs_node_batch, source_cmd+pbs_batch_file[pbs_node,pbs_proc] + ' &'
          endif else break
        endfor
      endif else printf, pbs_node_batch, '#PBS -l nodes=1"
      printf, pbs_node_batch, 'wait'
      close, pbs_node_batch
      free_lun, pbs_node_batch
      if pbs_batch_complete then break
    endfor
    
    cd, pbs_batch_dir
    for pbs_node=0L, pbs_nodes-1 do spawn, 'qsub '+pbs_node_file[pbs_node]
    
  endif
  
  return
end
