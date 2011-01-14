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
;   pbsnodes=pbsnodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_walltime=pbs_walltime, ember=ember]

; ARGUMENTS:
;  indir: raw data directory, defaulting to $RAWDATA_DIR
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
;   pbs_dir    - Write all pbs batch files to $pbs_dir/batch and 
;                all pbs script jobs to $pbs_dir/job
;                default to $PROCDATA_PBS
;   pbs_a      - If set, use #PBS -A pbs_a, otherwise
;                default to none
;   pbs_walltime   - If set, use #PBS -l walltime=pbs_walltime, otherwise
;                default to none
;   ember      - If set, then setup the defaults for the University of Utah:
;                pbs_nodes = 8 (for 8 nodes, without node sharing) 
;                pbs_ppn = 12 (12 processors per node)
;                pbs_a = 'bolton-em' (Bolton's account, limited to 8 nodes: ember253 - ember260)

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
;-

pro run_sdssproc, indir=indir, outdir=outdir, mjdlist=mjdlist, clobber=clobber, $
 gzip=gzip, minflat=minflat, maxflat=maxflat, $
 pbs_dir=pbs_dir, pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_walltime=pbs_walltime, ember=ember

  ; Defaults and prep:
  if (n_elements(minflat) ne 1) then minflat = 0.8
  if (n_elements(maxflat) ne 1) then maxflat = 1.2
  if (not keyword_set(indir)) then indir = getenv('RAWDATA_DIR')
  if (not keyword_set(outdir)) then if getenv('PROCDATA_DIR') ne '' then outdir=getenv('PROCDATA_DIR') else outdir = '.'
  if (not keyword_set(mjdlist)) then mjdlist = fileandpath(file_search(indir+'/?????'))
     mjdlist = fileandpath(mjdlist)
   if keyword_set(ember) then begin
     pbs_nodes=7
     pbs_ppn=12
     pbs_walltime='120:00:00'
     pbs_a = 'bolton-em'
   endif
  
  s_mjdlist = strtrim(string(mjdlist), 2)
  nmjd = n_elements(s_mjdlist)
  
  if (nmjd lt 1) then begin
     print, 'No MJDs specified or found.'
     return
  endif
  
  if keyword_set(pbs_nodes) then begin ; initialize the pbs node/processor hierarchy 
    source_cmd = 'source '
    linefeed = String(10B)
    if (not keyword_set(pbs_dir)) then pbs_dir=getenv('PROCDATA_PBS') else pbs_dir = 0L
    pbs_batch_dir = djs_filepath('batch',root_dir=pbs_dir)
    pbs_job_dir = djs_filepath('job/' + s_mjdlist,root_dir=pbs_dir)  
    pbs_job_dir = djs_filepath('job/' + s_mjdlist,root_dir=pbs_dir)  
    if keyword_set(pbs_dir) then begin
      if file_test(pbs_dir) and (file_test(pbs_batch_dir) or file_test(djs_filepath('job/',root_dir=pbs_dir))) then begin
         pos = strpos(pbs_dir,'/',strlen(pbs_dir)-1)
         if pos ge 0 then pbs_dir=strmid(pbs_dir,0,pos)
         shift_pbs_dir = pbs_dir+'.*'
         shift_pbs = file_search(shift_pbs_dir, count=nshift_pbs)
         max_shift = -1L
         for i=0,nshift_pbs-1 do begin
           pos0 = strpos(shift_pbs[i],'.',/reverse_search)+1
           pos1 = strlen(shift_pbs[i])
           next_shift = fix(strmid(shift_pbs[i],pos0,pos1-pos0))
           max_shift = (next_shift gt max_shift) ? next_shift : max_shift
        endfor
        shift_pbs_dir = pbs_dir + '.' + strtrim(max_shift+1,2)
        print, 'RUN_SDSSPROC: Renaming previous PBS directory to: '+shift_pbs_dir
        file_move, pbs_dir, shift_pbs_dir
        file_mkdir, pbs_dir
      endif
    endif else begin
      print, 'Please use the pbs_dir keyword, or set the $PROCDATA_PBS environment variable'
      return
    endelse
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
     if (file_test(out_this) eq 0) then file_mkdir, out_this
     if keyword_set(pbs_nodes) then if (file_test(pbs_job_dir[i]) eq 0) then file_mkdir, pbs_job_dir[i]
       
     sdr_full = file_search(indir + '/' + s_mjdlist[i] + '/' + 'sdR-[b,r][1,2]-????????.fit*')
     nsdr = n_elements(sdr_full)
 
     print, 'MJD = ' + s_mjdlist[i] + ': Found ' + strtrim(string(nsdr),2) + ' sdR files.'
  
     if keyword_set(pbs_nodes) then pbs_job_file = djs_filepath('script-'+strtrim(indgen(nsdr),2),root_dir=pbs_job_dir[i])
     ;if keyword_set(pbs_nodes) then pbs_job_file = djs_filepath('sdssproc'+strtrim(indgen(nsdr),2)+'.job',root_dir=pbs_job_dir[i])
                 
     for j = 0L, nsdr - 1 do begin
        sdr_file = fileandpath(sdr_full[j])
        sdr_base = (strsplit(sdr_file, '.', /extract))[0]
        ;print, ' File ' + strtrim(string(j),2) + ' (' + sdr_base + ')'
        outfile = out_this + '/' + sdr_base + '-IMAGE.fits'
        varfile = out_this + '/' + sdr_base + '-INVVAR.fits'
        
        if (keyword_set(clobber) or (file_test(outfile+'*') eq 0) or (file_test(varfile+'*') eq 0)) then begin
           spawn, 'rm -f ' + outfile+'*'
           spawn, 'rm -f ' + varfile+'*'
           
           if keyword_set(pbs_nodes) then begin ;batch mode
             openw, pbs_job, pbs_job_file[j] ,/get_lun
             s_outfile = "'" + outfile + "'"
             s_varfile = "'" + varfile + "'"
             s_sdr_full = "'" + sdr_full[j] + "'"
             printf, pbs_job, 'idl -e "sdssproc, '+s_sdr_full+', outfile='+s_outfile + $
             ', varfile='+s_varfile+', /applycrosstalk, /applypixflat, /applybias, minflat='+strtrim(minflat,2) + $
             ', maxflat='+strtrim(maxflat,2)+', /silent"'
             if keyword_set(gzip) then begin
                printf, pbs_job, 'gzip ' + outfile
                printf, pbs_job, 'gzip ' + varfile
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
