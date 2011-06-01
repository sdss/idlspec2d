;------------------------------------------------------------------------------
function bbspec_batch, idlcmd, pbsfile=pbsfile1, _EXTRA=Extra

   fq = "'"
   cmd = 'idl -e "'+idlcmd
   if (size(Extra,/tname) EQ 'STRUCT') then begin
      tags = tag_names(Extra)
      for i=0, n_elements(tags)-1 do begin
         qstring = size(Extra.(i),/tname) EQ 'STRING'
         thisquote = qstring ? fq : ''
         num = n_elements(Extra.(i))
         cmd += ','+tags[i]+'='+(num GT 1 ? '[':'') $
          +thisquote+strtrim(Extra.(i)[0],2)+thisquote
         for j=1, num-1 do $
          cmd += ','+strtrim(Extra.(i)[j],2)
         cmd += (num GT 1 ? ']':'')
      endfor
   endif
   cmd += '"'
   splog, cmd

   if (keyword_set(pbsfile1)) then begin
      pbsfile = pbsfile1
   endif else begin
      outfile = Extra.(where(tags EQ 'OUTFILE'))
      pbsfile = repstr(outfile,'.fits','')
   endelse
   openw, olun, pbsfile, /get_lun
   printf, olun, '# Auto-generated batch file '+systime()
   printf, olun, '#PBS -l nodes=1'
   printf, olun, '#PBS -l walltime=48:00:00'
   printf, olun, '#PBS -W umask=0022'
   printf, olun, '#PBS -V'
   printf, olun, '#PBS -j oe'
   printf, olun, 'cd $PBS_O_WORKDIR'
   printf, olun, 'set -o verbose'
;   printf, olun, 'setup idlspec2d '+(strsplit(idlspec2d_version(),/extract))[0]
   printf, olun, cmd
   close, olun
   free_lun, olun

   splog, 'Submitting file '+pbsfile
   spawn, 'qsub '+pbsfile, jobid

   return, jobid
end
;------------------------------------------------------------------------------
