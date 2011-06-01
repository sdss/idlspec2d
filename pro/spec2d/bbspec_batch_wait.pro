;------------------------------------------------------------------------------
; Wait for all PBS jobs to complete, then return
pro bbspec_batch_wait, jobid

   joblist = ''
   for i=0, n_elements(jobid)-1 do joblist += ' '+strtrim(jobid[i])
   retval = 1B
   while (keyword_set(retval)) do begin
      spawn, 'qstat'+joblist, retval, reterr
      wait, 5
   end

   return
end
;------------------------------------------------------------------------------
