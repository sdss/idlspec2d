;------------------------------------------------------------------------------
function spplan_create_spexp, expnum, plateid, mjd, mapname, flavor, exptime, $
 filename, cameras, minexp=minexp

   if (flavor NE 'flat' AND flavor NE 'arc' $
    AND flavor NE 'science' AND flavor NE 'smear') then $
    return, 0
   if (keyword_set(minexp)) then begin
      if (flavor EQ 'science' AND exptime LT minexp) then $
       return, 0 
   endif

   badname = 'UNKNOWN'
   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   spexp = {spexp, $
    plateid : long(plateid), $
    mjd     : long(mjd), $
    mapname : string(mapname), $
    flavor  : string(flavor), $
    exptime : float(exptime), $
    name    : strarr(ncam) }

   for icam=0, ncam-1 do begin
      ii = where(cameras EQ camnames[icam], ct)
      if (ct GT 1) then $
       message, 'Multiple files with EXPOSURE=' + string(expnum) $
        + ' CAMERAS=' + camnames[icam]
      if (ct EQ 1) then $
       spexp.name[icam] = filename[ii[0]] $
      else $
       spexp.name[icam] = badname
   endfor

   return, spexp
end
