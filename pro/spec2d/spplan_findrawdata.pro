;------------------------------------------------------------------------------
function spplan_findrawdata, inputdir, nfile

   fullnames = findfile(filepath('sdR*.fit', root_dir=inputdir), count=nfile)
   gzipnames = findfile(filepath('sdR*.fit.gz', root_dir=inputdir), count=n)

   if n EQ 0 then return, fullnames

   place = rstrpos(gzipnames,'.gz')  
   for i=0,n-1 do begin

     if place[i] GT 1 then begin
        tempname = strmid(gzipnames[i], 0, place[i])

        if fullnames[0] EQ '' then fullnames = tempname $
        else fullnames = [fullnames, tempname]
     endif
   endfor

   ss = sort(fullnames)
   fullnames = fullnames[ss]
   nfile = n_elements(fullnames)

   return, fullnames
end
