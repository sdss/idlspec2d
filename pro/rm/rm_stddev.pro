function rm_stddev,array,dim=dim
   if (NOT keyword_set(dim)) then begin 
       st = stddev(array)
   endif else begin
       si=size(array)
       ndim=si[0]
       if dim lt ndim then begin
          st=stddev(array)
       endif else begin
          if ndim eq 1 then begin
             st=stddev(array)
          endif
          if ndim eq 2 then begin
             if dim eq 1 then begin
                st=fltarr(si[2])
                for i=0L, si[2]-1 do begin
                   st[i]=stddev(array[*,i])
                endfor
             endif
             if dim eq 2 then begin
                st=fltarr(si[1])
                for i=0L, si[1]-1 do begin
                   st[i]=stddev(array[i,*])
                endfor
             endif
          endif
          if ndim eq 3 then begin
             if dim eq 1 then begin
                st=fltarr(si[2],si[3])
                for i=0L, si[2]-1 do begin
                   for j=0L, si[3]-1 do begin
                      st[i,j]=stddev(array[*,i,j])
                   endfor
                endfor
             endif
             if dim eq 2 then begin
                st=fltarr(si[1],si[3])
                for i=0L, si[1]-1 do begin
                   for j=0L, si[3]-1 do begin
                       st[i,j]=stddev(array[i,*,j])
                   endfor
                endfor
             endif
             if dim eq 3 then begin
                st=fltarr(si[1],si[2])
                for i=0L, si[1]-1 do begin
                   for j=0L, si[2]-1 do begin
                       st[i,j]=stddev(array[i,j,*])
                   endfor
                endfor
             endif
          endif
       endelse
   endelse
   return, st
end
