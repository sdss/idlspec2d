
function lookforgzip, filename, count=ct

     f = findfile(filename, count=ct)

     if ct GT 0 then return, f

     f = findfile(filename+'.gz', count=ct)

     if ct GT 0 then return, f

     return, ''
end
