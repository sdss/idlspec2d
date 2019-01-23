; supply a 2d array, return the quantile for the nth dimension

function quantile_2d, frac, arr, dim = dim

if not keyword_set(dim) then dim = 1L

arr_size = size(arr)
if dim eq 1L and arr_size[0] ne 1 then nn = arr_size[dim+1] else nn = arr_size[dim-1]
result = dblarr(nn)

for i=0L, nn - 1L do begin
  if dim eq 1L then src = arr[*,i] else src = arr[i,*]
  nnn = n_elements(src)
 
;  range = [min(src), max(src)]
;  result[i] = quantile(frac, src, range = range)
  ind = sort(src)
  src = src[ind]
  ind_quantile =  floor(nnn*frac)
 
  result[i] = src[ind_quantile]

endfor

return, result

end
