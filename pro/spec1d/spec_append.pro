;------------------------------------------------------------------------------
; Append the array ARG2 to the array ARG1.
; If the first dimension of these arrays is the same, then append
; as [[ARG1],[ARG2]].  If the first dimension of one array is larger
; than the other, then increase the first dimension of the smaller array
; (and zero-fill).

pro spec_append, arg1, arg2

   if (n_elements(arg1) EQ 0) then begin
      arg1 = arg2
      return
   endif
   if (n_elements(arg2) EQ 0) then return

   dims1 = size([arg1], /dimens)
   dims2 = size([arg2], /dimens)
   itype = size(arg1, /type)

   if (dims2[0] EQ dims1[0]) then begin
      ; Case where the new data has the same number of points
      arg1 = [[arg1], [arg2]]
   endif else if (dims2[0] LT dims1[0]) then begin
      ; Case where the new data has fewer points
      dims3 = dims2
      dims3[0] = dims1[0]
      arg3 = make_array(dimension=dims3, type=itype)
      arg3[0:dims2[0]-1,*] = arg2
      arg1 = [[arg1], [arg3]]
   endif else begin
      ; Case where the new data has more points
      dims3 = dims1
      dims3[0] = dims2[0]
      arg3 = make_array(dimension=dims3, type=itype)
      arg3[0:dims1[0]-1,*] = arg1
      arg1 = [[arg3], [arg2]]
   endelse

   return
end

;------------------------------------------------------------------------------
