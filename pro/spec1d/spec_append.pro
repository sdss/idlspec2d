;------------------------------------------------------------------------------
; Append the array ARG2 to the array ARG1.
; If the first dimension of these arrays is the same, then append
; as [[ARG1],[ARG2]].  If the first dimension of one array is larger
; than the other, then increase the first dimension of the smaller array
; (and zero-fill).

pro spec_append, arg1, arg2, pixshift

   if (n_elements(arg1) EQ 0) then begin
      arg1 = arg2
      return
   endif
   if (n_elements(arg2) EQ 0) then return

   if (size(arg1,/n_dimen) GT 2 OR size(arg2,/n_dimen) GT 2) then $
    message, 'Dimensions of input arrays must be <= 2'

   dims1 = size(arg1, /dimens)
   dims2 = size(arg2, /dimens)
   dims1 = dims1 > 1 ; Fix for new behaviour in IDL 5.4 and later.
   dims2 = dims2 > 1 ; Fix for new behaviour in IDL 5.4 and later.

   ; Promote the type of variable to the bigger (i.e., if one float
   ; and one double, output type double).
   itype = size(arg1[0]+arg2[0], /type)

   ;----------
   ; Decide how much to pre-padd each array based upon PIXSHIFT

   nadd1 = 0
   nadd2 = 0
   if (keyword_set(pixshift)) then begin
      if (pixshift LT 0) then begin
         ; Case where we first have to pre-pad the ARG1 spectra.
         nadd1 = -pixshift
         nadd2 = 0
      endif else if (pixshift GT 0) then begin
         ; Case where we first have to pre-pad the ARG2 spectra.
         nadd1 = 0
         nadd2 = pixshift
      endif
   endif

   ;----------
   ; Construct the output array

   nrow1 = n_elements(arg1) / dims1[0]
   nrow2 = n_elements(arg2) / dims2[0]

   dims3 = [max([dims1[0]+nadd1, dims2[0]+nadd2]), nrow1+nrow2]
   arg3 = make_array(dimension=dims3, type=itype)

   ;----------
   ; Fill the output array with the input arrays

   arg3[nadd1:nadd1+dims1[0]-1,0:nrow1-1] = arg1
   arg3[nadd2:nadd2+dims2[0]-1,nrow1:nrow1+nrow2-1] = arg2

   arg1 = temporary(arg3)

   return
end
;------------------------------------------------------------------------------
