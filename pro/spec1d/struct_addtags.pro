;+
; NAME:
;   struct_addtags
;
; PURPOSE:
;   Add tags from one structure (array) to another
;
; CALLING SEQUENCE:
;   outstruct = struct_addtags(astruct, bstruct)
;
; INPUTS:
;   astruct    - First structure, which can be an array
;   bstruct    - Second structure, which can be an array
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   outstruct  - Ouput structure array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function struct_addtags, astruct, bstruct

   if (N_elements(astruct) EQ 0) then $
    return, bstruct

   num1 = N_elements(astruct)
   num2 = N_elements(bstruct)
   if (num1 NE num2) then begin
      message, 'Both structures must have the same number of elements'
   endif

   ; Create an empty structure with all the tags from both structures
   obj1 = astruct[0]
   tname = tag_names(bstruct)
   ntag = N_tags(bstruct)
   for i=0, ntag-1 do $
    obj1 = create_struct(obj1, tname[i], bstruct[0].(i))
   outstruct = replicate(obj1, num1)

   ; Assign elements from NEWDAT into the new output structure
   newtags = tag_names(outstruct)
   for i=0, ntag-1 do begin
      j = (where(newtags EQ tname[i]))[0]
      outstruct.(j) = bstruct.(i)
   endfor

   return, outstruct
end
;------------------------------------------------------------------------------
