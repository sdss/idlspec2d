;+
; NAME:
;   sxcombinepar_v2
;
; PURPOSE:
;   Combine values of specified header cards from many FITS headers.
;
; CALLING SEQUENCE:
;   sxcombinepar, hdrarr, cardname, outhdr, [ func=, weights=, /zeros, $
;    outcard=, _EXTRA=KeywordsForSxaddpar ]
;
; INPUTS:
;   hdrarr     - Array of pointers to FITS headers
;   cardname   - Name(s) of header cards to average
;   outhdr     - Output header
;
; OPTIONAL KEYWORDS:
;   func       - Function to apply:
;                  'average'
;                  'median'
;                  'min'
;                  'max'
;                  'total'
;                Default to 'average'
;   weights    - If set, then weight each of the input headers by these weights;
;                only applicable when the function type is 'average'.
;   zeros      - If set, then include zero values when determining the
;                average or other function.  But never use the zeros
;                returned by SXPAR() if a header is missing that card
;                altogether.
;   outcard    - Card name(s) in output header; if not specified, then use
;                the same name as in CARDNAME.
;   comments   - Comment to replace or add a comment to the header, can use a dictionary
;                with the keys being the cardname, or a single string
;   _EXTRA     - Optional keywords for SXADDPAR (such as BEFORE,AFTER,FORMAT).
;
; OUTPUTS:
;   outhdr     - (Modified.)
;
; COMMENTS:
;
; BUGS:
;
; PROCEDURES CALLED:
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   31-Jan-2002  Written by D. Schlegel, Princeton
;   04-Feb-2022  Modified by S. Morrison to add comments

;-
;------------------------------------------------------------------------------
pro sxcombinepar_v2, hdrarr, cardname, outhdr, func=func, camnames=camnames, $
          weights=weights, zeros=zeros, outcard=outcard, comments=comments, $
          _EXTRA=KeywordsForSxaddpar

   if (n_params() LT 3) then begin
      print, 'Syntax - sxcombinepar, hdrarr, cardname, outhdr, [ func=, /zeros, outcard=, comments=]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'average'
   if (n_elements(zeros) EQ 0) then zeros = 0
   if (NOT keyword_set(outcard)) then outcard = cardname

   if (n_elements(outcard) NE n_elements(cardname)) then $
    message, 'Number of elements in OUTCARD and CARDNAME must agree'
   if keyword_set(weights) then begin
       if (n_elements(hdrarr) NE n_elements(weights)) then $
         message, 'Number of elements in HDRARR and WEIGHTS must agree'
   endif

   ;----------
   ; Call this routine recursively if CARDNAME has multiple elements

   ncard = n_elements(cardname)
   if (ncard EQ 0) then return
   if (ncard GT 1) then begin
      for icard=0, ncard-1 do begin
         sxcombinepar_v2, hdrarr, cardname[icard], outhdr, func=func, camnames=camnames, $
          weights=weights, zeros=zeros, outcard=outcard[icard], comments=comments
      endfor
      return
   endif

   allweights = []

   for ihdr=0, n_elements(hdrarr)-1 do begin
      thiscam = sxpar(*hdrarr[ihdr], 'CAMERAS')
      if keyword_set(camnames) then begin
        if not strmatch(camnames[0], '*'+strtrim(thiscam,2)+'*', /fold_case) then continue
      endif
      thisval = sxpar(*hdrarr[ihdr], cardname[0], count=thisct)
      if (thisct GT 0 AND (keyword_set(thisval) or keyword_set(zeros))) then begin
         if (NOT keyword_set(allval)) then allval = thisval $
         else allval = [allval, thisval]
         if keyword_set(weights) then begin
            allweights = [allweights, weights[ihdr]]
         endif else begin
            allweights = [allweights,1]
         endelse
      endif
   endfor

   nval = n_elements(allval)
   if (nval GT 0) then begin
      if not strmatch(typename(allval[0]),'string',/fold_case) then begin
        case strlowcase(func) of
            'average': begin
                if total(allweights) eq 0.0 then allweights[*] = 1.0
                outval = total(allval * allweights) / total(allweights)
                end
            'median' : outval = median(allval)
            'min'    : outval = min(allval)
            'max'    : outval = max(allval)
            'total'  : outval = total(allval)
            else     : message, 'Invalid FUNCTION'
        endcase
        if ~finite(outval) then outval = string(outval)
      endif else outval = allval[0]
      
      if keyword_set(comments) then begin
        if isa(comments,/str) then begin
           sxaddpar, outhdr, outcard[0], outval,comments, _EXTRA=KeywordsForSxaddpar
        endif else begin
           if not comments.haskey(key) then begin
               sxaddpar, outhdr, outcard[0], outval, _EXTRA=KeywordsForSxaddpar
           endif else begin
               sxaddpar, outhdr, outcard[0], outval, comments[outcard[0]], _EXTRA=KeywordsForSxaddpar
           endelse
        endelse
      endif else sxaddpar, outhdr, outcard[0], outval, _EXTRA=KeywordsForSxaddpar
   endif

   return
end
;------------------------------------------------------------------------------
