;+
; NAME:
;   zrefind
;
; PURPOSE:
;   Re-fit redshifts from a previous call to ZFIND(), but doing a local
;   fit around the previous answers.
;
; CALLING SEQUENCE:
;   result = zrefind( objflux, objivar, pwidth=, pspace=, width=, zold= ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;
; REQUIRED KEYWORDS:
;   zold       - Z structure from an initial call to ZFIND().
;
; OPTIONAL KEYWORDS:
;   pwidth     - Search width in pixels about the intial redshift; default to 5
;   pspace     - Keyword for ZCOMPUTE().
;   width      - Keyword for ZCOMPUTE().
;
; OUTPUTS:
;   result     - Structure with redshift-fit information, modified from ZOLD.
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
;   zfind()
;
; REVISION HISTORY:
;   17-Aug-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function zrefind, objflux, objivar, hdr=hdr, $
 pwidth=pwidth, pspace=pspace, width=width, zold=zold, _EXTRA=EXTRA

   if (NOT keyword_set(pwidth)) then pwidth = 5

   ;----------
   ; Copy results

   result = zold
   nres = n_elements(result)

   if (size(zold,/n_dim) EQ 1) then nfind = 1 $
    else nfind = (size(zold,/dimens))[0]

   ;----------
   ; Identify which redshift measurements are valid and will be re-measured

   qvalid = result.tfile NE '' AND result.tcolumn[0] NE -1 $
    AND result.dof NE 0
   ivalid = where(qvalid)
   if (ivalid[0] EQ -1) then return, result

   ;----------
   ; Find all eigen-spectra files being used
   ; The ISELECT array will give indices for unique (and valid) sort strings.

   ntcol = n_elements(result[0].tcolumn)
   tcolstring = strarr(nres)
   for ires=0, nres-1 do $
    tcolstring = string(result[ires].tcolumn, format='('+string(ntcol)+'i)')
   sortstring = result[*].tfile + tcolstring + string(result[*].npoly)

   isort = ivalid[ sort(sortstring[ivalid]) ]
   iselect = isort[ uniq(sortstring[isort]) ]

   for ifile=0, n_elements(iselect)-1 do begin

      i0 = iselect[ifile] ; Index of result specifying this TFILE,TCOLUMN,NPOLY

      ; Identify which redshifts are measured with this template file,
      ; TCOLUMN, and NPOLY.
      indx = where(sortstring EQ sortstring[i0])

      ; Find object numbers corresponding to these indices
      iobj = indx / nfind

      ; Re-fit the redshifts using the specified PWIDTH,PSPACE
      ncol = (reverse(where(result[i0].tcolumn NE -1)))[0] + 1
      res1 = zfind(objflux[*,iobj], objivar[*,iobj], hdr=hdr, $
       eigenfile=result[i0].tfile, columns=result[i0].tcolumn[0:ncol-1], $
       npoly=result[i0].npoly, zguess=result[indx].z, $
       pwidth=pwidth, pspace=pspace, nfind=1, width=width, _EXTRA=EXTRA)

      ; Copy the results into the output structure
      result[indx].z = res1[*].z
      result[indx].z_err = res1[*].z_err
      result[indx].rchi2 = res1[*].rchi2
      result[indx].dof = res1[*].dof
      result[indx].theta = res1[*].theta

   endfor

   return, result
end
;------------------------------------------------------------------------------
