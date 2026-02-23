;+
; NAME:
;   tracefibercrude
;
; PURPOSE:
;   Calling script to return full fiber traces using TRACE_CRUDE.
;
; CALLING SEQUENCE:
;   xset = tracefibercrude( fimage, invvar, [ ystart=, nmed=, $
;    xmask=, yset=, maxerr=, maxshifte=, maxshift0=, xerr=, maxdev=, ngrow=, $
;    fibermask=, cartid=, flathdr=, padding=, plottitle= ] )
;
; INPUTS:
;   fimage     - Image
;   cartid     - Cartridge ID from plugmap
;
; OPTIONAL INPUTS FOR TRACE_CEN:
;   ystart     - Y position in image to search for initial X centers; default
;                to the central row
;   nmed       - Number of rows to median filter around YSTART; default to 21
;   plottitle  -
;
; OPTIONAL INPUTS FOR TRACE_CRUDE:
;   flathdr    - FITS header for determining CARTID and MJD
;   invvar     - Inverse variance (weight) image
;   radius     - Radius for centroiding; default to 3.0
;   maxerr     - Maximum error in centroid allowed for valid recentering;
;                default to 0.2
;   maxshifte  - Maximum shift in centroid allowed for valid recentering;
;                default to 0.1
;   maxshift0  - Maximum shift in centroid allowed for initial row;
;                default to 0.5
;   padding    - number added to ndegree in loop to fix deviant centroids
;
; OPTIONAL INPUTS:
;   maxdev     - Maximum deviation of X in pixels; default to rejecting any
;                XPOS positions that deviate by more than 1.0 pixels from
;                a polynomial remapping of the centroids from other rows.
;   ngrow      - For each trace, replace all centroids within NGROW rows
;                of a bad centroid with the predicted centroid locations.
;                Default to 5.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;   xmask      - Mask set to 1 for good fiber centers, 0 for bad;
;                same dimensions as XSET.
;   fibermask  - (Modified.)
;   bundlefibers
;   nbundle
;
; COMMENTS:
;   Without djs_maskinterp, hot columns skew traces
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   fibermask_bits()
;   trace_crude()
;   trace_fweight()
;
; REVISION HISTORY:
;   13-Sep-1999  Written by David Schlegel, Princeton.
;    8-Jul-2001  Added djs_maskinterp call
;   05-Oct-2010  ASB added masking of rows with invvar all zero
;-
;------------------------------------------------------------------------------
function tracefibercrude, image, invvar, ystart=ystart, nmed=nmed, $
 xmask=xmask, radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshifte, $
 maxshift0=maxshift0, xerr=xerr, maxdev=maxdev, ngrow=ngrow, $
 fibermask=fibermask, cartid=cartid, flathdr=flathdr, padding=padding, $
 plottitle=plottitle, plates=plates, legacy=legacy, bundlefibers=bundlefibers, $
 nbundle=nbundle, flatname=flatname, opfdir = opfdir

   if (NOT keyword_set(maxdev)) then maxdev = 1.0
   if (NOT keyword_set(ngrow)) then ngrow = 5
   if (NOT keyword_set(radius)) then radius = 3.0
   if (NOT keyword_set(padding)) then padding=0

   if keyword_set(plates) then begin
      opfibersFile='opFibersP.par'
   endif else begin
      if keyword_set(legacy) then begin
        opfibersFile='opFibers.par'
      endif else begin
        fps=1
        opfibersFile='opFibersFPS.par'
      endelse
   endelse
 
   if not keyword_set(fps) then begin
        fiberparam = yanny_readone(djs_filepath(opfibersFile, $
            root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles'), 'FIBERPARAM')
        if (NOT keyword_set(fiberparam)) then message, opfibersFile+' file not found!'
 
        camname = strtrim(sxpar(flathdr, 'CAMERAS'),2)
        mjd = sxpar(flathdr, 'MJD')
        if (keyword_set(cartid) * keyword_set(camname) * keyword_set(mjd) EQ 0) $
            then message, 'Must set CARTID, CAMERAS, MJD in flat header!'
        i = where(fiberparam.cartid EQ cartid AND fiberparam.camname EQ camname $
                    AND fiberparam.mjd LE mjd, ct)
        if (ct EQ 0) then $
            message, 'No match for this CARTID + MJD in '+opfibersFile+'!'

        isort = i[reverse(sort(fiberparam[i].mjd))]
        fiberparam = fiberparam[isort[0]]
        nbundle = (long(total(fiberparam.fiberspace NE 0)))[0]
        BUNDLEFIBERS = make_array(25,/integer,value=20)
   endif else begin
        camname = strtrim(sxpar(flathdr, 'CAMERAS'),2)
        mjd = sxpar(flathdr, 'MJD')
        if not keyword_set(opfdir) then opfdir = filepath(getenv('RUN2D'), root_dir = getenv('BOSS_SPECTRO_REDUX'))
        opfile = filepath('opFibersFPS-'+strtrim(mjd,2)+'.par', root_dir=opfdir, subdir=['trace',strtrim(mjd,2)])
        if FILE_TEST(opfile) then begin
            splog, 'Using local opfiber file '+strtrim(opfile,2)
        endif else begin
            splog, 'Using default opfiber file'
            opfile = djs_filepath(opfibersFile, root_dir=getenv('IDLSPEC2D_DIR'),subdir='opfiles')
        end
        
        fiberparam=opfiber_read(opfile, cartid, camname, mjd)
 
        if (keyword_set(cartid) * keyword_set(camname) * keyword_set(mjd) EQ 0) $
                then message, 'Must set CAMERAS, MJD in flat header!'
        i = where(fiberparam.cartid EQ strtrim(cartid,2) AND fiberparam.camname EQ camname $
                    AND fiberparam.mjd LE mjd, ct)
        if (ct EQ 0) then message, 'No match for this MJD in '+opfibersFile+'!'
        nbundle = fiberparam.nbundles
    BUNDLEFIBERS = fiberparam.bundlefibers
   endelse
   ; Assume the fiber bundles used are the first NBUNDLE ones...
   if (total(fiberparam.fiberspace[0:nbundle-1] EQ 0) GT 0) then $
        message, 'Some FIBERSPACE parameters are zero!'

   ;----------
   ; If INVVAR is set, then start by interpolating over bad pixels

   if (keyword_set(invvar)) then $
    fimage = djs_maskinterp(image, (invvar LE 0), iaxis=0) $
   else $
    fimage = image

   ;----------
   ; Find the X-centers in the row specified by YSTART

   ; XGOOD=1 for fibers that were actually found, 0 otherwise
   xposition = trace_cen(fimage, xstart=fiberparam.bundlegap[0], $
                         ystart=ystart, nmed=nmed, $
                         nfiber=fiberparam.nfiber, nbundle=nbundle, $
                         fiberspace=fiberparam.fiberspace[0:nbundle-1], $
                         bundlespace=fiberparam.bundlegap[1:nbundle-1], $
                         bundlefibers=BUNDLEFIBERS, $
                         xgood=xgood, plottitle=plottitle, fluxvec=fluxvec, $
                         fmodel=fmodel, flatname=flatname)

   ntrace = n_elements(xposition)
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(ntrace)
;   forprint, xposition,  /NOCOMMENT, textout=repstr(repstr(flatname, 'sdR','spFlat_v0'),'.fit','.prt')
   ;----------
   ; Trace

   xset = trace_crude(fimage, invvar, xstart=xposition, ystart=ystart, $
                      radius=radius, yset=yset, maxerr=maxerr, $
                      maxshifte=maxshifte, maxshift0=maxshift0, xerr=xerr)
   xmask = xerr LT 990  ; =1 for good centers, =0 for bad

;   WRITE_CSV, repstr(repstr(flatname, 'sdR','spFlat_vp5'),'.fit','.prt'), xset

   ; Mask contributions from completely bad rows (ticket #1025)
   nullrow = total(invvar gt 0., 1) eq 0.
   wh_null = where(nullrow, n_null)
   if (n_null gt 0) then xmask[wh_null,*] = 0B

   ;--------------------------------------------------------------------
   ; Mark this trace as potentially bad (xgood[itrace] = 0)
   ; if either the initial extraction row had bad pixels,
   ; or the initial extraction row was off the left or right
   ; edge of the CCD.

   ncol = (size(invvar,/dimen))[0]
   if ncol GT 0 then begin
     for itrace=0, ntrace-1 do begin
       ix1 = floor(xposition[itrace]-radius) > 0L
       ix2 =  ceil(xposition[itrace]+radius) < (ncol-1)
       if (ix1 GT ncol-1 OR ix2 LT 0) then begin
          ; Case where the initial extraction position
          ; was off the left or right edge of the CCD.
          xgood[itrace] = 0
       endif else begin
          ; Check if any bad pixels at initial centroiding position.
          junk = where(invvar[ix1:ix2,ystart] LE 0, nbad)
          if (nbad GT 0) then xgood[itrace] = 0
       endelse
     endfor
   endif

   ;----------
   ; Compare the traces in each row to those in row YSTART.
   ; Our assumption is that those centers should be a polynomial mapping
   ; of the centers from row YSTART.  Centers that are deviant from this
   ; mapping are replaced with the position predicted by this mapping.

   ny = (size(fimage, /dimens))[1]
   ndegree = 4 ; Five terms

   ;----------
   ; Loop to find all deviant centroids, and add these to the mask XMASK.
   ; XMASK=1 for good.

   for iy=0, ny-1 do begin
      xcheck = xgood AND xmask[iy,*] ; Test for good fiber & good centroid
      if (total(xcheck) GT ((ndegree+2) > (0.2 * ntrace))) then begin
        res=mpfitfun('poly',xposition,reform(xset[iy,*]),xcheck,fltarr(ndegree+1),weights=xcheck,yfit=xfit,/quiet)
        xdiff = xfit - xset[iy,*]
        ibad = where(abs(xdiff) GT maxdev)
        if (ibad[0] NE -1) then xmask[iy,ibad] = 0
      endif else begin
        xmask[iy,*] = 0 ; Too few good centroids in this row; mark all as bad
      endelse
   endfor

   ;----------
   ; Smooth the bad centroids to NGROW adjacent rows (of the same trace)

   for itrace=0, ntrace-1 do begin
      xmask[*,itrace] = smooth( xmask[*,itrace]+0.0, 2*ngrow+1) EQ 1
   endfor

   ;----------
   ; Loop to fix deviant centroids

   for iy=0, ny-1 do begin
      ixbad = where(xmask[iy,*] EQ 0, nbad)
      if (nbad GT 0 AND nbad LT ntrace-1-(ndegree+padding)) then begin
         ixgood = where(xmask[iy,*] EQ 1)

         if (!version.release LT '5.4') then begin
            coeff = polyfitw(xposition, xset[iy,*], xmask[iy,*], $
             ndegree, xfit)
         endif else begin
            coeff = polyfitw(xposition, xset[iy,*], xmask[iy,*], $
             ndegree, xfit, /double)
         endelse

         xset[iy,ixbad] = xfit[ixbad]
      endif
   endfor



   ;----------------------------------------------------------------------
   ;  This is new code to replace traces which have been
   ;  affected by bad columns (or masked pixels in general) at the
   ;  starting row.  Any start position which is offset from the true
   ;  position will produce systematic errors in all of the centroids
   ;  which are corrected above and based on xposition.
   ;
   ;  This algorithm works as follows:
   ;  1) Select from the near eight neighbors (-4 to +4) traces
   ;       which have fewer than 100 bad pixels. (Checktrace)
   ;  2) Calculate the mean trace of the selected neighbors = meantrace
   ;  3) Calculate the offset of the problem trace with the selected neighbors
   ;           a) use only rows (goodrows) which have good centroids in
   ;                        all selected neighbors AND in the problem trace.
   ;           b) calculate median offset w.r.t. neighbors in goodrows only
   ;           c) Use this mean offset to correct the meantrace to the correct
   ;                zero point position

   problemtraces = where(xgood EQ 0, ct)
   nrow = (size(xset,/dimen))[1]
   tracenum = lindgen(ntrace)
   tmp_xpos = trace_fweight(fimage, xset, yset, $
            radius=radius, xerr=tmp_xerr, invvar=invvar)

   xorig = xset
   badpix = total(tmp_xerr EQ 999,1)

   if ct GT 0 then $
     splog, 'Warning: Fixing traces: '+strjoin(strtrim(fix(problemtraces)+1,2),' ')+'(1-indexed)'
   for ii=0, ct-1 do begin
      itrace = problemtraces[ii]

;----------------------------------------------------------------------------
;    Really simple minded loop to check for nearest 8 neighbors who might
;        be suitable for substitution
;
      checktrace = -1
      for icheck = itrace-4 > 0, (itrace+4) < (ntrace-1) do begin

      ;-------------------------------------------------------------------
      ; accept only good traces (xgood), which are not itself, and
      ; have less than 100 bad centroids returned from trace_fweight
      ;
         if xgood[icheck] AND icheck NE itrace AND $
               badpix[icheck] LT 100 then begin
           if checktrace[0] EQ -1 then checktrace = icheck $
           else checktrace = [checktrace, icheck]
         endif
      endfor

      ;-------------------------------------------------------------------
      ;  Checktrace contains the selected good neighbors
      ;
      ncheck = n_elements(checktrace)

      ;---------------------------------------------------------------------
      ;  Need at least two good neighbors to perform correction
      ;
      if ncheck GE 2 then begin
        clean = total(tmp_xerr[*,checktrace] EQ 999,2) EQ 0
        meantrace = total(xset[*,checktrace],2)/ncheck

        goodrows = where(clean AND tmp_xerr[*,itrace] NE 999, ngoodrows)

      ;---------------------------------------------------------------------
      ;  Require at least 100 common good rows to carry-on
      ;
        if ngoodrows GE 100 then begin
          xmask[goodrows,itrace] = 1
          xset_good = xset[goodrows,*]
          offset = tmp_xpos[goodrows,itrace] # replicate(1,ncheck) - $
                 xset_good[*,checktrace]
          shift =  mean(djs_median(offset,1))
          xset[*,itrace] = meantrace + shift
        endif else begin
           splog, 'Fiber '+ strtrim(fix(itrace)+1,2)+ ' Only '+ $
             strtrim(fix(ngoodrows),2)+ ' rows to adjust trace, skipping'
             xmask[itrace,*] = 0
        endelse

      endif else begin
           splog, 'Fiber '+ strtrim(fix(itrace)+1)+' Only '+ $
             strtrim(fix(ncheck),2)+' neighboring good fibers, skipping'
           xmask[itrace,*] = 0
      endelse

    endfor
 
   ;----------
   ; Replace XSET with a smooth trace-set

;   xy2traceset, yset, xset, tset, ncoeff=5, yfit=xnew, invvar=xmask, $
;    maxdev=maxdev, maxrej=1, /sticky
;   xset = xnew

   ;----------
   ; Set FIBERMASK bit for any fiber with more than 20% of its positions
   ; --->let's change it to 45%, there are far too many traces flagged as
   ;      BADTRACE with so many bad columns and a 3 pixel radius!) ???
   ; masked, which includes any positions off the CCD.  Do not pay any
   ; attention to XGOOD, since that may indicate that a fiber is only
   ; bad near YSTART.

   ibad = where(total(1-xmask, 1) GT 0.45*ny)
   if (ibad[0] NE -1) then $
    fibermask[ibad] = fibermask[ibad] OR fibermask_bits('BADTRACE')


   if tag_exist(fiberparam, 'deadfibermask') then begin
        positions = fiberparam.deadfibermask
        valid = where(positions ne -1, n_valid)
        if n_valid gt 0 then begin
            idead = positions[valid] - 1
            fibermask[idead] = fibermask[idead] OR fibermask_bits('NOPLUG')
        endif
    endif
   
   if tag_exist(fiberparam, 'bundlefibers') then begin
         bundlefibers = fiberparam[0].bundlefibers
   endif else bundlefibers = make_array(25,/integer,value=20)
   return, xset
end
;------------------------------------------------------------------------------
