;+
; NAME:
;   opfiber_generate
;
; PURPOSE:
;   Generate or append entries to the opFibers.par file
;
; CALLING SEQUENCE:
;   opfiber_generate, plate, mjd=, [ camname=, expnum= ]
;
; INPUTS:
;   plate      - Plate number(s)
;   mjd        - MJD(s) for plate
;
; OPTIONAL INPUTS:
;   camname    - Camera name; default to all 4 cameras
;   expnum     - Exposure number (do not use if PLATE is an array);
;                default to first science exposure for this plate+MJD
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Append entries in the opFibers.par file, or create the file if it
;   does not yet exist.  This file used to establish the search parameters
;   for the trace-finding in TRACE320CRUDE.
;
;  The MJD for each entry is set to the observation, but should be changed
;  by hand to the earliest MJD for which the solution would be valid.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/opfiles/opFibers.par
;
; PROCEDURES CALLED:
;   spframe_read
;   splog
;   yanny_readone
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   27-Oct-2010  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro opfiber_generate, plate, mjd=mjd, expnum=expnum1, camname=camname1

   nplate = n_elements(plate)
   if (NOT keyword_set(plate)) then begin
      print, 'Must specify PLATE'
      return
   endif
   if (n_elements(mjd) NE nplate) then $
    message, 'PLATE and MJD must have same number of elements!'

   ; If PLATE is a vector, then call recursively
   if (nplate GT 1) then begin
      for iplate=0L, nplate-1L do begin
         opfiber_generate, plate[iplate], mjd=mjd[iplate], camname=camname1
      endfor
      return
   endif

   fiberparam_blank = create_struct(name='FIBERPARAM', $
    'CARTID'    ,   0, $
    'CAMNAME'   ,'  ', $
    'MJD'       ,  0L, $
    'XSTART_TOL',  0., $
    'NFIBER'    ,  0L, $
    'FIBERSPACE', fltarr(25), $
    'BUNDLEGAP' , fltarr(25) )
   fiberparam = replicate(fiberparam_blank, 4)

   if (keyword_set(camname)) then camname = camname1 $
    else camname = ['b1','b2','r1','r2']

   readspec, plate, mjd=mjd, objhdr=objhdr
   if (NOT keyword_set(objhdr)) then $
    message, 'spPlate file not found for plate '+string(plate)+ ' MJD '+string(mjd)
   idlist = sxpar(objhdr, 'EXPID*')
   cartid = sxpar(objhdr, 'CARTID')
   nfiber = sxpar(objhdr, 'NAXIS2')
   nhalf = nfiber / 2

   for icam=0, n_elements(camname)-1 do begin
      spectroid = long(strmid(camname[icam],1,1))
      camcolor = strmid(camname[icam],0,1)
      if (keyword_set(expnum1)) then begin
         expnum = expnum1
      endif else begin
         indx = (where(strmatch(idlist,'*'+camname[icam]+'*'), ct))[0]
         if (ct EQ 0) then message, 'No matching camera!'
         expnum = long(strmid(idlist[indx],3,8))
      endelse

      filename = 'spCFrame-'+camname[icam]+'-'+string(expnum,format='(i8.8)') $
       + '.fits'
      fullname = filepath(filename, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
       subdir=[getenv('RUN2D'), string(plate,format='(i4.4)')])
      spframe_read, fullname, ximg=ximg, hdr=hdr
      ; test that this is a valid XIMG, otherwise read from spFlat file
      if (max(ximg) LT 1000) then begin
         flatname = 'spFlat-'+strmid(sxpar(hdr,'FLATFILE'),4,11)+'.fits'
         fullflatname = filepath(flatname, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
          subdir=[getenv('RUN2D'), string(plate,format='(i4.4)')])
         fullflatname = (findfile(fullflatname+'*'))[0]
         if (keyword_set(fullflatname)) then begin
            tset = mrdfits(fullflatname,1)
            traceset2xy, tset, xx, ximg
         endif else ximg = 0
      endif

      if (NOT keyword_set(ximg)) then $
       message, 'File not found: '+fullname

      iy = (size(ximg,/dimens))[0] / 2
      xpos = reform(ximg[iy,*])
      xdiff = shift(xpos,-1) - xpos ; last value is then nonsense

      fiberparam[icam].cartid = cartid
      fiberparam[icam].camname = camname[icam]
      fiberparam[icam].mjd = mjd
      fiberparam[icam].xstart_tol = 70.
      fiberparam[icam].nfiber = nhalf
      nbundle = nhalf / 20
      for j=0, nbundle-1 do $
       fiberparam[icam].fiberspace[j] = median(xdiff[j*20:j*20+18])
      fiberparam[icam].bundlegap[0] = xpos[0] ; XSTART
      for j=1, nbundle-1 do $
       fiberparam[icam].bundlegap[j] = xpos[j*20] - xpos[j*20-1] $
        - fiberparam[icam].fiberspace[j]
   endfor

   ; Append this to the opFibers.par file
   opfile = filepath('opFibers.par', root_dir=getenv('IDLSPEC2D_DIR'), $
    subdir='opfiles')
   opdata = yanny_readone(opfile)
   if (keyword_set(opdata)) then begin
      splog, 'Appending to exiting file '+opfile
      opdata = struct_append(opdata, fiberparam)
   endif else begin
      splog,' Creating new file '+opfile
      opdata = fiberparam
   endelse
   yanny_write, opfile, ptr_new(opdata)

   return
end
;------------------------------------------------------------------------------
