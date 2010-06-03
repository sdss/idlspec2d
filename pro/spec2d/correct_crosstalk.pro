;+
; NAME:
; CORRECT_CROSSTALK
;
; PURPOSE:
; If there is keyword_set(image) is true this procedure corrects an image for crosstalk.  The input image is overwritten by the corrected image.
; If there is keyword_set(invvar) is true this procedure returns the diagonal of the resulting covariance matrix of the crosstalk-corrected
; image.  The input invvar is overwritten by the corrected invvar.  For now covariance matrices are not being used and the current implementation
; assumes this.
; 
; The amplifier layout and crosstalk terms are provided.
;
; CALLING SEQUENCE:
; CORRECT_CROSSTALK, Image, Crosstalk, Ampgrid
;
; INPUTS:
; Image -   2-dimensional array containing the image to be corrected. Note that this array is overwritten with the OUTPUT
; Crosstalk - An N by N grid containing the amplifier cross-talk.  N is the total number of amplifiers.  If NOT keyword_set(Crosstalk)
; then no correction is applied.
; Ampgrid - A 2-element array describing how the amplifiers are distributed on the chip.  For now, the code only
; knows how to deal with [2,2]
;
; OPTIONAL INPUTS:
; Invvar - 2-dimensional array containing the inverse variance that corresponds with the input image. This array is
; overwritten with the output
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
; Image - The cross-talk corrected image.  The input image is overwritten
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
; The input image is overwritten
;
; EXAMPLE:
;  infile='/Users/akim/Work/sdss3/testData/realdata/rawdata/55089/sdR-r1-00101369.fit'
;  sdssproc, infile, image, indir=indir, /applybias, /silent
;  configuration=obj_new('configuration', 55120)
;  correct_crosstalk,image, configuration->crosstalk_crosstalk('r1'), configuration->crosstalk_ampgrid('r1')
;
; MODIFICATION HISTORY:
;   Written by: Alex Kim, Oct 16, 2009
;-


pro correct_crosstalk, image, crosstalk, ampgrid, invvar

  if (NOT keyword_set(crosstalk)) then return
  
  if (ampgrid[0] NE 2 or ampgrid[1] NE 2) then begin
    splog, 'ABORT: Not implemented for this amplifier grid'
    return
  endif
  
  invxtalk = invert(crosstalk)
  
  if (keyword_set(image)) then begin
    ;
    ; Break apart the image
    ;
  
    ;determine the size of the image
    dim=size(image,/dim)
    
    ;get each amplifier readout and flip to same readout order
    im0=image[0:dim[0]/2-1,0:dim[1]/2-1]
    im1=rotate(image[dim[0]/2:*,0:dim[1]/2-1],5)
    im2=rotate(image[0:dim[0]/2-1,dim[1]/2:*],7)
    im3=rotate(image[dim[0]/2:*,dim[1]/2:*],2)
    
    ;containers for the answer
    
    im0p=im0
    im1p=im1
    im2p=im2
    im3p=im3
    
    ; do every pixel
    for i=0l,n_elements(im0)-1 do begin
      yvec=[im0[i],im1[i],im2[i],im3[i]]
      ans=invxtalk##yvec
      im0p[i]=ans[0]
      im1p[i]=ans[1]
      im2p[i]=ans[2]
      im3p[i]=ans[3]
    endfor
    
    ;flip back
    im1p=rotate(temporary(im1p),5)
    im2p=rotate(temporary(im2p),7)
    im3p=rotate(temporary(im3p),2)
    
    image[0:dim[0]/2-1,0:dim[1]/2-1]=im0p
    image[dim[0]/2:*,0:dim[1]/2-1]=im1p
    image[0:dim[0]/2-1,dim[1]/2:*]=im2p
    image[dim[0]/2:*,dim[1]/2:*]=im3p
  endif
  
  if (keyword_set(invvar)) then begin
  
    invxtalktrans=transpose(invxtalk)
    
    ;
    ; Break apart the image
    ;
    
    ;determine the size of the image
    dim=size(invvar,/dim)
    
    ;get each amplifier readout and flip to same readout order
    im0=invvar[0:dim[0]/2-1,0:dim[1]/2-1]
    im1=rotate(invvar[dim[0]/2:*,0:dim[1]/2-1],5)
    im2=rotate(invvar[0:dim[0]/2-1,dim[1]/2:*],7)
    im3=rotate(invvar[dim[0]/2:*,dim[1]/2:*],2)
    
    ;containers for the answer
    
    im0p=im0
    im1p=im1
    im2p=im2
    im3p=im3
    
    ; do every pixel
    for i=0l,n_elements(im0)-1 do begin
      v=fltarr(4,4)
      v(indgen(4),indgen(4))=[im0[i],im1[i],im2[i],im3[i]]
      w=where(v(indgen(4),indgen(4)) eq 0,neq0)
      if (neq0 gt 0) then v[w,w]=1/65535.
      v=invert(v)
      ans = invxtalk ## v ## invxtalktrans
      ans = invert(ans)
      
      ;this works for now while we are only dealing with diagonals
      if (neq0 gt 0) then ans[w,w]=0
      im0p[i]=ans[0,0]
      im1p[i]=ans[1,1]
      im2p[i]=ans[2,2]
      im3p[i]=ans[3,3]
    endfor
    
    ;flip back
    im1p=rotate(temporary(im1p),5)
    im2p=rotate(temporary(im2p),7)
    im3p=rotate(temporary(im3p),2)
    
    invvar[0:dim[0]/2-1,0:dim[1]/2-1]=im0p
    invvar[dim[0]/2:*,0:dim[1]/2-1]=im1p
    invvar[0:dim[0]/2-1,dim[1]/2:*]=im2p
    invvar[dim[0]/2:*,dim[1]/2:*]=im3p
  endif

end
