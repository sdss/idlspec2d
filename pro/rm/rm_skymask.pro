;+
; NAME:
;
; PURPOSE:
;   Return a mask using David Kirby's skymask (dr9-sky-mask)
;   wh-sky-mask is the one generated using all sky spectra from the RM program by Yue Shen

function rm_skymask, wave, margin=margin, skyfile=skyfile, fmt=fmt

   if not keyword_Set(margin) then margin = 1L

   if ~keyword_Set(skyfile) then $
   ;skyfile = getenv('IDLRM_DIR')+'/etc/dr9-sky-mask.txt'
   skyfile=getenv('IDLRM_DIR')+'/etc/wh-sky-mask.txt'
   if ~keyword_set(fmt) then fmt='x,d' ; fmt='d' for dr9-sky-mask
   readcol,skyfile,format=fmt,/silent,skymask

   npix = n_elements(wave)
   mask = lonarr(npix)
   mask[*] = 1
   
   arr = abs( alog10( wave # (1./skymask) ) )*1d4
   dist_arr = round( min( arr, dim=2 ) ) 
   ind = where(dist_arr lt margin)
   if ind[0] ne -1 then mask[ind] = 0

   return, mask

end
