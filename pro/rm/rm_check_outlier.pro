; Check rate of outlier using standard stars

pro rm_check_outlier, maxdev=maxdev

  if not keyword_set(maxdev) then maxdev=0.1 ; dex 

  target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
  fibermap = mrdfits(target_file,1,/silent)
  plate=fibermap[0].plate & mjd=fibermap[0].mjd
  ind=where(plate gt 0)
  plate=plate[ind] & mjd=mjd[ind]

  ; this is for the pipeline
  ttt=rm_findstd(plate,mjd,calibdir='',synflux=synflux_pipe,calibflux=calibflux_pipe)
  ; this is for the custom reduction
  ttt=rm_findstd(plate,mjd,calibdir='recalib/',synflux=synflux_new,calibflux=calibflux_new)

  ratio_pipe=synflux_pipe / calibflux_pipe
  ratio_new=synflux_new / calibflux_new

  ; keep only gri band
  ratio_pipe=ratio_pipe[*,*,1:3]
  ratio_new=ratio_new[*,*,1:3]

  dims=size(ratio_pipe,/dim)
  nep=dims[0] & nstd=dims[1] & nband=dims[2]

  ratio_pipe=reform(ratio_pipe, nep*nstd, nband)
  ratio_new=reform(ratio_new, nep*nstd, nband)

  ; remove bad spectra
  badflag = lonarr(nep*nstd)
  badflag1=badflag
  for i=0L, nep*nstd - 1 do begin
     ind=where(ratio_pipe[i, *] gt 0, nnn)
     if nnn ne nband then badflag[i]=1L
     ind=where(ratio_new[i, *] gt 0, nnn)
     if nnn ne nband then badflag1[i]=1L
  endfor

  ind=where(badflag eq 0, ngood)
  ind1=where(badflag1 eq 0, ngood1)
  print, where(ind ne ind1)
  ratio_pipe=ratio_pipe[ind, *]
  ratio_new=ratio_new[ind,*]

  ; get outlier rate
  logratio_pipe=alog10(ratio_pipe)
  logratio_new=alog10(ratio_new)

  ind=where(abs(logratio_pipe[*,0]) gt maxdev or $
      abs(logratio_pipe[*,1]) gt maxdev or $
      abs(logratio_pipe[*,2]) gt maxdev, nout_pipe)
  ind=where(abs(logratio_new[*,0]) gt maxdev or $
      abs(logratio_new[*,1]) gt maxdev or $
      abs(logratio_new[*,2]) gt maxdev, nout_new )

  print, 'Outlier frac; |dev|>', maxdev, ' dex'
  print, 'pipeline: ', float(nout_pipe)/float(ngood)
  print, 'new: ', float(nout_new)/float(ngood)

  print,  logratio_new[ind,*]

end
