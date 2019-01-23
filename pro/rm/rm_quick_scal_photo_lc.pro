; This is a quick-and-dirty version to normalize the 
; difference img lc and the spec-synthetic lc using
; the sdss psf mag for each RM object

pro rm_quick_scal_photo_lc

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  epoch_info=mrdfits(file,2)

  ; diff LC dir
  topdir='/data3/yshen/ftp/sdssrm/collab/photo_lc/img_diff_lc/'

  ; output file -- by rmid
  outfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/final/'
