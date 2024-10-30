

;+
; NAME:
;   spSpec2spFullSky
;
; PURPOSE:
;   Build spFullsky files (all sky spField proxy) for target level (independent of field-mjd) coadds
;
; CALLING SEQUENCE:
;
; INPUTS:
;   coadd - the name of the coadd schema
;
; OPTIONAL KEYWORDS:
;   topdir - the daily coadd base directory
;   mjd    - the characteristic MJD of the coadds to include
;   runmjd - The characteristic MJD of the run
;   outdir - the output directory of the files
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   coaddhdr - the meta header of spFullsky
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
;
;-
;------------------------------------------------------------------------------


pro spSpec2spFullSky, coadd, topdir=topdir, mjd=mjd, runmjd=runmjd, outdir=outdir, coaddhdr=coaddhdr

RESOLVE_ROUTINE,'sdss_maskbits',/EITHER,/SKIP_EXISTING, /quiet
RESOLVE_ALL, /SKIP_EXISTING, /quiet, /CONTINUE_ON_ERROR;, class='COMMON'
CPU, TPOOL_NTHREADS = 1

  if not keyword_set(topdir) then begin
      topdir=getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D')
  endif
  if not keyword_set(outdir) then begin
      outdir = filepath(strtrim(coadd,2),root_dir = topdir)
  endif
  

  
  if not keyword_set(mjd) then begin
    mjd = '*'
    spFieldname = filepath('spFullsky-'+strtrim(coadd,2)+'.fits',root_dir=get_field_dir(topdir,'',coadd,/custom))
    plotsnfile = filepath('spSN2d-'+strtrim(coadd,2)+'.ps',root_dir=get_field_dir(topdir,'',coadd,/custom))
    logfile = filepath('spSpec2spFullsky-'+strtrim(coadd,2)+'.log',root_dir=get_field_dir(topdir,'',coadd,/custom))
  endif else begin
    spFieldname = filepath('spFullsky-'+strtrim(coadd,2)+'-'+strtrim(mjd,2)+'.fits',root_dir=get_field_dir(topdir,'',coadd,/custom))
    plotsnfile = filepath('spSN2d-'+strtrim(coadd,2)+'-'+strtrim(mjd,2)+'.ps',root_dir=get_field_dir(topdir,'',coadd,/custom))
    logfile = filepath('spSpec2spFullsky-'+strtrim(coadd,2)+'-'+strtrim(mjd,2)+'.log',root_dir=get_field_dir(topdir,'',coadd,/custom))
  endelse
  print, spFieldname
  if not keyword_set(runmjd) then runmjd = mjd

  ;cpbackup, logfile
  ;splog, filename=logfile
  ;splog, 'Log file ' + logfile + ' opened ' + systime()

  splog, filepath('spSpec*.fits',root_dir = get_field_dir(topdir,'',coadd,/custom), subdirectory=['coadd', strtrim(mjd,2)])
  spSpecfiles = findfile(filepath('spSpec*.fits',root_dir = get_field_dir(topdir,'',coadd,/custom), $
                                    subdirectory=['coadd', strtrim(mjd,2)]), count = nspSpec)

  if nspSpec eq 0 then begin
    splog, 'No spSpec files found'
    return
  endif else splog, nspSpec, ' spSpec files found'
    
  FLUX    = 0
  IVAR    = 0
  ANDMASK = 0
  ORMASK  = 0
  WAVEDISP= 0
  PLUGMAP = 0
  SKY     = 0
  SPECRESL= 0
  binsz   = 0
  wavemin = 0
  foreach spf, spSpecfiles, idx do begin
        ;read spSpec file
    if idx mod 100 eq 0 then splog, idx+1 , ' of ',nspSpec, ' spSpec files read'
    coadd   = mrdfits(spf, 1,/silent)
    plugmap_temp = mrdfits(spf, 2,/silent)
    hdr = headfits(spf, /silent)
    plugmap_temp[0].TARGET_INDEX = idx+1
    
    if not keyword_set(FLUX) then begin
        FLUX     = coadd.FLUX
        IVAR     = coadd.IVAR
        ANDMASK  = coadd.AND_MASK
        ORMASK   = coadd.OR_MASK
        WAVEDISP = coadd.WDISP
        PLUGMAP  = plugmap_temp
        SKY      = coadd.SKY
        SPECRESL = coadd.WRESL
        binsz = sxpar(hdr, 'COEFF1')
        wavemin = sxpar(hdr, 'COEFF0')
    endif else begin
        FLUX     = [[FLUX],    [coadd.FLUX]]
        IVAR     = [[IVAR],    [coadd.IVAR]]
        ANDMASK  = [[ANDMASK], [coadd.AND_MASK]]
        ORMASK   = [[ORMASK],  [coadd.OR_MASK]]
        WAVEDISP = [[WAVEDISP],[coadd.WDISP]]
        PLUGMAP  = struct_append(PLUGMAP, plugmap_temp)
        SKY      = [[SKY],     [coadd.SKY]]
        SPECRESL = [[SPECRESL],[coadd.WRESL]]
        binsz    = [binsz,sxpar(hdr, 'COEFF1')]
        wavemin  = [wavemin,sxpar(hdr, 'COEFF0')]
    endelse


  endforeach
  splog, idx , ' of ',nspSpec, ' spSpec files read'
  ;splog,'wavemin', min(wavemin), max(wavemin)
  ;splog,'binsz',  binsz

  wavemin = max(wavemin)
  binsz = max(binsz)
  finalwave = dindgen(n_elements(coadd.FLUX)) * binsz[0] + wavemin[0]

  hdr0 =[' ']
  sxaddpar, hdr0, 'COEFF0', wavemin, 'Central wavelength (log10) of first pixel'
  sxaddpar, hdr0, 'COEFF1', binsz, 'Log10 dispersion per pixel'
  sxaddpar, hdr0, 'CRVAL1', wavemin, 'Central wavelength (log10) of first pixel'
  sxaddpar, hdr0, 'CD1_1', binsz, 'Log10 dispersion per pixel'
  sxaddpar, hdr0, 'CRPIX1', 1,'Starting pixel (1-indexed)'
  sxaddpar, hdr0, 'CTYPE1', 'LINEAR'


  ;sxaddpar, hdr0, 'CPADD
  sxaddpar, hdr0, 'MJD', mjd, 'Modified Julian Date'
  sxaddpar, hdr0, 'RUNMJD', runmjd, 'Modified Julian Date of RUN'
  sxaddpar, hdr0, 'RUN2D', sxpar(hdr,'RUN2D'), 'IDLSpec2D Run2d'
  sxaddpar, hdr0, 'TILEID', '', "Tile ID for SDSS BOSS plates (-1 for SDSS)"

    
  jdtemp=mjd+2400000.5
  mphase,jdtemp,mfrac
  sxaddpar, hdr0, 'MOON_FRAC',mfrac, 'Moon Phase', after='EXPTIME'

    

  
  platesn, FLUX, IVAR, ANDMASK, PLUGMAP, finalwave, obs='apo', hdr=hdr0, plotfile=plotsnfile
  mwrfits_named, FLUX, spFieldname, hdr=hdr0, name='FLUX', /create, /silent

  ; HDU #1 IVAR
  mkhdr, hdr1, IVAR, /image, /extend
  add_iraf_keywords, hdr1, wavemin, binsz
  sxaddpar, hdr1, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
  mwrfits_named, IVAR, spFieldname, hdr=hdr1, name='IVAR', desc=' Inverse Variance', /silent

  ; HDU #2 ANDMASK
  mkhdr, hdr2, ANDMASK, /image, /extend
  add_iraf_keywords, hdr2, wavemin, binsz
  mwrfits_named, ANDMASK, spFieldname, hdr=hdr2,name='ANDMASK', desc=' AND Mask', /silent
 
  ; HDU #3 ORMASK
  mkhdr, hdr3, ORMASK, /image, /extend
  add_iraf_keywords, hdr3, wavemin, binsz
  mwrfits_named, ORMASK, spFieldname, hdr=hdr3, name='ORMASK', desc=' OR Mask',/silent
 
  ; HDU #4 WAVEDISP
  mkhdr, hdr4, WAVEDISP, /image, /extend
  add_iraf_keywords, hdr4, wavemin, binsz
  sxaddpar, hdr4, 'BUNIT', 'pixels'
  mwrfits_named, WAVEDISP, spFieldname, hdr=hdr4, name='WAVEDISP', desc=' Wavelength dispersion', /silent

  ; HDU #5 PLUGMAP
  mwrfits_named, PLUGMAP, spFieldname, hdr=hdr5, name='PLUGMAP', desc=' Plugmap structure', /silent

  ; HDU #6 SKY
  mkhdr, hdr6, SKY, /image, /extend
  add_iraf_keywords, hdr6, wavemin, binsz
  mwrfits_named, SKY, spFieldname, hdr=hdr6, name='SKY', desc=' Subtracted sky flux',/silent

  ; HDU #7 SPECRESL
  mkhdr, hdr7, SPECRESL, /image, /extend
  add_iraf_keywords, hdr7, wavemin, binsz
  sxaddpar, hdr7, 'BUNIT', 'angstroms'
  mwrfits_named, SPECRESL, spFieldname, hdr=hdr7, name='SPECRESL', desc=' Spectral resolution',/silent

  coaddhdr = hdr0
  splog, 'Successful completion of spSpec2spField at ' + systime()
  ;splog, /close
end
