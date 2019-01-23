;
;+
; NAME:
;	subtractOH.pro

;
; PURPOSE:
;	Remove OH sky signatures beyond 6700\AA from Sloan Digital Sky
;	Survey spectra using a Principal Component Analysis.
;
;
; CALLING SEQUENCE:
;
;	newflux = SUBTRACTOH(Flux, Error, Wave, Z, Plate, Rms=Rms,
;	            Nrecon=Nrecon, [...optional input, keywords])
;
;
; INPUTS:
;	Flux:	SDSS flux array, e.g. data[*,0]
;       Error:  Error array, e.g. data[*,2]  
;       Wave:   Wavelength array in log(angstroms),
;               e.g. (dindgen(npixels)*dispersion + wave_pix1)
;       Z:      Redshift of spectrum
;       Plate:  SDSS plate number
;         
; OPTIONAL INPUTS:
;       Eclass:   If the object is a galaxy, eclass must be provided.
;                 Eclass<-0.05 are classified absorption line
;                 galaxies. Eclass>-0.05 are classified emission line
;                 galaxies (e.g. 'eclass' parameter in fits
;                 headers).
; 
;       Filter:   Specify choice if filter scale, default 55
;                 pixels. Ensure that scale large enough to be
;                 uneffected by OH lines.
;
;	Linefile_prim: Filename in which a list of lines can be found which
;	               should be masked during the procedure. Follow format
;	               of supplied line mask files.
;
;       Linefile_sec: Second filename containing list of lines of a
;                     secondary object in the spectrum
;                     (e.g. DLA). NOTE optional input 'zsec' must also
;                     be given.
;
;       Weights:  Produced by program getweights.pro to give
;                 approximate Poisson errors from SDSS noise array. It
;                 is recommended to pass these to the program if
;                 multiple files are being calculated to greatly
;                 improve speed.
;
;       Zsec:     Redshift of a secondary object in the spectrum. NOTE
;                 optional input 'linefile_sec' must also be given.
;
; KEYWORD PARAMETERS:
;       CB1:    Set this keyword to indicate the use of a common block
;               to pass the eigenspectra array into the program. This
;               improves speed greatly if many files are being
;               computed at once. 
;               common subtractOH1, espec, lambda, pix_sky, pix_nosky

;
;       EXTGAL: Set this keyword to indicate that this object is an
;               extreme emission line galaxy. A different line list
;               will be used to mask the emission lines. Recommended
;               criteria for extreme emission line galaxies:
;               EW(Ha)>200,z<0.4; EW(OIII)>200,z<0.84 
;       
;       PLOTSPEC: Set this keyword to cause the function to produce
;                 plots of the sky-residual subtraction process. There
;                 must be a suitable plotting device open!
;
;       QSO:    Set this keyword to cause median filtering to be
;               performed after line masking (better for broad
;               emission lines), and QSO line file to be
;               used if optional input 'Linefile' is not
;               specified. The default is a galaxy.
;
;       STAR:   Set this keyword to indicate that this object is a
;               star and the stellar line file is to be
;               used if optional input 'Linefile' is not
;               specified. The default is a galaxy.
;
;       SILENT: Set this keyword to stop the function reporting anything.
;
; OUTPUTS:
;	This function returns the sky subtracted SDSS spectrum.
;        
;
; OPTIONAL OUTPUTS:
;       RMS:    An array containing the rms values of the spectrum
;               [non-sky pixels, sky pixels before, sky pixels after]
;       NRECON: The number of components used during the reconstruction.
;	
; NOTE:
;       This routine has been slightly modified to use new eigenvectors
; created for the RM program (BOSS spectra); the sky residual eigenvectors 
; were created on the set of
; sky spectra unnormalized by the scaled errors, which is likely different
; from Wild&Hewett in which the eigenvectors were created on the normalized
; sky spectra. However, the two approaches give negligible difference, so
; I am sticking to my approach.

;----------------------------------------------------------------------------------------
;;  EXAMPLE PROGRAM TO MAKE NEW FITS FILES WITH SKY-RESIDUAL SUBTRACTED SPECTRA:
;
;; SET DIR and DIR2 before running.
;
;; FOR QSOS remember to set QSO keyword, and don't use ECLASS
;
; PRO runsubtract
;
; DIR = ''                                 ;where your fits data files are
; DIR2 = 'newfiles/'                       ;where you want the new fits files
;
; galspc = ['0586/spSpec-52023-0586-324.fit','0385/spSpec-51877-0385-449.fit',$ ;example Galaxies
;          '0507/spSpec-52353-0507-399.fit','0412/spSpec-52258-0412-312.fit',$
;          '0581/spSpec-52356-0581-289.fit']
;; galspc = ['0434/spSpec-51885-0434-177.fit','0360/spSpec-51816-0360-495.fit',$ ;example QSOs
;;          '0499/spSpec-51988-0499-059.fit','0268/spSpec-51633-0268-160.fit',$
;;          '0451/spSpec-51908-0451-174.fit']
;
; ngal = n_elements(galspc)
; plate = uintarr(ngal)
; for i=0L,ngal-1 do plate[i] = uint((strsplit(galspc[i],'-',/extract))[2]) ;plate numbers
;
; ;Postscript plotting:
; set_plot,'ps'
; device,file='runsubtract.ps',/color,/portrait,xoffset=2,yoffset=2,ysize=25,xsize=18
; !p.multi=[0,1,4]
; for i=0,ngal-1 do begin
; 
;     ;read in data files and wave array
;     data = readfits(dir+galspc[i],header,/silent)
;     min_l = sxpar(header,'coeff0') ;central wavelength (log10) of 1st pix
;     disp = sxpar(header, 'coeff1') ;dispersion per pixel
;     npix  = sxpar(header, 'naxis1') ;no. of pixels
;     wave = (dindgen(npix)*disp +min_l) ;wave array
;     z = sxpar(header,'z')       ;redshift of galaxy
;     eclass = sxpar(header,'eclass') ;eigenclass 
;
;     final_spec = subtractoh(data[*,0], data[*,2], wave, z, plate[i], rms=rms, nrecon=nrecon,$
;                                         eclass=eclass, /plotspec)
; ;;;now do what you want with final_spec....
; ;e.g. make new fits files
; 
;     galspc2 = (strsplit(galspc[i],'/',/extract))[1]
;     file_copy,DIR+galspc[i],DIR2+galspc2
;     data[*,0]=final_spec     ;replace spectrum with sky-subtracted version
;
; ; Suggested modifications to fits header
;     SXADDPAR, header, 'SKYVAR0',rms[0],'Non-Sky pixel variance'
;     SXADDPAR, header, 'SKYVAR1',rms[1],'Bad-Sky pixel variance (before skysub)'
;     SXADDPAR, header, 'SKYVAR2',rms[2],'Bad-Sky pixel variance (after skysub)'
;     SXADDPAR, header, 'NRECON',nrecon,'Number of components in reconstruction'
; 
;     modfits,DIR2+galspc2,data,header    ;modfits must be the NEW IDLastro version
;     
;
; endfor
;
; cleanplot,/silent
; device,/close
; set_plot,'x'
;
; END
;----------------------------------------------------------------------------------------
;
; MODIFICATION HISTORY: Vivienne Wild, vw@ast.cam.ac.uk, 22/12/04
;     29/01/05 Included RMS2 and RMS5 functions which were accidentally removed
;     05/06/05 changed keyword linefile to linefile_prim
;              set stellar line file to be same as absorption galaxy line file
;
;****************************************************************************************
FUNCTION rms2,x
;calculate variance
RETURN, SQRT((MOMENT(x,/nan))[1])
END

FUNCTION rms5,x
;calculate 67th percentile
y = MEAN(x,/nan)
index = SORT(ABS(x-y))
a = N_ELEMENTS(x)*.67
percentile = ABS((x-y)[index[a]])

RETURN,percentile
END

;****************************************************************************************

FUNCTION subtractOH, Flux, Error, Wave, Z, Plate, Rms=Rms, Nrecon=Nrecon, $
                     Eclass=eclass, filter=filter, linefile_prim=linefile, linefile_sec = linefile_sec, zsec=zsec, weights=weights_in, $
                     CB1=CB1,QSO=QSO, STAR=STAR, EXTGAL=EXTGAL, sky=sky, $
            PLOTSPEC=PLOTSPEC, SILENT=SILENT, maxnrecon=maxnrecon

;**** Set this string to the path of the data files 
;DIR = '~/idl/pro/SDSS/SKY/PUBLIC/'
DIR=getenv('IDLRM_DIR') + '/pro/vwsky/'

;------------------------------------ CHECKING INPUTS ---------------------------
if (SIZE(flux,/dim))[0] ne (SIZE(wave,/dim))[0] and (SIZE(flux,/dim))[0] ne (SIZE(error,/dim))[0] then $
  MESSAGE,'FLUX and WAVE must be arrays of same size'

;what type of object?
if KEYWORD_SET(qso) then QSO=1 else QSO=0 ;QSO
if KEYWORD_SET(star) then STAR=1 else STAR=0 ;stellar spectrum
if KEYWORD_SET(extgal) then EXTGAL=1 else EXTGAL=0 ;extreme emission line galaxy

if N_ELEMENTS(eclass) eq 0 and NOT(QSO) and NOT(STAR) and NOT(SKY) then $
  MESSAGE, 'A galaxy needs an ECLASS: <-0.05 absorption line, >0.05 emission line' 
if (QSO) or (STAR) or (SKY) then eclass = -9999

if NOT(QSO) and NOT(STAR) and NOT(EXTGAL) and eclass lt -0.05 then ABSGAL=1 else ABSGAL=0 ;abs line galaxy
if NOT(QSO) and NOT(STAR) and NOT(EXTGAL) and eclass ge -0.05 then EMGAL=1 else EMGAL=0 ;em line galaxy

;if (QSO) then maxnrecon = 200 else maxnrecon = 150
if not keyword_set(maxnrecon) then begin
  if (QSO) then maxnrecon = 500 else maxnrecon = 200
endif

;bits and pieces
if N_ELEMENTS(filter) eq 0 then filter = 55 ;size of median filter in pixels.
if KEYWORD_SET(plotspec) then plotspec = 1 else plotspec=0 ;to plot or not
if KEYWORD_SET(silent) then silent=1 else silent=0 ;to make comments or not
if n_elements(CB1) eq 0 then CB1=1 ;common block for eigenspectra array

if N_ELEMENTS(zsec) ne 0 then SEC=1 else SEC=0 ;Secondary object in spectrum
if (SEC) and N_ELEMENTS(linefile_sec) eq 0 then message,'Please provide filename for file containing lines of secondary object'

;**** Read in eigen-spectra file
; contains parameters: espec, lambda, pix_sky, pix_nosky
; ESPEC = fltarr(no. espectra, no. pixels)
; LAMBDA = fltarr(no. pixels)
; PIX_SKY = lonarr to identify those pixels with sky signatures to correct
; PIX_NOSKY = lonarr to identify those pixels used as control

;if CB1 then $
;  common subtractOH1, espec, lambda, pix_sky, pix_nosky $
;else RESTORE, DIR+'espec_OH.sav'
common subtractOH1, espec, lambda, pix_sky, pix_nosky

; populate the common block
if n_elements(espec) eq 0 then begin
  ;file='/data3/quasar/yshen/ftp/bossredux/v5_7_1/wh_skysub/pca.fits'
  file=getenv('IDLRM_DIR')+'/template/sky_residual_pca.fits'

  pca=mrdfits(file,1)
  lambda=pca.loglam & espec=pca.eigenspec & pix_sky=pca.ind_sky & pix_nosky=pca.ind_nonsky
  ; cull out a red-wavelength region
  ind_tmp=where(10.^lambda gt 5000., ncomplement=ncomp)
  lambda=lambda[ind_tmp]
  pix_sky=pix_sky - ncomp & pix_nosky=pix_nosky - ncomp
  ; remove the nosky pixels below 5000 A
  ind_tmp=where(pix_nosky ge 0)
  pix_nosky=pix_nosky[ind_tmp]
endif

npix = N_ELEMENTS(lambda)       ;number of pixels we are working with
disp = 0.0001                   ;SDSS spectrum dispersion

;**** Find pixels we are correcting
red_ind = WHERE(wave ge lambda[0]-disp/3. and wave le lambda[npix-1]+disp/3.,count)

if NOT(count) then MESSAGE, 'Input wavelength array not right - observed frame? log(wavelength)?'

if count lt npix then begin     ;sometimes a few red pixels are missing 
    WAVEFLAG = 1
    flux_in = [flux,fltarr(npix-count)] ;fill gap with zeros
    error_in = [error,fltarr(npix-count)] ;set error=0 so not used in anything
    wave_in = [wave,lambda[count:npix-1]] ;fill up wavelength array
    red_ind = WHERE(wave ge lambda[0]-disp/3. and wave le lambda[npix-1]+disp/3.) ;reset red indexing
endif else begin
    WAVEFLAG = 0
    flux_in = flux
    error_in = error
    wave_in = wave
endelse

wave_red = wave_in[red_ind]


;-------------------------- PREPARE SPECTRUM -----------------------------------------
;**** Median filter the spectrum
;pass the whole flux array to allow correct filtering in blue end
;**** Create pixel-mask to ID emission/absorption features
;file containing wavelengths of potential line features to be masked
if N_ELEMENTS(linefile) eq 0 then begin
    if QSO then linefile = DIR+'qso_lines.dat'
    if STAR then begin
        linefile = DIR+'absgal_lines.dat'
        if NOT(SILENT) then print,'star: currently using absorption line galaxy line list'
    endif
    if ABSGAL then linefile = DIR+'absgal_lines.dat'
    if EMGAL then linefile = DIR+'emgal_lines.dat'
    if EXTGAL then linefile = DIR+'extgal_lines.dat'
endif

if NOT(QSO) then begin          ;GAL, STAR: median filter after masking
   
    flux2 = 1
    if not keyword_set(sky) then begin
      pix_mask = MASKLINES(flux_in,(10^wave_in)/(1+z), linefile, flux2, silent=silent)
    endif else begin
      flux2=flux_in
      pix_mask = lonarr(n_elements(flux2))
    endelse

    medflux = MEDFILT(flux2,red_ind,filter) ;median filter using flux with abs features masked
    flux_med = flux2[red_ind]-medflux
    pix_mask = pix_mask[red_ind]

endif else begin                ;QSO: median filter before masking

    if NOT(silent) then PRINT, 'SUBTRACTOH: masking QSO lines'

    medflux = MEDFILT(flux_in,red_ind,filter)
    flux_med = flux_in[red_ind]-medflux
    
    pix_mask = MASKLINESQSO(flux_med,(10^wave_red), z, linefile, filter, silent=silent) ;max mask size set by filter

endelse

if SEC then begin               ;secondary object in spectrum (e.g. DLA)
    if NOT(silent) then PRINT, 'SUBTRACTOH: masking secondary lines'
    pix_mask2 = MASKLINES(flux_med, (10^wave_red)/(1+zsec),linefile_sec, silent=silent) 
    pix_mask = pix_mask+pix_mask2 < 1 ; < "minimum operator"
endif

;**** ERROR=0 pixels - for some reason flux is undefined here
;add these pixels into pixel mask
ind = WHERE(error_in[red_ind] eq 0,count)
if count ne 0 then pix_mask[ind] = 1

;**** Create pixel arrays of sky/non-sky pixels, containing only pixels without line features
pix_sky_good = pix_sky
pix_espec = WHERE(NOT(pix_mask[pix_sky]))

if pix_espec[0] eq -1 then begin
    rms = fltarr(3)
    nrecon = 0
    print, 'problem with spectrum, returning'
    return,flux
endif

pix_sky_good = pix_sky_good[pix_espec]

pix_nosky_good = pix_nosky
ind = WHERE(NOT(pix_mask[pix_nosky]))

if ind[0] eq -1 then begin
    rms = fltarr(3)
    nrecon = 0
    print, 'problem with spectrum, returning'
    return,flux
endif

pix_nosky_good = pix_nosky_good[ind]

;-------------------------- RECONSTRUCT SKY FEATURES ------------------------------------
;**** Divide flux by mean plate error, weighted to correct SDSS noise.
;if N_ELEMENTS(weights_in) ne 0 then weights = weights_in else weights = GETWEIGHTS(plate, DIR)
; Note that the weights_in is in the same dimension as the input spectrum
weights = weights_in[red_ind] ; this is the reverse of the rescaled plate error array

;; check if any of the to-be-corrected sky pixels have bad weights (INFINITE)
;ttt=where(finite(weights[pix_sky_good]) eq 0, n_bad_weights)
;print, '# of bad weights: ', n_bad_weights
flux_wgt = flux_med*weights

;**** Project prepared spectrum onto sky eigenspectra
pcs = flux_wgt[pix_sky_good] ## (TRANSPOSE(espec[pix_espec,0:maxnrecon-1]))

;**** Calculate reconstruction based on rms of sky and non-sky pixels
RMS_nosky = RMS5(flux_wgt[pix_nosky_good])
RMS_b4sky = RMS2(flux_wgt[pix_sky_good])

recon = FLTARR(N_ELEMENTS(pix_sky))
final_spec = flux_in            ;outgoing spectrum = ingoing spectrum for now
nrecon = 0                      ;minimum no. of components
ii=0

while(1) do begin

    if ii eq 0 then RMS_afsky = RMS_b4sky else RMS_afsky = RMS2(flux_sub[pix_sky_good]) ;RMS of sky pixels

    if RMS_afsky lt RMS_nosky or ii eq maxnrecon then begin
        nozero = WHERE(weights[pix_sky] ne 0)
        final_spec[red_ind[pix_sky[nozero]]] = $
          flux_in[red_ind[pix_sky[nozero]]] - recon[nozero]/weights[pix_sky[nozero]] ;correct outgoing spectrum
        nrecon = ii             ;final number of components
        break
    endif

;if RMS not satisfied then carry on reconstructing
    recon = recon+pcs[ii]*espec[*,ii]

    flux_sub = flux_wgt
    flux_sub[pix_sky] = flux_wgt[pix_sky]-recon
    ii = ii+1

endwhile

rms = [RMS_nosky, RMS_b4sky, RMS_afsky]

if nrecon eq maxnrecon and NOT(silent) then $
  PRINT, 'SUBTRACTOH WARNING: reached max no. of components set. Recommend inspecting output for problems.'

;-------------------- MAKE SOME PLOTS ------------------------------------------------------
if (PLOTSPEC) then begin
    ;LOADCT, 23,ncolors=30,silent=silent

    !P.multi=[0,1,2]
    yrange=[-5, max(medflux)*1.5]
    xrange=[min(10.^wave_red), max(10.^wave_red)]
    PLOT, 10^wave_red, flux_in[red_ind],/ynozero,xstyle=1,title='black: before; red: after; blue: median filter '+string(rms,form='(3F7.3)'), yrange=yrange, xrange=xrange
    OPLOT, 10^wave_red, final_spec[red_ind],color=cgcolor('red')
    OPLOT, 10^wave_red, medflux,color=cgcolor('blue')

    ind = WHERE(pix_mask eq 1,count)
    if count ne 0 then OPLOT, 10^wave_red[ind],fltarr(count)+yrange[0], psym=3,color=cgcolor('cyan')


    PLOT, 10^wave_red[pix_sky_good], flux_med[pix_sky_good], xstyle=1,title='black: before, median filtered; red: reconstruction '+string(nrecon,form='(I3)'),psym=5, xrange=xrange
    OPLOT, 10^wave_red[pix_sky], recon/weights[pix_sky],color=cgcolor('red'),psym=5
    if count ne 0 then OPLOT, 10^wave_red[ind],fltarr(count), psym=3,color=cgcolor('cyan')
    !P.multi=0

endif


if WAVEFLAG then final_spec = final_spec[0:N_ELEMENTS(flux)-1] ;remove added on pixels

;message,'stop and diag'

; set the ivar=0 pixels to the original flux
ind=where(error eq 0)
if ind[0] ne -1 then final_spec[ind] = flux[ind]

return, final_spec

END
