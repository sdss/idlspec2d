PRO runsubtract

DIR = ''                                 ;where your fits data files are
DIR2 = 'newfiles/'                       ;where you want the new fits files

galspc = ['0586/spSpec-52023-0586-324.fit','0385/spSpec-51877-0385-449.fit',$ ;example Galaxies
         '0507/spSpec-52353-0507-399.fit','0412/spSpec-52258-0412-312.fit',$
         '0581/spSpec-52356-0581-289.fit']
; galspc = ['0434/spSpec-51885-0434-177.fit','0360/spSpec-51816-0360-495.fit',$ ;example QSOs
;          '0499/spSpec-51988-0499-059.fit','0268/spSpec-51633-0268-160.fit',$
;          '0451/spSpec-51908-0451-174.fit']

ngal = n_elements(galspc)
plate = uintarr(ngal)
for i=0L,ngal-1 do plate[i] = uint((strsplit(galspc[i],'-',/extract))[2]) ;plate numbers

;Postscript plotting:
set_plot,'ps'
device,file='runsubtract.ps',/color,/portrait,xoffset=2,yoffset=2,ysize=25,xsize=18
!p.multi=[0,1,4]
for i=0,ngal-1 do begin

    ;read in data files and wave array
    data = readfits(dir+galspc[i],header,/silent)
    min_l = sxpar(header,'coeff0') ;central wavelength (log10) of 1st pix
    disp = sxpar(header, 'coeff1') ;dispersion per pixel
    npix  = sxpar(header, 'naxis1') ;no. of pixels
    wave = (dindgen(npix)*disp +min_l) ;wave array
    z = sxpar(header,'z')       ;redshift of galaxy
    eclass = sxpar(header,'eclass') ;eigenclass 

    final_spec = subtractoh(data[*,0], data[*,2], wave, z, plate[i], rms=rms, nrecon=nrecon,$
                                        eclass=eclass, /plotspec)
;;;now do what you want with final_spec....
;e.g. make new fits files

    galspc2 = (strsplit(galspc[i],'/',/extract))[1]
    file_copy,DIR+galspc[i],DIR2+galspc2
    data[*,0]=final_spec     ;replace spectrum with sky-subtracted version

; Suggested modifications to fits header
    SXADDPAR, header, 'SKYVAR0',rms[0],'Non-Sky pixel variance'
    SXADDPAR, header, 'SKYVAR1',rms[1],'Bad-Sky pixel variance (before skysub)'
    SXADDPAR, header, 'SKYVAR2',rms[2],'Bad-Sky pixel variance (after skysub)'
    SXADDPAR, header, 'NRECON',nrecon,'Number of components in reconstruction'

    modfits,DIR2+galspc2,data,header    ;modfits must be the NEW IDLastro version
    

endfor

cleanplot,/silent
device,/close
set_plot,'x'

END
