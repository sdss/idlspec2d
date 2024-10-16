;+
; NOTE: reprocess the single spectra spSepc files
; NAME:
;   reformat_spec
;
; PURPOSE:
;   Reformat spectra into a single fits file per object, combining all necessary
;   pieces from spField, spZall, spZline and spZbest
;
;   For each object, there is one file per target,
;   containing the finall coadded spectrum, and each blue-red merged single spectrum per exposure
;   frames. The file contains the next headers:
;
;     HDU 0  : Header info
;     HDU 1  : Coadded spectrum
;     HDU 2  : Summary metadata copied from spZbest
;     HDU 3  : Summary metadata copied from spZall
;     HDU 4  : Line fitting metadata from spZline
;     HDU 5  : Individual frame spectra per exposure
;
;   The format of each spec file is:
;
;      HDU 0 :
;      Header : from input spField with additional keywords from specObj:
;
;      Keyword     specObjColumn   Comment
;      -------     -----------   -------
;      [*] PLUG_RA     PLUG_RA       RA of object [deg]
;      [*] PLUG_DEC    PLUG_DEC      dec of object [deg]
;      TARGET_ID       TARGET_ID      Unique object identifier
;
;      [*] Note that RA and DEC already exist in the spPlate headers
;      but they are the telescope boresite RA and DEC, not the
;      RA and DEC of the object.
;
;      Data   : None
;      HDU 1 : Coadded Spectrum
;      Header : Minimum to define table
;      Data : Binary table with columns taken from the original spSpec and spZbest:
;      flux      : flux in 10^-17 ergs/s/cm^2/A
;      loglam    : log10(lambda[A])
;      ivar      : inverse variance
;      and_mask  : mask bits which affect every spectrum in coadd
;      or_mask   : mask bits which affect at least one spectrum in coadd
;      wdisp     : wavelength dispersion in dloglam units
;      sky       : subtracted sky flux in 10^-17 ergs/s/cm^2/A
;      wresl     : spectral resolution in A units
;      model     : best fit model for classification & redshift (from spZbest)
;
;      HDU 2 : Copy of row for this object from spZbest table + plugmap info of spFrame
;      SPECTRO_REDUX/RUN2D/FIELD/RUN1D/spZbest*.fits
;
;      HDU 3 : Copy of rows for this object from spZall table
;      SPECTRO_REDUX/RUN2D/FIELD/RUN1D/spZline*.fits
;
;      HDU 4 : Copy of rows for this object from spZline table
;      SPECTRO_REDUX/RUN2D/FIELD/RUN1D/spZline*.fits
;
;      HDU 5 .. +n_exp : Individual frames.
;      For each exposure there is one HDU.
;      These are in the order of the FPS confiuration_id, and
;
;      Header : Minimum to define table
;      Data: Binary table with columns taken from the original spSpec:
;      flux       : flux in 10^-17 ergs/s/cm^2/A
;      loglam     : log10(lambda[A])
;      ivar       : inverse variance
;      and_mask   : mask bits which affect every spectrum
;      or_mask    : mask bits which affect at least one spectrum
;      wdisp      : wavelength dispersion in dloglam units
;      sky        : subtracted sky flux in 10^-17 ergs/s/cm^2/A
;      wresl      : spectral resolution in A units
;
; CALLING SEQUENCE:
;   reformat_spec, [ platefile, fiberid=, run1d=, /doplot, /debug, chop_data= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fieldfile  - Field file(s) from spectro-2D; default to all files
;                matching 'spField*.fits'
;   run1d      - Optional override value for the environment variable $RUN1D
;   doplot     - If set, then generate plots.  Send plots to a PostScript
;                file spDiagDebug1d-$PLATE-$MJD.ps unless /DEBUG is set.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Input files are read from the current directory.
;   Output files are written to the subdirectory $RUN2D/SPECTRA/$MJD.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   The procedure requires these files:
;   ZALLFILE = 'spZall-Field-MJD.fits'
;   ZBESTFILE = 'spZbest-Field-MJD.fits'
;   ZLINEFILE = 'spZline-Field-MJD.fits'
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Jan-2020  Written by Hector Ibarra at UIUC
;------------------------------------------------------------------------------
;
pro struct_delete_field, struct, tag
  ;Delete an existing field from a structure.
  ;
  ;Inputs:
  ; tag (string) Case insensitive tag name describing structure field to
  ;  delete. Leading and trailing spaces will be ignored. If the requested
  ;  field does not exist, the structure is returned unchanged and without
  ;  error.
  ;
  ;Input/Output:
  ; struct (structure) structure to be modified.
  ;
  ;Examples:
  ;
  ; Delete sme.wave from structure:
  ;
  ;   IDL> struct_delete_field, sme, 'wave'
  ;
  ;History:
  ; 2003-Jul-26 Valenti  Adapted from struct_replace_field.pro.
  ; 2020-Jan-27 HJIM  Adapted for bhm.
  if n_params() lt 2 then begin
    print, 'syntax: struct_delete_field, struct, tag'
    return
  endif

  ;Check that input is a structure.
  if size(struct, /tname) ne 'STRUCT' then begin
    message, 'first argument is not a structure'
  endif

  ;Get list of structure tags.
  tags = tag_names(struct)
  ntags = n_elements(tags)

  ;Check whether the requested field exists in input structure.
  ctag = strupcase(strtrim(tag, 2))   ;canoncial form of tag
  itag = where(tags eq ctag, nmatch)
  if nmatch eq 0 then return
  itag = itag[0]        ;convert to scalar
  ;print,itag
  ;Copy any fields that precede target field.
  if itag gt 0 then begin     ;target field occurs first
    ;new = create_struct(tags[0], struct.(0))  ;initialize structure
    valt=struct.(0);initialize structure
    nd=SIZE(valt, /N_DIMENSIONS)
    if nd eq 1 then begin
      if n_elements(valt) eq 5 or n_elements(valt) eq 10 then begin
        valf=valt
        siz_str=1
      endif else begin
        valf=valt[0]
        siz_str=n_elements(valt)
      endelse
    endif else begin
      if nd eq 2 then begin
        if n_elements(valt) eq 100 then begin
          valf=valt
          siz_str=1
        endif else begin
          valf=valt[*,0]
          siz_str=n_elements(valt[0,*])
        endelse
      endif else begin
        valf=valt[*,0]
        siz_str=n_elements(valt[0,*])
      endelse
    endelse
    if tags[0] eq 'NEXP' then begin
      valf=0
    endif
    if tags[0] eq 'TARGET_INDEX' then begin
      valf=0
    endif
    new=replicate(create_struct(tags[0],valf),siz_str)
    new.(0)=valt
    for i=1, itag-1 do begin      ;insert leading unchange
      valt=struct.(i)
      nd=SIZE(valt, /N_DIMENSIONS)
      if nd eq 1 then begin
        if n_elements(valt) eq 5 or n_elements(valt) eq 10 then begin
          valf=valt
          siz_str=1
        endif else begin
          valf=valt[0]
          siz_str=n_elements(valt)
        endelse
      endif else begin
        if nd eq 2 then begin
          if n_elements(valt) eq 100 then begin
            valf=valt
            siz_str=1
          endif else begin
            valf=valt[*,0]
            siz_str=n_elements(valt[0,*])
          endelse
        endif else begin
          valf=valt[*,0]
          siz_str=n_elements(valt[0,*])
        endelse
      endelse
      if tags[i] eq 'NEXP' then begin
        valf=0
      endif
      if tags[i] eq 'TARGET_INDEX' then begin
        valf=0
      endif
      values_t=replicate(create_struct(tags[i],valf),siz_str)
      new=struct_addtags(new,values_t)
      new.(i)=valt
    endfor
  endif

  ;Replicate remainder of structure after desired tag.
  for i=itag+1, ntags-1 do begin
    valt=struct.(i)
    nd=SIZE(valt, /N_DIMENSIONS)
    if nd eq 1 then begin
      if n_elements(valt) eq 5 or n_elements(valt) eq 10 or n_elements(valt) eq 35 or n_elements(valt) eq 4 or n_elements(valt) eq 3 or n_elements(valt) eq 2  then begin; the magnitude vector has 5 elements
        valf=valt
        siz_str=1
      endif else begin
        valf=valt[0]
        siz_str=n_elements(valt)
      endelse
    endif else begin
      if nd eq 2 then begin
        if n_elements(valt) eq 100 then begin
          valf=valt
          siz_str=1
        endif else begin
          valf=valt[*,0]
          siz_str=n_elements(valt[0,*])
        endelse
      endif else begin
        if nd eq 3 then begin
          if n_elements(valt[*,*,0]) eq 100 then begin
            valf=valt[*,*,0]
            siz_str=n_elements(valt[0,0,*])
          endif ;else begin
          ;  valf=valt[*,0]
          ;  siz_str=n_elements(valt[0,*])
          ;endelse
        endif else begin
          valf=valt[*,0]
          siz_str=n_elements(valt[0,*])
        endelse
        ;print,valt[*,*,0]
        ;valf=valt[*,0]
        ;siz_str=n_elements(valt[0,*])
      endelse
    endelse
    if tags[i] eq 'NEXP' then begin
      valf=0
    endif
    if tags[i] eq 'TARGET_INDEX' then begin
      valf=0
    endif    
    ;print,tags[i]
    ;print,valt
    ;print,'--'
    ;print,valf
    ;print,nd
    ;print,n_elements(valt)
    ;print,siz_str,n_elements(new)
    values_t=replicate(create_struct(tags[i],valf),siz_str)
    new=struct_addtags(new,values_t)
    new.(i-1)=valt
  endfor

  ;Replace input structure with new structure.
  struct = new

end

pro snr_spec, outdir=outdir, run1d=run1d1, doplot=doplot, spectradir=spectradir, $
  run2d=run2d, plates=plates, legacy=legacy, sky=sky, lite=lite

CPU, TPOOL_NTHREADS = 1  

  spectro_redux = getenv('BOSS_SPECTRO_REDUX')
  if (NOT keyword_set(platefile)) then begin
    platefile = findfile('spField*.fits*', count=nplate)
  endif else begin
    if (size(platefile,/tname) NE 'STRING') then $
      message, 'FieldFILE must be a file name'
    if (keyword_set(platefile)) then nplate = n_elements(platefile) $
    else nplate = 0
  endelse
  if (keyword_set(run1d1)) then run1d = strtrim(run1d1,2) $
  else run1d = getenv('RUN1D')
  
  if (keyword_set(run2d)) then run2d = strtrim(run2d,2) $
  else run2d = getenv('RUN2D')
  comb=spectro_redux+'/'+run2d+'/'
  if (NOT keyword_set(spectradir)) then spectradir=comb
  rsky=1
  if (keyword_set(sky)) then begin 
    rsky = sky
  endif
  if (NOT keyword_set(outdir)) then outdir="outdir"

  logfile = djs_filepath('SNR-' + run2d + '.log', root_dir=outdir)
  if (keyword_set(logfile)) then begin
    cpbackup, logfile
    splog, filename=logfile
    splog, 'Log file ' + logfile + ' opened ' + systime()
  endif


  path=spectro_redux+'/'+run2d+'/conflist.fits'
  data=mrdfits(path,1)
  plates=data.plate
  mjds=data.mjd
  expN=data.nexp_b1
  NexpF=ulong64(total(expN))*500
  finalvalues=replicate(create_struct('SNR_r',0.0),NexpF)
  values_t=replicate(create_struct('ID',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('FIRSTCARTON',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('OBJTYPE',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  snvec_r=finalvalues.snr_r
  idvec=finalvalues.id
  cartonvec=finalvalues.firstcarton
  objtypevec=finalvalues.objtype

  ct=ulong64(0)
  spawn,'mkdir -p '+outdir
  outfile=outdir+'/spSNR-'+run2d+'.fits'
  for itt=0, n_elements(plates)-1 do begin
  platefile='spField-'+strtrim(string(plates[itt]),2)+'-'+strtrim(string(mjds[itt]),2)+'.fits'
  
  platemjd = strmid(fileandpath(platefile), 8, 11)
  zallfile = djs_filepath('spZall-' + platemjd + '.fits', root_dir=run1d)
  zbestfile = djs_filepath('spZbest-' + platemjd + '.fits', root_dir=run1d)
  zlinefile = djs_filepath('spZline-' + platemjd + '.fits', root_dir=run1d)
  logfile = djs_filepath('spec-' + platemjd + '.log', root_dir='')
  fieldid=strmid(platemjd,0,5)
  thismjd=strmid(platemjd,6,5)

  rootdir=comb+'/'+fieldid+'p/'

  splog,'Writing the final output files'
  
  
  plugmap = mrdfits(rootdir+platefile,5)
  zall = mrdfits(rootdir+zallfile,1)
  zbest = mrdfits(rootdir+zbestfile,1)
  zline = mrdfits(rootdir+zlinefile,1)
  zmodel = mrdfits(rootdir+zbestfile,2)
  
  if keyword_set(legacy) or keyword_set(plates) then begin
    fieldid=fieldid+'p'
  endif
  
  target_ind=plugmap.target_index
  single_basefile='coadd/'+thismjd+'/spSpec-'+platemjd+'-'
  single_out_basefile='spec-'+platemjd+'-'
  base_id=platemjd+'-'
  if keyword_set(lite) then begin
    lit_p='/lite'
  endif else begin
    lit_p='/full'
  endelse
  
  for itarget=0, n_elements(target_ind)-1 do begin
    plug_target=plugmap[itarget]
    otype=strtrim(plug_target.objtype,2)
    spec=1
    carton=plug_target.FIRSTCARTON
    if spec eq 1 then begin
      if keyword_set(plates) or keyword_set(legacy) then begin
        indx0 = (where((zbest.fiberid EQ itarget+1)))
        zbest_target=zbest[indx0]
      endif else begin
        zbest_target=zbest[itarget]
      endelse
      indx = (where((zall.fiberid EQ itarget+1)))
      tags = tag_names(zall)
      ntags = n_elements(tags)
      zall_targ=zall[indx]
      indx2 = (where((zline.fiberid EQ itarget+1)))
      tags = tag_names(zline)
      ntags = n_elements(tags)
      zline_targ=zline[indx2]
      if keyword_set(plates) or keyword_set(legacy) then begin;
        if keyword_set(legacy) then begin
           single_file=single_basefile+string(plug_target.fiberid,format='(i4.4)')+'.fits'
           id_str=base_id+string(plug_target.fiberid,format='(i4.4)')
        endif else begin
           if plug_target.catalogid eq 0 then begin
                single_file=single_basefile+string(plug_target.fiberid,format='(i11.11)')+'.fits'
                id_str=base_id+string(plug_target.fiberid,format='(i11.11)')
           endif else begin
                if plug_target.program.contains('offset', /FOLD_CASE ) then begin
                    single_file=single_basefile+strtrim(string(plug_target.catalogid),1)+'.fits'
                    id_str=base_id+strtrim(string(plug_target.catalogid),1)
                endif else begin
                    single_file=single_basefile+string(plug_target.catalogid,format='(i11.11)')+'.fits'
                    id_str=base_id+string(plug_target.catalogid,format='(i11.11)')
                endelse
           endelse
        endelse
      endif else begin
        single_file=single_basefile+plug_target.targetid+'.fits'
      endelse
      single_fileF=rootdir+single_file
      struct_delete_field,plug_target,'objtype'
      struct_delete_field,zbest_target,'fiberid'
      struct_delete_field,zall_targ,'fiberid'
      struct_delete_field,zline_targ,'fiberid'
      if keyword_set(plates) or keyword_set(legacy) then begin
         struct_delete_field,zbest_target,'field'
      endif
      fin_plug=struct_addtags(plug_target,zbest_target)
      nexp=plug_target.nexp
      if keyword_set(plates) or keyword_set(legacy) then begin
         if keyword_set(legacy) then begin
            file_name=single_out_basefile+string(plug_target.fiberid,format='(i4.4)')+'.fits'
         endif else begin
            if plug_target.catalogid eq 0 then begin
                file_name=single_out_basefile+string(plug_target.fiberid,format='(i11.11)')+'.fits'
            endif else begin
                if plug_target.program.contains('offset', /FOLD_CASE ) then begin
                    file_name=single_out_basefile+strtrim(string(plug_target.catalogid),1)+'.fits'
                endif else begin
                    file_name=single_out_basefile+string(plug_target.catalogid,format='(i11.11)')+'.fits'
                endelse
            endelse
         endelse
      endif
      if not keyword_set(lite) then begin
         for i=0, nexp-1 do begin
           single = mrdfits(single_fileF,3+i,hdri)
           flux=single.flux
           ivar=single.ivar
           loglam=single.loglam
           snimg = flux * sqrt(ivar)
           ;gwave = where(loglam GT alog10(4000) AND loglam LT alog10(5500))
           rwave = where(loglam GT alog10(5600) AND loglam LT alog10(6900))
           ;iwave = where(loglam GT alog10(6910) AND loglam LT alog10(8500))
           filtsz=25
           ig = where(ivar[rwave] GT 0, nwave)
           if (nwave GT filtsz) then $
             sntemp = djs_median(snimg[rwave[ig]], $
             width=filtsz, boundary='reflect')
           snr = djs_mean(sntemp)
           snvec_r[ct]=snr
           idvec[ct]=id_str
           cartonvec[ct]=carton
           objtypevec[ct]=otype
           ct=ct+1
         endfor
      endif
      splog,'File '+file_name+' was created, target '+strtrim(string(itarget+1),2)+' of '+strtrim(string(n_elements(target_ind)),2)+' targets'
    endif
      ;stop
  endfor
  endfor
  finalvalues.snr_r=snvec_r
  finalvalues.id=idvec
  finalvalues.firstcarton=cartonvec
  finalvalues.objtype=objtypevec
  ;sxaddpar, coadd_val, 'EXTNAME', 'COADD', ' Coadded spectrum'
  mwrfits, finalvalues, outfile, /create, /silent
  

  splog, 'Successful completion of REFORMAT_SPEC at ' + systime()
  if (keyword_set(logfile)) then splog, /close
 end
