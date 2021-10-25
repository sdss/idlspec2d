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

pro sky_spec, outdir=outdir, run1d=run1d1, doplot=doplot, spectradir=spectradir, $
  run2d=run2d, plates=plates, legacy=legacy, sky=sky, lite=lite

CPU, TPOOL_NTHREADS = 1  

  spectro_redux = getenv('BOSS_SPECTRO_REDUX')
  ;if (NOT keyword_set(platefile)) then begin
  ;  platefile = findfile('spField*.fits*', count=nplate)
  ;endif else begin
  ;  if (size(platefile,/tname) NE 'STRING') then $
  ;    message, 'FieldFILE must be a file name'
   ; if (keyword_set(platefile)) then nplate = n_elements(platefile) $
  ;  else nplate = 0
  ;endelse
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
  spawn,'mkdir -p '+outdir
  logfile = djs_filepath('SKYSpec-' + run2d + '.log', root_dir=outdir)
  if (keyword_set(logfile)) then begin
    cpbackup, logfile
    splog, filename=logfile
    splog, 'Log file ' + logfile + ' opened ' + systime()
  endif
  
  path=spectro_redux+'/'+run2d+'/conflist.fits'
  data=mrdfits(path,1)
  plates=data.plate
  mjds=data.mjd
  program=data.programname
  nt=where(program.contains('RM'))
  plates=plates[nt]
  mjds=mjds[nt]
  program=program[nt]
  NexpF=ulong64(n_elements(mjds))*500
  finalvalues=replicate(create_struct('FLUX_g',0.0),NexpF)
  values_t=replicate(create_struct('FLUX_r',0.0),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('FLUX_i',0.0),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('SKYFLUX_g',0.0),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('SKYFLUX_r',0.0),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('SKYFLUX_i',0.0),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('ID',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('FIRSTCARTON',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('OBJTYPE',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  values_t=replicate(create_struct('PROGRAM',""),NexpF)
  finalvalues=struct_addtags(finalvalues,values_t)
  snvec_g=finalvalues.flux_g
  snvec_r=finalvalues.flux_r
  snvec_i=finalvalues.flux_i
  skyvec_g=finalvalues.skyflux_g
  skyvec_r=finalvalues.skyflux_r
  skyvec_i=finalvalues.skyflux_i
  idvec=finalvalues.id
  cartonvec=finalvalues.firstcarton
  objtypevec=finalvalues.objtype
  programvec=finalvalues.program

  ct=ulong64(0)
  spawn,'mkdir -p '+outdir
  outfile=outdir+'/spSKYspec-'+run2d+'.fits'
  for itt=0, n_elements(plates)-1 do begin
  programF=program[itt]
  platefile='spField-'+strtrim(string(plates[itt]),2)+'-'+strtrim(string(mjds[itt]),2)+'.fits'
  
  platemjd = strmid(fileandpath(platefile), 8, 11)
  fieldid=strmid(platemjd,0,5)
  thismjd=strmid(platemjd,6,5)

  rootdir=comb+'/'+fieldid+'p/'

  splog,'Analysing the spec files'
  
  if FILE_TEST(rootdir+platefile) then begin
  plugmap = mrdfits(rootdir+platefile,5)
  
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
         ;for i=0, nexp-1 do begin
           single = mrdfits(single_fileF,1,hdri)
           flux=single.flux
           ivar=single.ivar
           loglam=single.loglam
           skyfkux=single.sky
           snimg = flux ;* sqrt(ivar)
           
           gwave = where(loglam GT alog10(4000) AND loglam LT alog10(5500))
           rwave = where(loglam GT alog10(5600) AND loglam LT alog10(6900))
           iwave = where(loglam GT alog10(6910) AND loglam LT alog10(8500))
           
           filtsz=25
           ig = where(ivar[gwave] GT 0 AND flux[gwave] GE 0, nwave)
           sntemp=0
           sntemp2=0
           if (nwave GT filtsz) then begin
             sntemp = djs_median(snimg[gwave[ig]], width=filtsz, boundary='reflect')
             sntemp2 = djs_median(skyfkux[gwave[ig]], width=filtsz, boundary='reflect')
           endif
           snr = djs_mean(sntemp)
           skyf= djs_mean(sntemp2)
           snvec_g[ct]=snr
           skyvec_g[ct]=skyf
           
           filtsz=25
           ig = where(ivar[rwave] GT 0 AND flux[rwave] GE 0, nwave)
           sntemp=0
           sntemp2=0
           if (nwave GT filtsz) then begin
             sntemp = djs_median(snimg[rwave[ig]], width=filtsz, boundary='reflect')
             sntemp2 = djs_median(skyfkux[rwave[ig]], width=filtsz, boundary='reflect')
           endif
           snr = djs_mean(sntemp)
           skyf= djs_mean(sntemp2)
           snvec_r[ct]=snr
           skyvec_r[ct]=skyf
           
           filtsz=25
           ig = where(ivar[iwave] GT 0 AND flux[iwave] GE 0, nwave)
           sntemp=0
           sntemp2=0
           if (nwave GT filtsz) then begin
             sntemp = djs_median(snimg[iwave[ig]], width=filtsz, boundary='reflect')
             sntemp2 = djs_median(skyfkux[iwave[ig]], width=filtsz, boundary='reflect')
           endif
           snr = djs_mean(sntemp)
           skyf= djs_mean(sntemp2)
           snvec_i[ct]=snr
           skyvec_i[ct]=skyf
           
           
           idvec[ct]=id_str
           cartonvec[ct]=carton
           objtypevec[ct]=otype
           programvec[ct]=programF
           ct=ct+1
         ;endfor
      endif
      splog,'Spectrum '+file_name+' analysed, target '+strtrim(string(itarget+1),2)+' of '+strtrim(string(n_elements(target_ind)),2)+' targets'
    endif
      ;stop
  endfor
  endif
  endfor
  
  finalvalues.flux_g=snvec_g
  finalvalues.flux_r=snvec_r
  finalvalues.flux_i=snvec_i
  finalvalues.skyflux_g=skyvec_g
  finalvalues.skyflux_r=skyvec_r
  finalvalues.skyflux_i=skyvec_i
  finalvalues.id=idvec
  finalvalues.firstcarton=cartonvec
  finalvalues.objtype=objtypevec
  finalvalues.program=programvec
  
  mwrfits, finalvalues, outfile, /create, /silent
  

  splog, 'Successful completion of REFORMAT_SPEC at ' + systime()
  if (keyword_set(logfile)) then splog, /close
 end
