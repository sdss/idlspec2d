;+
; NAME:
;       read_gaia
; PURPOSE:
;       Read GAIA stars near a given RA, dec with the implementation of the DR2 and DR1
;
; INPUTS:
;   racen     - RA of region center (J2000)    [degrees]
;   deccen    - dec of region center (J2000)   [degrees]
;   ra        - RA of targets
;   dec       - DEC of targets
;
; OPTIONAL INPUTS:
;   gaia_dr2        - data release name (default 'dr2')
;
; OUTPUTS:
;   ra        - new Gaia RA of targets
;   dec       - new Gaia DEC of targets
;   pmra  -proper motions from Gaia
;   pmdec -proper motions from Gaia
;   mjdepoch -MJD of the Gaia DR epoch
;
; COMMENTS:
;   Assumes data in $GAIA_DATA/[dr]/gaia_source/fits_sorted
;   produced by gaia_setup.py.
;
; REVISION HISTORY
;   2016-Oct-20   MRB, NYU
;
;------------------------------------------------------------------------
;
pro rm_read_gaia, ra0,dec0,stars,gaia_dr2=gaia_dr2,dist_std=dist_std;,mjdepoch=mjdepoch,pmra=pmra,pmdec=pmdec,ra=ra,dec=dec,stars=stars

  if not keyword_set(gaia_dr2) then begin
    gaia_dr2=1
  endif
  if keyword_set(gaia_dr2) then begin
    gaiadir=getenv('GAIA_DATA')+'/dr2/gaia_source/csv/'
    gaia_list=getenv('IDLSPEC2D_DIR')+'/opfiles/list_gaia'
    readcol,gaia_list,format='A',gaiafile_l
    dpf_nest2ring, 4096, lindgen(12*long(4090)*long(4096)), ipring
    gaia_id1=lon64arr(n_elements(gaiafile_l))
    gaia_id2=lon64arr(n_elements(gaiafile_l))
    for it=0l, n_elements(gaiafile_l)-1 do begin
      gaia_id_s=strsplit(strjoin(strsplit(gaiafile_l[it],".",/extract),"_"),"_",/extract)
      gaia_id1[it]=long64(gaia_id_s[1]);/34359738368
      gaia_id2[it]=long64(gaia_id_s[2]);/34359738368
    endfor
  endif else begin
    gaiadir=getenv('GAIA_DATA')+'/dr1/gaia_source/fits_sorted/'
  endelse
  ra_std=stars.ra
  dec_std=stars.dec
  dist_std=stars.dec*0.0
  ;epoch_std=stars.epoch
  ;rapm_std=stars.pmra
  ;decpm_std=stars.pmdec
  ;to=0
  for i=0, 6 do begin
    for j=0, 6 do begin
      if keyword_set(gaia_dr2) then begin
        ra0t=(ra0-3+i)
        dec0t=dec0-3+j
        phi=ra0t*!pi/180.0
        the=(90-dec0t)*!pi/180.0
        ang2pix_ring,4096,the,phi,pixi_r
        idt=where(ipring eq pixi_r[0],indt)*34359738368
        if indt gt 0 then begin
          ga_id=idt[0]
          nt=where((gaia_id1 le ga_id) and (gaia_id2 ge ga_id), nga)
          readcol,gaiadir+gaiafile_l[nt[0]],delimiter=",",format="L,A,L,I,F,D,D,D,D,D,D,D,D,D,D,D,D", $
            a,b,c,d,epch_g,ra_g,ra_e_g,dec_g,de_g,dec_e_g, $
            pr,pr_e,pr_ve,pmra_g,pmra_g_e,pmdec_g,pmdec_g_e,/compress,/silent
          spherematch, ra_g, dec_g, ra_std, dec_std, 1./3600, m1, m2
          if m2[0] ne -1 then begin
            ra_std[m2]=ra_g[m1];tmpdt0.ra
            dec_std[m2]=dec_g[m1];tmpdt0.dec
            dist_std[m2]=1.0/abs((pr[m1]-0.0)*1e-3);zero point parallax
            ;rapm_std[m2]=pmra_g[m1];tmpdt0.pmra*0.0
            ;decpm_std[m2]=pmdec_g[m1];tmpdt0.pmdec*0.0
            ;epoch_std[m2]=epch_g[m1];tmpdt0.ref_epoch
            ;to=to+n_elements(m1)
          endif
          ;exit
        endif

      endif else begin
        gaiafile=gaiadir+'gaia-dr1-'+strtrim(string(uint(ra0-3+i)),2)+'-'+strtrim(string(uint(dec0-3+j+90)),2)+'.fits'
        tmp_data = mrdfits(gaiafile, 1, hdr2)
        spherematch, tmp_data.ra, tmp_data.dec, ra_std, dec_std, 1./3600, m1, m2
        if m2[0] ne -1 then begin
          tmpdt0=tmp_data(m1)
          ra_std[m2]=tmpdt0.ra
          dec_std[m2]=tmpdt0.dec
          dist_std[m2]=0.0
          ;rapm_std[m2]=0.0;tmpdt0.pmra*0.0
          ;decpm_std[m2]=0.0;tmpdt0.pmdec*0.0
          ;epoch_std[m2]=tmpdt0.ref_epoch
          ;to=to+n_elements(m1)
        endif
      endelse
    endfor
  endfor
  ;stars.ra=ra_std
  ;stars.dec=dec_std
  ;stars.pmra=rapm_std
  ;stars.pmdec=decpm_std
  ;stars.epoch=epoch_std
  ;nt1=where(epoch_std gt 2008.0,np)
  ;if np gt 0 then begin
  ;  ra=ra[nt1]
  ;  dec=dec[nt1]
  ;  pmra=pmra[nt1]
  ;  pmdec=pmdec[nt1]
  ;  epoch=epoch[nt1]
  ;  stars=stars(nt1)
  ;endif
  
  ;mjdepoch=(1461.0*(epoch+4800.0 +(0.0-14.)/12.))/4.+(367*(0.0-2-12.* ((0.0-14.)/12.)))/12.-(3.0*((epoch+4900.+(0.0-14.0)/12.)/100.))/4.-32075.0-2400000.5D
  ;return
end