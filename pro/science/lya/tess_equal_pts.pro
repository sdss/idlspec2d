;+
; NAME:
;   tess_equal_pts
;
; PURPOSE:
;   Tesselate a point distribution of points on the sky into chunks
;   with an equal number per chunk
;
; CALLING SEQUENCE:
;   tess_equal_pts, infile, [ outfile, nchunk=, ntry=, niter=, $
;    srange=, /debug, plotfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   infile     - Input FITS binary table file containing RA,DEC
;   outfile    - Output FITS binary table with identical format at INFILE
;                but with the addition of an ICHUNK tag running from 0 to
;                NCHUNK-1;
;                default to INFILE replacing '.fits' with '-ichunk.fits'
;   nchunk     - Number of tessellation chunks; default to 20
;   ntry       - Number of offsets to test in each step in the annealing
;                of the solution; default to 20
;   niter      - Number of iterations of annealing; default to 200
;   srange     - Range of mean offsets for testing tessellation center offsets;
;                default to [0.01,3] deg; begin the iterations with offsets
;                with an RMS distribution of SRANGE[1]; each iteration
;                scales this down by 10% (if not improvement found) or up
;                by 10% (if improvement is found) but always bounded by
;                the SRANGE[0],SRANGE[1]
;   debug      - If set, then make sky plot at each iteration
;   plotfile   - If set, then make final PostScript plot; either string
;                value for the file name, or set /PLOTFILE for default name
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;   30-May-2013  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
; The annealing code improves the following metric in each step which
; has three terms:
;   1. Standard dev. of the number of objects assigned to each chunk
;   2. Maximum number of objects in any chunk
;   3. Total of all object-to-chunk-center distances
; The first half of the iterations include the 3rd term, which is then
; tuned to zero for the second half of iterations
function tess_equal_val, cnum, dmin, nchunk=nchunk, detune=detune
   vhist = lonarr(nchunk)
   populate_image, vhist, cnum
   term1 = stddev(vhist)
   term2 = max(vhist)
   term3 =  total(dmin)*detune
   return, term1 + term2 + term3
end

;------------------------------------------------------------------------------
pro tess_equal_plot, objs, cobj, cnum, nchunk=nchunk, isort=isort, qngc=qngc, $
 psfile=psfile

   if (qngc) then begin
      aitoff, 180-objs.ra, objs.dec, xplot, yplot
      aitoff, 180-cobj.ra, cobj.dec, xplotc, yplotc
      title = 'NGC Aitoff projection'
   endif else begin
      aitoff, 270-objs.ra, objs.dec, xplot, yplot
      aitoff, 270-cobj.ra, cobj.dec, xplotc, yplotc
      title = 'SGC Aitoff projection'
   endelse

   ; Plot each chunk, sorted by RA of the chunk center
   colorvec = ['red','green','blue','cyan','magenta','grey', $
    'yellow','light red','light green','light blue', $
    'light cyan','light magenta']
   if (keyword_set(psfile)) then $
    djs_plot, xplot, yplot, /nodata, xstyle=5, ystyle=5, title=title $
   else $
    splot, xplot, yplot, /nodata, xstyle=5, ystyle=5, title=title
   for j=0L, nchunk-1L do begin
      indx = where(cnum EQ isort[j], ct)
      if (ct GT 0) then $
       if (keyword_set(psfile)) then $
        djs_oplot, xplot[indx], yplot[indx], psym=3, $
         color=djs_icolor(colorvec[j MOD n_elements(colorvec)]) $
       else $
        soplot, xplot[indx], yplot[indx], psym=3, $
         color=djs_icolor(colorvec[j MOD n_elements(colorvec)])
   endfor
   if (keyword_set(psfile)) then $
    djs_oplot, xplotc, yplotc, psym=2 $
   else $
    soplot, xplotc, yplotc, psym=2

   return
end

;------------------------------------------------------------------------------
pro tess_equal_pts, infile, outfile, nchunk=nchunk1, ntry=ntry1, niter=niter1, $
 srange=srange1, debug=debug, plotfile=plotfile1

   if (NOT keyword_set(infile)) then $
    message, 'Must specify INFILE'
   if (NOT keyword_set(outfile)) then begin
      thisfile = fileandpath(infile, path=thispath)
      ipos = strpos(thisfile, '.fits', /reverse_search)
      if (ipos EQ -1) then ipos = strlen(thisfile)
      outfile = djs_filepath(strmid(thisfile,0,ipos)+'-ichunk.fits', $
       root_dir=thispath)
   endif
   if (keyword_set(plotfile1)) then begin
      if (size(plotfile1,/tname) EQ 'STRING') then begin
         plotfile = plotfile1
      endif else begin
         thisfile = fileandpath(outfile, path=thispath)
         ipos = strpos(thisfile, '.fits', /reverse_search)
         if (ipos EQ -1) then ipos = strlen(thisfile)
         plotfile = djs_filepath(strmid(thisfile,0,ipos)+'.ps', $
          root_dir=thispath)
      endelse
   endif

   if (keyword_set(nchunk1)) then nchunk = nchunk1 $
    else nchunk = 20
   if (keyword_set(ntry1)) then ntry = ntry1 $
    else ntry = 20
   if (keyword_set(niter1)) then niter = niter1 $
    else niter = 200
   if (n_elements(srange1) EQ 2) then srange = srange1 $
    else srange = [0.01, 3.]

   objs = mrdfits(infile,1)
   nobj = n_elements(objs) * keyword_set(objs)
   if (nobj LT nchunk) then $
    message, 'Fewer objects than chunks!'
   tags = tag_names(objs)
   if (total(tags EQ 'RA') + total(tags EQ 'DEC') NE 2) then $
    message, 'Input file must contain RA,DEC!'

   if (keyword_set(debug) OR keyword_set(plotfile)) then begin
      euler, objs.ra, objs.dec, ll, bb, 1
      qngc = total(bb GT 0) GT 0
      ll = 0
      bb = 0
   endif

   ; Start with a list of chunks centered on a random subset
   ; of the object positions
   ; Start with a list of chunks centered on a random subset
   ; of the object positions
   indx = (sort(randomu(1234,nobj)))[0:nchunk-1]
   cobj = replicate(create_struct('ra', 0d, 'dec', 0d), nchunk)
   cobj.ra = objs[indx].ra
   cobj.dec = objs[indx].dec

   ; Compute the initial chunks and the metric
   darr = dblarr(nchunk, nobj)
   for ichunk=0L, nchunk-1L do $
    darr[ichunk,*] = djs_diff_angle(objs.ra, objs.dec, $
     cobj[ichunk].ra, cobj[ichunk].dec)
   cnum = lonarr(nobj) ; chunk number for each object (closest center)
   dmin = dblarr(nobj) ; distance to this
   for i=0L, nobj-1L do begin
      dmin[i] = min(darr[*,i], j)
      cnum[i] = j
   endfor
   detune0 = 1000. * float(nchunk) / float(nobj)
   detune = detune0
   tvalue = tess_equal_val(cnum, dmin, nchunk=nchunk, detune=detune)

   ; Make a list of possible offsets
   ;   inum - Chunk center number to change
   ;   ra, dec - New location for that center
   ;   dthis[NOBJ] - Distance to that center
   ;   dnew[NOBJ] - Best new distance
   ;   cnew[NOBJ] - Best new chunk number
   ;   vnew - New value
   tobj = replicate(create_struct('inum', 0L,'ra', 0d, 'dec', 0d, $
    'dthis', dblarr(nobj), 'dnew', dblarr(nobj), 'cnew', lonarr(nobj), $
    'vnew', 0d), ntry)

   ; Iterate
   scale = srange[1]
   for iiter=0L, niter-1L do begin
      tobj.inum = long(randomu(2345+iiter, ntry) * nchunk) < (nchunk-1)
      for itry=0L, ntry-1L do begin
         tobj[itry].ra = cobj[tobj[itry].inum].ra + randomn(3456+itry) * scale
         tobj[itry].dec = cobj[tobj[itry].inum].dec + randomn(4567+itry) * scale
         tobj[itry].dthis = djs_diff_angle(tobj[itry].ra, tobj[itry].dec, $
          objs.ra, objs.dec)

         ; The objects that could change are those where the chunk number
         ; is already this chunk, or where the distance to this chunk is
         ; smaller than the existing best distance
         tobj[itry].cnew = cnum ; start with the existing chunk assignments
         tobj[itry].dnew = dmin ; start with the existing chunk assignments
         indx = where(cnum EQ tobj[itry].inum $
          OR tobj[itry].dthis LT dmin, nn)
         for i=0L, nn-1L do begin
            dtmp = darr[*,indx[i]]
            dtmp[tobj[itry].inum] = tobj[itry].dthis[indx[i]]
            tobj[itry].dnew[indx[i]] = min(dtmp, itmp)
            tobj[itry].cnew[indx[i]] = itmp
         endfor

         detune = detune0*float(0.5*niter-iiter)/niter > 0
         tobj[itry].vnew = tess_equal_val(tobj[itry].cnew, $
          tobj[itry].dnew, nchunk=nchunk, detune=detune)
      endfor

      ; See if the best of these possible offsets improves the solution,
      ; and if so then move to it
      vbest = min(tobj.vnew, ibest)
      if (vbest LT tvalue) then begin
         cobj[tobj[ibest].inum].ra = tobj[ibest].ra
         cobj[tobj[ibest].inum].dec = tobj[ibest].dec
         darr[tobj[ibest].inum,*] = tobj[ibest].dthis
         dmin = tobj[ibest].dnew
         cnum = tobj[ibest].cnew
         tvalue = vbest
         scale = (scale*1.10) < srange[1]
      endif else begin
         scale = (scale/1.10) > srange[0]
      endelse
      print, 'Iteration ', iiter, ' of ', niter, ' metric=', tvalue, $
       ' scale=', scale

      if (keyword_set(debug)) then $
       tess_equal_plot, objs, cobj, cnum, nchunk=nchunk, isort=sort(cobj.ra), $
        qngc=qngc
   endfor

   if (keyword_set(plotfile)) then begin
      print, 'Writing plot file '+plotfile
      dfpsplot, plotfile, /square, /encap, /color
      tess_equal_plot, objs, cobj, cnum, nchunk=nchunk, isort=sort(cobj.ra), $
       qngc=qngc, /psfile
      dfpsclose
   endif

   for ichunk=0L, nchunk-1L do $
    print, 'Chunk #', ichunk,' contains ', $
     long(total(cnum EQ ichunk)), ' objects'

   ; Write output file
   print, 'Writing output file '+outfile
   outdat = struct_addtags(objs, replicate(create_struct('ichunk',0L),nobj))
;   outdat = replicate(create_struct('ra', 0d, 'dec', 0d, 'ichunk', 0L), nobj)
;   outdat.ra = objs.ra
;   outdat.dec = objs.dec
   outdat.ichunk = cnum
   mwrfits, outdat, outfile, /create

   return
end
;------------------------------------------------------------------------------
