;+
; NAME:
;   doscatter
;
; PURPOSE:
;   Generate a scattering surface for a given image.
;
;
; CALLING SEQUENCE:
;   surf = doscatter(camera, image, ivar, wset)
;
; INPUTS:
;   camera    - camera name, e.g. 'r1'
;   image     - image to generate the scattering surface for. 
;   ivar      - inverse variance of image.
;   wset      - wavelength solution 
;
; OPTIONAL INPUTS:
;   mjd       - of the observation (default=99999)
;   bandrows  - how many rows to process at once, with a
;               single set of parameters. (default=5)
;   params    - if set, the parameters to use.
;   configfile- if set, the filename to get all the parameters from. 
;   configdir - if set, the directory to find configfile in.
;
; OUTPUT:
;   the scattered light surface.
;
; INTERNAL SUPPORT ROUTINES:
;

function doscatter, camera, image, ivar, wset, mjd=mjd, $
                    bandrows=bandrows, params=params, sigma=sigma, $
                    configfile=configfile, configdir=configdir

    if not keyword_set(sigma) then sigma = 0.9

    color = strmid(camera,0,1)
    if not keyword_set(mjd) then mjd=99999

                                ; This is an important but imperfect
                                ; step: the image flux drives the
                                ; scattering surface, so saturated
                                ; pixels and bad columns need to have
                                ; some reasonable values. Because most
                                ; the the masked pixels are in bad
                                ; columns, interpolate from the
                                ; sides. We could do much better, since we
                                ; know where the traces are and
                                ; basically what their shapes are.
    fitimg = djs_maskinterp(image, ivar EQ 0, iaxis=0)
    imgsize = size(fitimg,/dim)
    nr = imgsize[1]
    nc = imgsize[0]

    traceset2xy, wset, wx, wy

                                ; Fetch the correct parameters
    if (NOT keyword_set(configfile)) then begin
       if not keyword_set(configdir) then $
          configdir = filepath('', $
                               root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
       splog, 'Looking for scattering file in ', configdir
       configfile = findopfile(string('scattering-', camera, '-*.fits'), $
                               mjd, configdir, $
                               /abort_notfound, silent=silent)
       configfile = configdir + '/' + configfile
    end

    if not keyword_set(params) then begin
       params = mrdfits(configfile, 0, hdr1)
                                ; Prepend zero DC term and request that a
                                ; core be notched out.
       if size(params, /n_dim) eq 1 then begin
          params = [0, sigma, params]
       end else begin
          sz = size(params,/dim)
          nparams = dblarr(sz[0]+2, sz[1])
          nparams[1,*] = params[0,*]
          nparams[2:sz[0]+1,*] = params
          params = nparams
       end
    end

    if color EQ 'r' then begin
       evalfunc = 'redhalo'
       maxkernr = 256
       addrows = 64
       bandrows = keyword_set(bandrows) ? bandrows : 8
    end else if color eq 'b' then begin
       evalfunc = 'bluescat'
       maxkernr = 64
       addrows = 32
       bandrows = keyword_set(bandrows) ? bandrows : 4
    end else message, "unknown color, camera: ", color, camera


                                ; Pass the evaluator a single
                                ; wavelength per row, in microns. 
    if size(wy, /n_dim) GT 1 then $
       wy = djs_median(wy,2)

                                ; Scattered light comes in from beyond
                                ; the detector region of the chip. We
                                ; only _really_ care about the red end
                                ; of red, and slightly care about the
                                ; blue end of blue. Fabricate a lie
                                ; for the evaluator by extending the
                                ; wavelength solution and copying the
                                ; last row of the detector. Actually,
                                ; we want to extend on the sides as
                                ; well, and probably attenuate the
                                ; extended flux.
    if keyword_set(addrows) then begin
       splog, "adding rows: ", addrows

       dlam = wy[nr-1] - wy[nr-2]
       elam = findgen(addrows)/(addrows-1) * (dlam * addrows) + (wy[nr-1] + dlam)
       wy = [wy, elam]

       fitimg = [[fitimg], [rebin(fitimg[*,nr-1], [nc, addrows])]]

       if size(params, /n_dim) GT 1 then begin
          nparams = (size(params, /dim))[0]
          params = [[params], [rebin(params[*,nr-1], [nparams, addrows])]]
       end
    end

    fitmicrons = 10^wy / 10000
    ymodel = gen_scatter(fitimg, fitmicrons, params, $
                         evalfunc=evalfunc, bandrows=bandrows, $
                         maxkernr=maxkernr, keeprows=[0,nr-1])
    return, ymodel
 end
