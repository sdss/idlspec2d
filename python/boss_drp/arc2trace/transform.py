from boss_drp.arc2trace import html, Plan, plot, Report
from boss_drp.field import field_to_string as f2s
from boss_drp.sos.report_err import report
from pyvista import imred, stars, image

import copy
import numpy as np
import os
import os.path as ptt
import multiprocessing as mp

def transform(mjd,obs='lco',thresh=400,nskip=40, rad=2,
                  clobber=False, outdir=None, planfile = None,
                  threads=8, cams=None, sos=False, designMode=None,
                  vers=None, outdir = None) :
    """ Get transformations from first arc for all arcs on a given MJD
        Make plots and HTML page

        Parameters
        ----------
        mjd  : int
               MJD to process
        obs  : str
               observatory = 'lco'|'apo', default='lco'
        nskip : int
               use every nskip line, default=False
        clobber : bool
        outdir : str
        vers : str
               string for BOSS version, with planfile
        planfile : str
               name of plan file with full path
               
    """

    plan = Plan(mjd, obs, planfile=planfile, outdir=outdir, vers= vers)
    plan.read()
    reports = Reports(plan)
    boss=imred.Reducer('BOSS',dir=ptt.join(os.environ[plan.data_env], f'{mjd}'))
    for channel,cam in zip(plan.channels,plan.cams) :

        flatfile = plan.flat2traceflat(channel)
        tracetab = fits.open(ptt.join(plan.outdir,flatfile)[5].data
        col0 = np.arange(len(tracetab[0]))

        reffile = plan.raw2out(plan.tracearc, channel)

        im0=boss.reduce(reffile,channel=channel)
        print('finding lines in reference: ', reffile)
        lines0=stars.find(im0.data,thresh=thresh,
                          sharp=[0,0.5],round=[-0.25,0.75])[::nskip]
        print('automarking lines in reference: ', reffile)
        lines=stars.automark(im0.data,lines0,rad=rad,dx=0,dy=0,
                             background=False,func='marginal_gfit')

        # process all of the arcs in parallel with multiprocessing pool if threads>0
        pars=[]
        for arc in plan.arcs :
            # build up input parameters
            arcfile = plan.raw2out(arc['name'], channel, exten='.fit.gz')
            
            hard = ptt.join(outdir,arcfile.replace('.fit.gz',''))
            if clobber or ptt.exists(hard.replace('.png','_compare.png'))==False :
                print(arcfile,channel)
                im=boss.reduce(arcfile+'.gz',channel=channel)
                pars.append((im0,im,lines,hard))

        # run the solutions for each arc
        if threads > 0 :
            pool = mp.Pool(threads)
            outputs = pool.map_async(transform_thread, pars).get()
            pool.close()
            pool.join()
        else :
            outputs=[]
            for par in pars :
                outputs.append(transform_thread(par))

        # create the outputs
        iout=0
        for arc in plan.arcs:
            arcfile = plan.raw2out(arc['name'], channel, exten='.fit.gz')
            if outputs[iout] is None :
                message = f"Failure in Fitting Arc to Traces with {arcfile}"
                print(message)
                if sos:
                    message = f"BOSS_ARCS_TO_TRACES: {arcfile}: ABORT: Failure in Fitting Arc to Traces (Please retake calibrations for this field)"
                    reports.add(arcfile, cam, message, designMode)
                continue

            linfit = outputs[iout][0]
            im = pars[iout][1]
            iout+=1

            # calculate derived tracetable
            out=copy.copy(tracetab)
            for irow,row0 in enumerate(tracetab) :
                fit = linfit(np.vstack([row0,col0]).T)
                out[irow,:]= np.interp(np.arange(len(row0)),fit[:,1],fit[:,0])

            # write derived trace table
            hdu=fits.PrimaryHDU(out)
            hdu[0].header.set('EXTNAME','XCEN', 'Updated XCEN Trace Positions')
            hdu[0].header.set('TRACFLAT', plan.flat2traceflat(channel), 'spTraceFlat used in transform')
            outfile = plan.arc2tracetab(arc['name'], channel)
            hdu.writeto(ptt.join(plan.outdir,outfile), overwrite=True)
            
            plot.compare_one(plan, arcfile, out)
    html(plan)
    if sos:
        reports.finalize()

def transform_thread(pars) :
    try :
        return image.transform(pars[0],pars[1],pars[2],hard=pars[3])
    except :
        return None
