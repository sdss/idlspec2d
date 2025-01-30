import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')
import os.path as ptt
def compare_one(plan, arcfile, out):
    hard = ptt.join(plan.outdir,arcfile.replace('.fit.gz',''))

    fig,ax = plt.subplots(1,2,squeeze=True)
    ax[0].set_title(f'difference between derived trace and {plan.flatfile}')
    for a,b in zip(out,tracetab) :
        ax[0].plot(a-b)
    ax[0].set_ylim(-5,5)
    ax[0].set_xlim(0,len(a))

    # if we have a flat at the field, compare it to derived trace
    field = arc['fieldid']
    match = plan.flats[plan.flats['fieldid'] == field]
    if len(match) > 0 :
        flat=match[-1]['name'][channel%2].astype(str)
        try:
            rff = ptt.join(os.environ['BOSS_SPECTRO_REDUX'],
                           plan.vers, f2s(field),
                           (flat.replace('sdR','spFlat')
                                .replace('.gz','')
                                .replace('.fit','.fitz.gz')))
            with fits.open(rff) as flat_data:
                ax[1].set_title('difference between derived trace and {:s}'.format(flat))
                x=np.arange(len(out[0]))
                ax[1].set_xlim(0,len(x))
                for a,b,c in zip(out,flat_data[5].data,flat_data[0].data) :
                    p=ax[1].plot(x,a-b,ls=':')
                    gd =np.where(c>0)[0]
                    ax[1].plot(x[gd],a[gd]-b[gd],color=p[-1].get_color())
                ax[1].set_ylim(-5,5)
                dt=(flat_data[0].header['TAI-BEG']-im.header['TAI-BEG'])/3600.
                ax[1].text(0.5,0.9,'dt: {:.2f}'.format(dt),transform=ax[1].transAxes)
        except : pass
    fig.tight_layout()
    fig.savefig(hard.replace('.png','_compare.png'))
    plt.close()

def compare(mjd, obs, planfile = None, outfile = None, vers=None):
    
    plan = Plan(mjd, obs, planfile=planfile, outdir=outdir, vers= vers)
    plan.read()

    for channel,cam in zip(plan.channels,plan.cams) :

        for arc in plan.arcs:
            arcfile = arc['name'][channel%2].astype(str)

            outfile = ptt.join(plan.outfile,(arcfile.replace('sdR','spTraceTab')
                                                    .replace('.fit','.fits')))
            out = fits.getdata(outfile,0)
            plot.compare_one(plan, arcfile, out)
