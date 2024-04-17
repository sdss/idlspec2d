import matplotlib
matplotlib.use('agg')
import matplotlib.patheffects as PathEffects
from matplotlib import pyplot as plt
from numpy import ma
import numpy as np
import os.path as ptt
from astropy.io import fits
from astropy.table import Table, unique

def plot_sky_locations(topdir, flist_file, splog):
    splog.info('Producing Field Location Plots')
    with fits.open(ptt.join(topdir, flist_file), memmap=True) as hdul:
        flist = hdul[1].data
        idx = np.where(np.char.strip(flist['STATUS1D'].data) == 'Done')[0]

        RA   = flist['RACEN'][idx]
        DEC  = flist['DECCEN'][idx]
        prog = np.char.upper(flist['PROGRAMNAME'][idx])
        fcad = np.char.lower(flist['FIELD_CADENCE'][idx])
        fid  = flist['FIELD'][idx]
        fsur = np.char.lower(flist['SURVEY'][idx])
        status = np.char.lower(flist['STATUS1D'][idx])
        nexp = flist['NEXP'][idx]
        flist = None
# plot the RA/DEC in an area-preserving projection
# convert coordinates to degrees
    RA *= np.pi / 180
    DEC *= np.pi / 180
    phi = np.linspace(0, 2.*np.pi, 36)  #36 points
    r = np.radians(1.5)
    C0=C1=C2=C3=C4=C5=C6=C7=0

    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.0), layout='constrained', subplot_kw={'projection': 'mollweide'})
    plt.grid(True)
    plt.title('SDSS plate/field locations')
    for i in range(0, len(RA)):
        
        if RA[i] < np.pi:
            x = RA[i] + r*np.cos(phi)
        else:
            x = RA[i] + r*np.cos(phi)-2*np.pi
        y = DEC[i] + r*np.sin(phi)
        if 'dark' in fcad[i].lower():
            label = 'FPS DARK'  if C0 == 0 else None
            ax.plot(x, y, color = 'blue', label = label, alpha=.5)
            C0 = 1
        elif 'bright' in fcad[i]:
            label = 'FPS BRIGHT'  if C1 == 0 else None
            ax.plot(x, y, color = 'orange', label = label, alpha=.5)
            C1 = 1
        elif 'plate' in fcad[i]:
            label = 'PLATE'  if C2 == 0 else None
            ax.plot(x, y, color = 'green', label = label, alpha=.5)
            C2 = 1
        else:
#            pass
            label = 'MANUAL'  if C3 == 0 else None
            ax.plot(x, y, color = 'saddlebrown', label = label, alpha=.5)
            C3 = 1

    plt.legend(loc=1,fontsize=10)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    txt = ax.set_xticklabels(xlab)
    for t in txt:
        t.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
    plt.savefig(ptt.join(topdir,'SDSSVc_s.png'),dpi=50,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSVc.png'),dpi=500,bbox_inches='tight')
    plt.close()

    C0=C1=C2=C3=C4=C5=C6=C7=C8=C9=0

####################################################################################
    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.0), layout='constrained', subplot_kw={'projection': 'mollweide'})
    plt.grid(True)
    plt.title('SDSS-V plate/field locations')
    for i in range(0, len(RA)):
        if RA[i] < np.pi:
            x = RA[i] + r*np.cos(phi)
        else:
            x = RA[i] + r*np.cos(phi)-2*np.pi
        y = DEC[i] + r*np.sin(phi)
        if 'RM' in prog[i]:
            label = 'RM'  if C0 == 0 else None
            ax.plot(x, y, color = 'blue', label = label)
            C0 = 1
        elif 'MWM' in prog[i]:
            label = 'MWM'  if C1 == 0 else None
            ax.plot(x, y, color = 'red', label = label)
            C1 = 1
        elif 'AQMES-Medium'.upper() in prog[i]:
            label = 'AQMES-Medium'  if C2 == 0 else None
            ax.plot(x, y, color = 'pink', label = label)
            C2 = 1
        elif 'AQMES-Wide'.upper() in prog[i]:
            label = 'AQMES-Wide'  if C3 == 0 else None
            ax.plot(x, y, color = 'orange', label = label)
            C3 = 1
        elif 'AQMES-Bonus'.upper() in prog[i]:
            label = 'AQMES-Bonus'  if C4 == 0 else None
            ax.plot(x, y, color = 'magenta', label = label)
            C4 = 1
        elif 'eFEDS'.upper() in prog[i]:
            label = 'eFEDS' if C5 == 0 else None
            ax.plot(x, y, color = 'green', label = label)
            C5 = 1
        elif 'OFFSET'.upper() in prog[i]:
            label = 'OFFSET'  if C6 == 0 else None
            ax.plot(x, y, color = 'black', label = label)
            C6 = 1
        elif 'mwm-bhm' in fsur[i]:
            label = 'MWM-BHM FPS'  if C7 == 0 else None
            ax.plot(x, y, color = 'navy', label = label)
            C7 = 1
        elif 'bhm-mwm' in fsur[i]:
            label = 'BHM-MWM FPS'  if C8 == 0 else None
            ax.plot(x, y, color = 'lightcoral', label = label)
            C8 = 1
        else:
            ax.plot(x, y, color = 'saddlebrown')
            C9 = 1
    plt.legend(loc=1,fontsize=10)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    txt = ax.set_xticklabels(xlab)
    for t in txt:
        t.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
    plt.savefig(ptt.join(topdir,'SDSSV_s.png'),dpi=50,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSV.png'),dpi=500,bbox_inches='tight')
    plt.close()


####################################################################################
    allpointings = Table()
    if len(idx) == 0:
        allpointings['RA'] = [np.NaN]
        allpointings['DEC'] = [np.NaN]
        allpointings['Nexp'] = [np.NaN]
    else:
        allpointings['RA']  = RA * 180/np.pi
        allpointings['DEC'] = DEC * 180/np.pi
        allpointings['Nexp'] = nexp
        RA   = DEC = prog = fcad = fsur = status = None
        
    pointings = unique(allpointings, keys=['RA','DEC'])
    pointings.add_column(0, name='Nobs')
    pointings['Nexp'] = 0

    pointings.add_column(0.0, name='x')
    pointings.add_column(0.0, name='y')
    for row in pointings:
        idx = np.where((allpointings['RA'].data == row['RA']) &
                       (allpointings['DEC'].data == row['DEC']))[0]
        row['Nobs'] = len(idx)
        row['Nexp'] = np.sum(allpointings[idx]['Nexp'].data)
        if row['RA']*np.pi / 180 < np.pi:
            row['x'] = row['RA']*np.pi / 180
        else:
            row['x'] = row['RA']*np.pi / 180 - 2*np.pi
        row['y'] = row['DEC']*np.pi / 180
                
    
    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.25),layout='constrained', subplot_kw={'projection': 'mollweide'})
    with np.errstate(invalid="ignore"):
        plt.scatter(pointings['x'].data, pointings['y'].data, s=30,
                    c=pointings['Nobs'].data, marker= 'x',
                    cmap=plt.cm.jet,alpha=.5)

    plt.grid(True)
    plt.title('SDSS field locations')
    cb = plt.colorbar(location='bottom')
        
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    txt = ax.set_xticklabels(xlab)
    for t in txt:
        t.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
    cb.set_label('Number of Field-MJD ')
    plt.savefig(ptt.join(topdir,'SDSSV3.png'),dpi=500,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSV3_s.png'),dpi=50,bbox_inches='tight')
    plt.close()


####################################################################################
def plot_sky_targets(topdir, spall_file, splog, nobs=False, maxn=1000):
    splog.info('Producing Observed Target Density Plots')
    RA1 = None
    if not nobs:
        if not ptt.exists(ptt.join(topdir,'SDSSV2.png')):
            RA1  = np.asarray([np.NaN])
            DEC1 = np.asarray([np.NaN])
            Nobs = np.zeros_like(DEC1)

    if RA1 is None:
        if ptt.exists(spall_file):
            with fits.open(spall_file,memmap=True) as hdul:
                table_data = hdul[1].data
                RA1  = table_data.field('RACAT')
                DEC1 = table_data.field('DECCAT')
                Nobs = table_data.field('NSPECOBS')
                table_data = None
            if len(RA1) == 0:
                RA1 =np.asarray([np.NaN])
                DEC1=np.asarray([np.NaN])
                Nobs=np.zeros_like(DEC1)
        else:
            RA1  = np.asarray([np.NaN])
            DEC1 = np.asarray([np.NaN])
            Nobs = np.zeros_like(DEC1)
            
    RA1 *= np.pi / 180
    DEC1 *= np.pi / 180
    for i in range(0, len(RA1)):
        if RA1[i] >= np.pi:
            RA1[i]=RA1[i]-2*np.pi

    fig, ax = plt.subplots(figsize=(9.5*1.0,5.5*1.25),layout='constrained', subplot_kw={'projection': 'mollweide'})
    
    idx = np.argsort(Nobs)
    maxn = np.nanmax(ma.where(Nobs < np.min([maxn, max(Nobs)]), Nobs, np.nan))
    with np.errstate(invalid="ignore"):
        plt.scatter(RA1[idx], DEC1[idx], s=.05, c=Nobs[idx], cmap=plt.cm.jet,
                    edgecolors='none', linewidths=0, vmin = 1, vmax = maxn)
    plt.grid(True)
    plt.title('SDSS Observed Targets')
    cb = plt.colorbar(location='bottom')
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    txt = ax.set_xticklabels(xlab)
    for t in txt:
        t.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
    cb.set_label('Nobs')
    plt.savefig(ptt.join(topdir,'SDSSV2_s.png'),dpi=50,bbox_inches='tight')
    plt.savefig(ptt.join(topdir,'SDSSV2.png'),dpi=500,bbox_inches='tight')
    plt.close()
