#!/usr/bin/env python
import sys
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import SkyCoord
import astropy.units as units
import pandas as pd
from pathlib import Path
import argparse
from astropy.table import Table


try: from dustmaps.bayestar import BayestarQuery
except: print('WARNING: dustmaps is not installed')

from sdssdb.peewee.sdss5db.catalogdb import database
database.set_profile('operations')

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
  
def get_mags(row):
    try:
        from sdssdb.peewee.sdss5db.catalogdb import AllWise, Gaia_DR2, GUVCat
        wise_best=AllWise.select(AllWise.ra,AllWise.dec,AllWise.w1mpro,AllWise.w2mpro, AllWise.w3mpro,AllWise.w4mpro,AllWise.j_m_2mass,AllWise.h_m_2mass, AllWise.k_m_2mass)
        gaia_best=Gaia_DR2.select(Gaia_DR2.ra,Gaia_DR2.dec,Gaia_DR2.parallax, Gaia_DR2.pmra,Gaia_DR2.pmdec)
        guvcat_best=GUVCat.select(GUVCat.ra,GUVCat.dec,GUVCat.fuv_mag,GUVCat.nuv_mag)
        ra=row.ra
        dec=row.dec
        try:
            tp=wise_best.where(AllWise.cone_search(ra, dec, 2.0/3600.0))
            if len(tp) > 0:
                row.w1mpro=float(tp[0].w1mpro)
                row.w2mpro=float(tp[0].w2mpro)
                row.w3mpro=float(tp[0].w3mpro)
                row.w4mpro=float(tp[0].w4mpro)
                try:
                    row.j2mass=float(tp[0].j_m_2mass)
                    row.h2mass=float(tp[0].h_m_2mass)
                    row.k2mass=float(tp[0].k_m_2mass)
                except: pass
            else: pass
        except: pass
        try:
            tp=gaia_best.where(Gaia_DR2.cone_search(ra, dec, 2.0/3600.0))
            if len(tp) > 0:
                row.parallax=float(tp[0].parallax)
                row.pmra=float(tp[0].pmra)
                row.pmdec=float(tp[0].pmdec)
            else: pass
        except: pass
        try:
            tp=guvcat_best.where(GUVCat.cone_search(ra, dec, 2.0/3600.0))
            if len(tp) > 0:
                row.fuv=float(tp[0].fuv_mag)
                row.nuv=float(tp[0].nuv_mag)
            else: pass
        except: pass
    except: pass
    return(row)


def get_RJCE_ext(row):
    try:
        from sdssdb.peewee.sdss5db.catalogdb import CatalogToAllWise, AllWise
        catid=np.long(row.catid)
        wise_best=AllWise.select().join(CatalogToAllWise).where(CatalogToAllWise.catalogid == catid).where(AllWise.designation == CatalogToAllWise.target.designation)
        if len(wise_best) > 0:
            row.EB_rjce=0.918 * (float(wise_best_v2[0].h_m_2mass)-float(wise_best_v2[0].w2mpro) - 0.05)/0.302
        else: pass
    except: pass
    return(row)
 
def get_gaia_red(data):
    ll = data.ll.values
    bb = data.bb.values
    rr = data.rr.values
    with HiddenPrints():
        try:
            bayestar = BayestarQuery(version='bayestar2015')
        except FileNotFoundError:
            try:
                import dustmaps.bayestar
                dustmaps.bayestar.fetch(version='bayestar2015')
                bayestar = BayestarQuery(version='bayestar2015')
            except ImportError: return(data)
        coords = SkyCoord(ll*units.deg, bb*units.deg,distance=rr*units.pc, frame='galactic')
        reddening = bayestar(coords, mode='median')

    reddening[np.where(np.asarray(reddening)<=0.0)[0].tolist()]=np.nan
    data.reddening_gaia=reddening
    return(data)

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description='Runs Plugmap Suppliments')
    parser.add_argument('catalogfile', type=str, help='Catalog file')
    parser.add_argument('--log', '-l',  help='log file', type=str)
    parser.add_argument('--mags', '-m', help='Extra Magnitudes', action='store_true', default=False)
    parser.add_argument('--rjce', '-r', help='RJCE extintion method', action='store_true', default=False)
    parser.add_argument('--gaia', '-g', help='GAIA extintion method', action='store_true', default=False)
    
    args = parser.parse_args()
    fp=open(args.catalogfile)
    lines = fp.readlines()

    data=pd.DataFrame()
    for line in lines :
        row=pd.Series({'ra':float(line.split()[0]),
                       'dec':float(line.split()[1]),
                       'catid':line.split()[2],
                       'stdflag':int(line.split()[3]),
                       'll':float(line.split()[4]),
                       'bb':float(line.split()[5]),
                       'rr':float(line.split()[6])})
        data=data.append(row, ignore_index=True)


    magcols= {'w1mpro':-999,'w2mpro':-999,'w3mpro':-999,'w4mpro':-999,'j2mass':-999,
              'h2mass':-999,'k2mass':-999,'fuv':-999,'nuv':-999,'parallax':-99999,
              'pmra':-99999,'pmdec':-99999}
    for col in magcols.keys(): data[col]=magcols[col]
    data['EB_rjce']=-100
    data['reddening_gaia']=np.NaN

    if args.mags is True:
        print("Obtaing the WISE, TWOMASS, GUVCAT and GAIA parallax and pm")
        if args.log is not None:
            os.system('echo "Obtaing the WISE, TWOMASS, GUVCAT and GAIA parallax and pm" >> '+args.log)
        data=data.apply(get_mags,axis=1)
    if args.rjce is True:
        print("Defining the Extintion using the RJCE extintion method")
        if args.log is not None:
            os.system('echo "Defining the Extintion using the RJCE extintion method" >> '+args.log)
        data=data.apply(get_RJCE_ext,axis=1)
    if args.gaia is True:
        print("Defining the Extintion using the Bayestar 3D dust extintion maps")
        if args.log is not None:
            os.system('echo "Defining the Extintion using the Bayestar 3D dust extintion maps" >> '+args.log)
        data=get_gaia_red(data)
    filename = Path(Path(args.catalogfile).stem+'_supp')

    t = Table.from_pandas(data)
    t.write(filename.with_suffix('.fits'), overwrite=True)
