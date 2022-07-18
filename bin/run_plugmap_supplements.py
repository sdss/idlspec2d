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
  
def get_FieldCadence(Field_id, rs_plan):
    try:
        from sdssdb.peewee.sdss5db.targetdb import Field,Version
        tp = (Field.select().join(Version)
                            .where(Field.field_id == Field_id)
                            .where(Version.plan == rs_plan))
        print([(('Fieldid', 'Version_pk', 'RS_tag', 'RS_plan'),
                (t.field_id, t.version.pk, t.version.tag, t.version.plan)) for t in tp])
        if len(tp) > 0: return(tp[0].cadence.label)
        else: return('')
    except: return('')
    return('')

def get_CartonInfo(row):
    try:
        from sdssdb.peewee.sdss5db.targetdb import CartonToTarget, Carton
    except: return(row)
    tp = CartonToTarget.select().join(Carton).where(CartonToTarget.pk == int(row.carton_to_target_pk))
    if len(tp) > 0:
        try: row.program = tp[0].carton.program
        except: pass
        try: row.carton = tp[0].carton.carton
        except: pass
        try: row.CatVersion = tp[0].carton.version.plan #0.5.0
        except: pass
        try: row.mapper = tp[0].carton.mapper.label #"BHM"/"MWM"
        except:
            if row.program == 'open_fiber': row.mapper = 'open_fiber'
            elif 'ops' in row.program: row.mapper = 'ops'
            else: row.mapper = ''
    return (row)
    
def get_mags(row, mags=True, astr=True):
    try:
        from sdssdb.peewee.sdss5db.catalogdb import AllWise, Gaia_DR2, GUVCat
        wise_best=AllWise.select(AllWise.ra,AllWise.dec,AllWise.w1mpro,AllWise.w2mpro, AllWise.w3mpro,AllWise.w4mpro,AllWise.j_m_2mass,AllWise.h_m_2mass, AllWise.k_m_2mass)
        gaia_best=Gaia_DR2.select(Gaia_DR2.ra,Gaia_DR2.dec,Gaia_DR2.parallax, Gaia_DR2.pmra,Gaia_DR2.pmdec)
        guvcat_best=GUVCat.select(GUVCat.ra,GUVCat.dec,GUVCat.fuv_mag,GUVCat.nuv_mag)
        ra=row.ra
        dec=row.dec
        if mags is True:

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
                tp=guvcat_best.where(GUVCat.cone_search(ra, dec, 2.0/3600.0))
                if len(tp) > 0:
                    row.fuv=float(tp[0].fuv_mag)
                    row.nuv=float(tp[0].nuv_mag)
                else: pass
            except: pass
        if astr is True:
            try:
                tp=gaia_best.where(Gaia_DR2.cone_search(ra, dec, 2.0/3600.0))
                if len(tp) > 0:
                    row.parallax=float(tp[0].parallax)
                    row.pmra=float(tp[0].pmra)
                    row.pmdec=float(tp[0].pmdec)
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

    reddening[np.where(np.asarray(rr)<=0.0)[0].tolist()]=np.nan
    data.reddening_gaia=reddening
    return(data)

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description='Runs Plugmap Suppliments')
    parser.add_argument('catalogfile', type=str, help='Catalog file')
    parser.add_argument('--log', '-l',  help='log file', type=str)
    parser.add_argument('--mags', '-m', help='Extra Magnitudes', action='store_true', default=False)
    parser.add_argument('--astrometry', help='Gaia astrometry', action='store_true', default=False)
    parser.add_argument('--rjce', '-r', help='RJCE extintion method', action='store_true', default=False)
    parser.add_argument('--gaia', '-g', help='GAIA extintion method', action='store_true', default=False)
    parser.add_argument('--cart', '-c', help='Get Carton Meta Data', action='store_true', default=False)
    parser.add_argument('--fieldid', help='SDSS-V FieldID', type=int, default=None)
    parser.add_argument('--rs_plan', help='Robostrategy Plan', type=str, default=None)
    
    args = parser.parse_args()
    fp=open(args.catalogfile)
    lines = fp.readlines()

    if args.fieldid is not None and args.rs_plan is not None:
        logstr= "Obtaining Field Cadence"
        print(logstr)
        if args.log is not None: os.system('echo "'+logstr+'" >> '+args.log)
        fieldCadence = get_FieldCadence(args.fieldid, args.rs_plan)
    else: fieldCadence = ''
    data=pd.DataFrame()
    cols=pd.Series({'w1mpro':np.NaN,'w2mpro':np.NaN,'w3mpro':np.NaN,'w4mpro':np.NaN,'j2mass':np.NaN,
                       'h2mass':np.NaN,'k2mass':np.NaN,'fuv':np.NaN,'nuv':np.NaN,'parallax':np.NaN,
                       'pmra':np.NaN,'pmdec':np.NaN, 'EBV_rjce':np.NaN,'reddening_gaia':np.NaN,
                       'program':'', 'carton':'', 'CatVersion':'', 'mapper':'', 'fieldCadence':fieldCadence})
    for line in lines :
        row=pd.Series({'ra':float(line.split()[0]),
                       'dec':float(line.split()[1]),
                       'catid':line.split()[2],
                       'carton_to_target_pk':int(line.split()[3]),
                       'stdflag':int(line.split()[4]),
                       'll':float(line.split()[5]),
                       'bb':float(line.split()[6]),
                       'rr':float(line.split()[7])})
        row=pd.concat([row,cols])
        
        data=data.append(row, ignore_index=True)

    if args.cart is True:
        logstr = "Obtaining the Carton Meta Data"
        print(logstr)
        if args.log is not None: os.system('echo "'+logstr+'" >> '+args.log)
        data=data.apply(get_CartonInfo, axis=1)
        
    if args.mags is True:
        if args.astrometry is True: logstr= "Obtaining the WISE, TWOMASS, GUVCAT Mag and GAIA parallax and pm"
        else: logstr = "Obtaining the WISE, TWOMASS, and GUVCAT Mag"
        print(logstr)
        if args.log is not None: os.system('echo "'+logstr+'" >> '+args.log)
        data=data.apply(get_mags,axis=1,mags=args.mags, astr=args.astrometry)
    elif args.astrometry is True:
        logstr = "Obtaining Gaia parallex and pm"
        print(logstr)
        if args.log is not None: os.system('echo "'+logstr+'" >> '+args.log)
        data=data.apply(get_mags,axis=1,mags=args.mags, astr=args.astrometry)
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
