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
except: print('ERROR: dustmaps is not installed')

try: from dustmaps.sfd import SFDQuery
except: print('ERROR: dustmaps is not installed')

from sdssdb.peewee.sdss5db.targetdb import database
import sdssdb
database.set_profile('operations')


SDSSDBVersion=sdssdb.__version__


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def get_FieldCadence(designID, rs_plan):
    try:
        from sdssdb.peewee.sdss5db.targetdb import Design, Field, Version
        from sdssdb.peewee.sdss5db.targetdb import DesignToField as d2f
                
        field = Field.select().join(d2f).join(Design).switch(Field)\
                     .join(Version).where(Design.design_id == designID)\
                     .where(Version.plan==rs_plan)
        
        print([(('Fieldid', 'Version_pk', 'RS_tag', 'RS_plan'),
                (t.field_id, t.version.pk, t.version.tag, t.version.plan)) for t in field])
        if len(field) > 0: return(field[0].cadence.label)
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

def corrections(row):
    try:
        from sdssdb.peewee.sdss5db.targetdb import RevisedMagnitude as v05_rev_mag
    except: return(row)
    tp = v05_rev_mag.select().where(v05_rev_mag.carton_to_target_pk == int(row.carton_to_target_pk))

    if len(tp) == 1:
        row.mag_g = tp[0].g
        row.mag_r = tp[0].r
        row.mag_i = tp[0].i
        row.mag_z = tp[0].z
        row.mag_j = tp[0].j
        row.mag_h = tp[0].h
        row.mag_k = tp[0].k
        row.gaia_g = tp[0].gaia_g
        row.gaia_bp = tp[0].bp
        row.gaia_rp = tp[0].rp
        row.optical_prov = tp[0].optical_prov
        row.v05_rev_mag = 1
    return(row)

def get_mags(row, mags=True, astr=True, gaia_id=True):
    try:
        from sdssdb.peewee.sdss5db.catalogdb import CatalogToGUVCat, GUVCat
        from sdssdb.peewee.sdss5db.catalogdb import CatalogToAllWise, AllWise

        from sdssdb.peewee.sdss5db.catalogdb import CatalogToTIC_v8
        from sdssdb.peewee.sdss5db.catalogdb import TIC_v8, Gaia_DR2

        from sdssdb.peewee.sdss5db.catalogdb import Gaia_DR2 as Gaia



        wise_best=AllWise.select(AllWise.ra,AllWise.dec,AllWise.w1mpro,AllWise.w2mpro, AllWise.w3mpro,AllWise.w4mpro,AllWise.j_m_2mass,AllWise.h_m_2mass, AllWise.k_m_2mass)
        gaia_best=Gaia_DR2.select(Gaia_DR2.ra,Gaia_DR2.dec,Gaia_DR2.parallax, Gaia_DR2.pmra,Gaia_DR2.pmdec,Gaia_DR2.source_id)
        guvcat_best=GUVCat.select(GUVCat.ra,GUVCat.dec,GUVCat.fuv_mag,GUVCat.nuv_mag)
        ra=row.ra
        dec=row.dec
        if mags is True:

            try:
                tp = wise_best.join(CatalogToAllWise).where(CatalogToAllWise.catalogid == row.catid)
                if len(tp) == 0:
                    tp = wise_best.where(AllWise.cone_search(ra, dec, 2.0/3600.0))
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
                tp = guvcat_best.join(CatalogToGUVCat).where(CatalogToGUVCat.catalogid == row.catid)
                if len(tp) == 0:
                    tp=guvcat_best.where(GUVCat.cone_search(ra, dec, 2.0/3600.0))
                if len(tp) > 0:
                    row.fuv=float(tp[0].fuv_mag)
                    row.nuv=float(tp[0].nuv_mag)
                else: pass
            except: pass
        if astr is True:
            try:
                tp = TIC_v8.select().join(CatalogToTIC_v8).where(CatalogToTIC_v8.catalogid == row.catid)\
                                    .join(Gaia_DR2, on=(TIC_v8.gaia == Gaia_DR2.source_id)).switch(Gaia_DR2)
                
                if len(tp) == 0:
                    tp=gaia_best.where(Gaia_DR2.cone_search(ra, dec, 2.0/3600.0))
                else:
                    tp = [tp[0].gaia]
                if len(tp) > 0:
                    row.parallax=float(tp[0].parallax)
                    row.pmra=float(tp[0].pmra)
                    row.pmdec=float(tp[0].pmdec)
                    if gaia_id is True:
                        row.gaia_id = tp[0].source_id
                else: pass
            except: pass
        elif gaia_id is True:
            try:
                tp = TIC_v8.select().join(CatalogToTIC_v8).where(CatalogToTIC_v8.catalogid == row.catid)\
                                    .join(Gaia_DR2, on=(TIC_v8.gaia == Gaia_DR2.source_id)).switch(Gaia_DR2)
                
                if len(tp) == 0:
                    tp=gaia_best.where(Gaia_DR2.cone_search(ra, dec, 2.0/3600.0))
                else:
                    tp = [tp[0].gaia]
                if len(tp) > 0:
                    row.gaia_id = tp[0].source_id
                else: pass
            except: pass
    except: pass
    return(row)


def get_RJCE_ext(row):
    try:
        from sdssdb.peewee.sdss5db.catalogdb import CatalogToAllWise, AllWise
        catid=np.long(row.catid)
        wise_best=AllWise.select().join(CatalogToAllWise)\
                        .where(CatalogToAllWise.catalogid == catid)\
                        .where(AllWise.designation == CatalogToAllWise.target.designation)
        if len(wise_best) > 0:
            row.EB_rjce=0.918 * (float(wise_best[0].h_m_2mass)-float(wise_best[0].w2mpro) - 0.05)/0.302
        else: pass
    except: pass
    return(row)
 
def get_gaia_red(data, fps=False):
    ll = data.ll.values
    bb = data.bb.values
    
    if fps is False:
        PARALLAX = data.parallax
        PARALLAX[np.where(PARALLAX <= -999.0)[0]] = 0
        PARALLAX[np.where(np.isnan(PARALLAX))[0]] = 0
        dist_std=1.0/np.abs((PARALLAX-0.0)*1e-3) #zero point parallax
        dist_std[np.where(PARALLAX <= -999.0)[0]] = 0
        dist_std[np.where(np.isnan(PARALLAX))[0]] = 0
        rr = dist_std
    else:
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


def get_sfd_red(data):
    ll = data.ll.values
    bb = data.bb.values
    with HiddenPrints():
        try:
            sfd = SFDQuery()
        except ImportError:
            return(data)
        coords = SkyCoord(ll*units.deg, bb*units.deg, frame='galactic')
        reddening = sfd(coords)

    data.EBV_sfd=reddening
    return(data)



def run_plugmap_supplements(catalogfile, log=None, lco=False, mags=False,
                            astrometry=False, rjce=False, gaia=False,
                            sfd = False, gaia_id=False,
                            cart=False, designID=None, rs_plan=None):

    fp=open(catalogfile)
    lines = fp.readlines()
    fps= False
    if designID is not None and rs_plan is not None:
        logstr= "Obtaining Field Cadence"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        fieldCadence = get_FieldCadence(designID, rs_plan)
        fps = True
    else: fieldCadence = ''
    data=pd.DataFrame()
    cols=pd.Series({'w1mpro':np.NaN,'w2mpro':np.NaN,'w3mpro':np.NaN,'w4mpro':np.NaN,'j2mass':np.NaN,
                    'h2mass':np.NaN,'k2mass':np.NaN,'fuv':np.NaN,'nuv':np.NaN,'parallax':np.NaN,
                    'pmra':np.NaN,'pmdec':np.NaN, 'gaia_id':-1,
                    'EBV_rjce':np.NaN,'reddening_gaia':np.NaN,'EBV_sfd':0,
                    'program':'', 'carton':'', 'CatVersion':'', 'mapper':'', 'fieldCadence':fieldCadence,
                    'mag_g':np.NaN,'mag_r':np.NaN, 'mag_i':np.NaN,'mag_z':np.NaN,
                    'mag_j':np.NaN,'mag_h':np.NaN, 'mag_k':np.NaN,
                    'gaia_g':np.NaN,'gaia_bp':np.NaN,'gaia_rp':np.NaN,'optical_prov':'', 'v05_rev_mag':0,
                    })

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

    if cart is True:
        logstr = "Obtaining the Carton Meta Data"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=data.apply(get_CartonInfo, axis=1)
        
    if mags is True:
        if astrometry is True: logstr= "Obtaining the WISE, TWOMASS, GUVCAT Mag and GAIA ID, parallax and pm"
        else: logstr = "Obtaining the WISE, TWOMASS, and GUVCAT Mag"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=data.apply(get_mags,axis=1,mags=mags, astr=astrometry, gaia_id=gaia_id)
    elif astrometry is True:
        logstr = "Obtaining Gaia ID, parallex and pm"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=data.apply(get_mags,axis=1,mags=mags, astr=astrometry, gaia_id=gaia_id)
    elif gaia_id is True:
        logstr = "Obtaining Gaia ID"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=data.apply(get_mags,axis=1,mags=mags, astr=astrometry, gaia_id=gaia_id)

    if rjce is True:
        logstr = "Defining the Extintion using the RJCE extintion method"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=data.apply(get_RJCE_ext,axis=1)
    if gaia is True:
        logstr = "Defining the Extintion using the Bayestar 3D dust extintion maps"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=get_gaia_red(data, fps=fps)
    if sfd is True:
        logstr = "Defining the Extintion using the SFD dust extintion maps"
        print(logstr)
        if log is not None: os.system('echo "'+logstr+'" >> '+log)
        data=get_sfd_red(data)

    data = data.apply(corrections, axis=1)
    filename = Path(Path(catalogfile).stem+'_supp')

    t = Table.from_pandas(data)
    t.write(filename.with_suffix('.fits'), overwrite=True)
    print("run_plugmap_supplements Complete")



if __name__ == '__main__' :
    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description='Runs Plugmap Suppliments')
    parser.add_argument('catalogfile', type=str, help='Catalog file')
    parser.add_argument('--log', '-l',  help='log file', type=str)
    parser.add_argument('--lco', help='LCO Observatory', action='store_true', default=False)
    parser.add_argument('--mags', '-m', help='Extra Magnitudes', action='store_true', default=False)
    parser.add_argument('--astrometry', help='Gaia astrometry', action='store_true', default=False)
    parser.add_argument('--id_gaia', '-i',help='Get Gaia ID', action='store_true', default=False)
    parser.add_argument('--rjce', '-r', help='RJCE extintion method', action='store_true', default=False)
    parser.add_argument('--gaia', '-g', help='GAIA extintion method', action='store_true', default=False)
    parser.add_argument('--sfd', '-s', help='SFD extintion method', action='store_true', default=False)
    parser.add_argument('--cart', '-c', help='Get Carton Meta Data', action='store_true', default=False)
    parser.add_argument('--designID', '-d', help='SDSS-V DesignID', type=int, default=None)
    parser.add_argument('--rs_plan', help='Robostrategy Plan', type=str, default=None)
    
    args = parser.parse_args()

    run_plugmap_supplements(args.catalogfile, log=args.log, lco=args.lco, mags=args.mags,
                            astrometry=args.astrometry, rjce=args.rjce, gaia=args.gaia,
                            gaia_id=args.id_gaia, sfd=args.sfd, 
                            cart=args.cart, designID= args.designID, rs_plan=args.rs_plan)
