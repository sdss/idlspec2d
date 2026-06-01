from boss_drp.utils.splog import splog
from boss_drp import MOUNTAIN, database_profile
from sdss_access.path import Path
from sdss_access import Access

from astropy.table import Table, join, Column, unique

import warnings 
import os.path as ptt
import duckdb
import os
import numpy as np
import time



if not MOUNTAIN:
    try:
        from sdssdb.peewee.sdss5db.targetdb import database
        import sdssdb
        splog.add_external_handlers(sdssdb.log.name)
        test = database.set_profile(database_profile)

        if not test:
            splog.info('WARNING: No SDSSDB access - Defaulting to no_db')
            no_db_poss = True
        else:
            SDSSDBVersion = os.getenv('SDSSDB_VER',None)
            if SDSSDBVersion is None:
                SDSSDBVersion=sdssdb.__version__
    except:
        splog.info('WARNING: No SDSSDB access - Defaulting to no_db')
        no_db_poss = True
    else:
        no_db_poss = False
    #try:
    from sdss_semaphore.targeting import TargetingFlags
    try:
        from  sdss_semaphore.targeting import logger as sem_log
        splog.add_external_handlers(sem_log.name)
    except:
        pass

else:
    no_db_poss=True

if not no_db_poss:
    from sdssdb.peewee.sdss5db.catalogdb import CatalogToGUVCat, GUVCat
    from sdssdb.peewee.sdss5db.catalogdb import CatalogToAllWise, AllWise
    from sdssdb.peewee.sdss5db.catalogdb import CatalogToTIC_v8
    from sdssdb.peewee.sdss5db.catalogdb import TIC_v8, Gaia_DR2
    from sdssdb.peewee.sdss5db.catalogdb import Gaia_DR3
    from sdssdb.peewee.sdss5db.catalogdb import CatalogToGaia_DR3 as CatToGaia_DR3
    from sdssdb.peewee.sdss5db.catalogdb import CatalogToTwoMassPSC as C2TM, TwoMassPSC
    from sdssdb.peewee.sdss5db.catalogdb import SDSS_ID_flat
    from sdssdb.peewee.sdss5db.catalogdb import SDSS_ID_stacked


    from sdssdb.peewee.sdss5db.targetdb import Design, Field, Version
    from sdssdb.peewee.sdss5db.targetdb import DesignToField as d2f
    from sdssdb.peewee.sdss5db.targetdb import RevisedMagnitude
    from sdssdb.peewee.sdss5db.targetdb import CartonToTarget, Carton, Version, Mapper, Target

def get_Catalog(catalog, no_remote=False, release='sdsswork', **kwrds):
    kwrds['ftype'] = 'parquet'
    kwrds['num'] = '*' #left here because it will not break anything, but was required at one point
    kwrds['v_targ'] = kwrds.pop('V_TARG','*')
    path   = Path(release=release, preserve_envvars=True)
    access = Access(release=release)#, preserve_envvars=True)

    cats = []
    if 'v_targ' in kwrds:
        if kwrds['v_targ'] == '*':
            max_version = '*'
            try:
                versions = [path.extract(catalog, x)['v_targ'] for x in path.expand(catalog, **kwrds)]
                max_version = max(versions, key=lambda v: tuple(map(int, v.split("."))))
            except:
                max_version = path.extract(catalog, path.expand(catalog, **kwrds)[-1])['v_targ']
            kwrds['v_targ'] = max_version
        
    for pt in path.expand(catalog, **kwrds):
        tkwrds = path.extract(catalog, pt)
        try:
            path.exists(catalog, **tkwrds)
        except:
            print(catalog, tkwrds)
            if 'num' not in tkwrds: tkwrds['num'] = '*'   #left here because it will not break anything, but was required at one point    
        if path.exists(catalog, **tkwrds):
            cats.append(path.full(catalog, **tkwrds))
        elif (not no_remote):
            if (path.exists(catalog, **tkwrds, remote=True)):
                tcat = path.full(catalog, **tkwrds)
                access.remote()
                access.add(catalog, **tkwrds)
                access.set_stream()
                valid = access.commit()
                if valid:
                    cats.append(tcat)
                else:
                    splog.info('ERROR: Cannot find/get'+ptt.basename(tcat))
                    exit()
            else:
                tcat = path.full(catalog, **tkwrds)
                splog.info('ERROR: Cannot find/get'+ptt.basename(tcat))
                exit()
        else:
            tcat = path.full(catalog, **tkwrds)
            splog.info('ERROR: Cannot find/get'+ptt.basename(tcat))
            exit()
    return(cats)


def get_mags_astrom(search_table, db = True, fps=False, fast=False, release='sdsswork', no_remote=False,V_TARG='*'):
    gaia = False
    GUV = False
    allwise = False
    twomass = False

    splog.info('Getting Magnitudes, IDs, and Astrometry')
    if len(search_table[search_table['icatalogid']!= -999]) == 0:
        return(search_table)
    u_s_table = search_table[search_table['icatalogid']!= -999].group_by('icatalogid')
    u_s_table = u_s_table[u_s_table.groups.indices[:-1]]
    catalogids = np.unique(search_table['icatalogid'].data).tolist()
    a_catalogids = np.asarray(catalogids)
    sdssids = np.unique(u_s_table['SDSS_ID'].data).tolist()

    while True:
        try:
            catalogids.remove(0)
        except:
            break


    if fps is True:
        results = Table(names = ('icatalogid','gaia_id','j2mass','h2mass','k2mass'),
                        dtype = (int,int,float,float,float))
        gaia_cols = ['gaia_id']
    else:
        results = Table(names = ('icatalogid','parallax','pmra','pmdec','gaia_id','j2mass','h2mass','k2mass'),
                        dtype = (int,float,float,float,int,float,float,float))
        gaia_cols = ['parallax','pmra','pmdec','gaia_id']

    if db:
        # Get Gaia and Twomass
        tp = SDSS_ID_flat.select(CatToGaia_DR3.catalogid, SDSS_ID_flat.sdss_id, \
                                Gaia_DR3.parallax,Gaia_DR3.pmra,Gaia_DR3.pmdec,Gaia_DR3.source_id.alias('gaia_id'), \
                                TwoMassPSC.j_m.alias('j2mass'),TwoMassPSC.h_m.alias('h2mass'),TwoMassPSC.k_m.alias('k2mass'))\
                         .join(CatToGaia_DR3, on=(SDSS_ID_flat.catalogid == CatToGaia_DR3.catalogid)).join(Gaia_DR3).switch(SDSS_ID_flat)\
                         .join(C2TM, on=(SDSS_ID_flat.catalogid == C2TM.catalogid)).join(TwoMassPSC).switch(SDSS_ID_flat)\
                         .where(SDSS_ID_flat.sdss_id.in_(sdssids)).dicts()
        
    else:
        splog.warning('Getting GAIA DR3 and TWOMASS mags from SDSS-V MOS Targeting Product')

        sdss_id_flat_path = get_Catalog('mos_target_sdss_id_flat', release=release, V_TARG=V_TARG, no_remote=no_remote)
        cat_to_gaia_path = get_Catalog('mos_target_catalog_to_gaia_dr3_source', release=release, V_TARG=V_TARG, no_remote=no_remote)
        gaia_dr3_path  = get_Catalog('mos_target_gaia_dr3_source', release=release, V_TARG=V_TARG, no_remote=no_remote)
        c2tm_path = get_Catalog('mos_target_catalog_to_twomass_psc', release=release, V_TARG=V_TARG, no_remote=no_remote)
        twomass_path  = get_Catalog('mos_target_twomass_psc', release=release, V_TARG=V_TARG, no_remote=no_remote)
        tp = duckdb.execute("""
                SELECT
                    cgaia.catalogid, s.sdss_id, g.parallax, g.pmra, g.pmdec, g.source_id AS gaia_id,
                    tm.j_m AS j2mass, tm.h_m AS h2mass, tm.k_m AS k2mass
                FROM read_parquet(?) s
                JOIN read_parquet(?) cgaia
                    ON s.catalogid = cgaia.catalogid
                JOIN read_parquet(?) g
                    ON cgaia.target_id = g.source_id
                JOIN read_parquet(?) c2tm
                    ON s.catalogid = c2tm.catalogid
                JOIN read_parquet(?) tm
                    ON c2tm.target_id = tm.pts_key
                WHERE s.sdss_id IN ?
            """, [
                sdss_id_flat_path,
                cat_to_gaia_path,
                gaia_dr3_path,
                c2tm_path,
                twomass_path,
                sdssids
            ]).fetchdf().to_dict('records')

    for t in tp:
        for key in t.keys():
            if t[key] is None:
                if key in ['parallax','pmra','pmdec','j2mass','h2mass','k2mass']:
                    t[key] = np.nan
                elif key in ['gaia_id']:
                    t[key] = -999
        cid = u_s_table[u_s_table['SDSS_ID'] == t['sdss_id']]['icatalogid'][0]
        if fps is True:
            results.add_row((cid,int(t['gaia_id']),float(t['j2mass']),float(t['h2mass']),float(t['k2mass'])))
        else:
            results.add_row((cid,float(t['parallax']),float(t['pmra']),float(t['pmdec']),int(t['gaia_id']),
                                float(t['j2mass']),float(t['h2mass']),float(t['k2mass'])))
                
        # Get Gaia and Twomass for pre-v1 targets
    if db:
        tp = CatalogToTIC_v8.select(CatalogToTIC_v8.catalogid, CatalogToTIC_v8.best, \
                                    Gaia_DR2.parallax,Gaia_DR2.pmra,Gaia_DR2.pmdec,Gaia_DR2.source_id.alias('gaia_id'),\
                                    TIC_v8.jmag.alias('j2mass'), TIC_v8.hmag.alias('h2mass'), TIC_v8.kmag.alias('k2mass'))\
                        .join(TIC_v8).join(Gaia_DR2, on=(TIC_v8.gaia == Gaia_DR2.source_id)).switch(CatalogToTIC_v8)\
                        .where(CatalogToTIC_v8.catalogid.in_(catalogids)).dicts()

    else:
        splog.warning('Getting TIC_v8 mags from SDSS-V MOS Targeting Product')
        catalog_to_tic_path = get_Catalog('mos_target_catalog_to_tic_v8', release=release, V_TARG=V_TARG, no_remote=no_remote)
        tic_path = get_Catalog('mos_target_tic_v8', release=release, V_TARG=V_TARG, no_remote=no_remote)
        gaia_path = get_Catalog('mos_target_gaia_dr2_source', release=release, V_TARG=V_TARG, no_remote=no_remote)
        tp = duckdb.execute("""
                SELECT
                    c.catalogid, c.best, g.parallax, g.pmra, g.pmdec, g.source_id AS gaia_id,
                    t.jmag AS j2mass, t.hmag AS h2mass, t.kmag AS k2mass
                FROM read_parquet(?) c
                JOIN read_parquet(?) t
                    ON c.target_id = t.id
                JOIN read_parquet(?) g
                    ON t.gaia = g.source_id
                WHERE c.catalogid IN ?
            """, [catalog_to_tic_path, tic_path, gaia_path, catalogids]).fetchdf().to_dict('records')
        
    for t in tp:
        if t['best'] is False: continue
        for key in t.keys():
            if t[key] is None:
                if key in ['parallax','pmra','pmdec','j2mass','h2mass','k2mass']:
                    t[key] = np.nan
                elif key in ['gaia_id']:
                    t[key] = -999
        try:
            t['catalogid']
        except:
            t['catalogid'] = t['catalog']
        if t['catalogid'] in results['icatalogid'].data:
            # Check if the is a matching row and update missing values
            row = results[results['icatalogid'] == t['catalogid']]
            if np.isnan(row['j2mass'][0]) and np.isnan(row['h2mass'][0]) and np.isnan(row['k2mass'][0]):
                for key in ['j2mass','h2mass','k2mass']:
                    results[results['icatalogid'] == t['catalogid']][key] = t[key]
            if row['gaia_id'][0] == -999:
                for key in gaia_cols:
                    results[results['icatalogid'] == t['catalogid']][key] = t[key]
            continue
        # catalogid is not in results yet
        if fps is True:
            results.add_row((t['catalogid'],int(t['gaia_id']),
                            float(t['j2mass']),float(t['h2mass']),float(t['k2mass'])))
        else:
            results.add_row((t['catalogid'],float(t['parallax']),float(t['pmra']),float(t['pmdec']),int(t['gaia_id']),
                            float(t['j2mass']),float(t['h2mass']),float(t['k2mass'])))
    

    if len(results) > 0:
        gaia = True
        twomass = True
        search_table = join(search_table,results,keys='icatalogid',join_type='left')

    if not fast:
        if db:
            tp = CatalogToGUVCat.select(CatalogToGUVCat.catalogid, CatalogToGUVCat.best, GUVCat.fuv_mag, GUVCat.nuv_mag)\
                            .join(GUVCat).switch(CatalogToGUVCat)\
                            .where(CatalogToGUVCat.catalogid.in_(catalogids)).dicts()
        else:
            splog.warning('Getting Allwise mags from SDSS-V MOS Targeting Product')

            catalog_to_guv_path = get_Catalog('mos_target_catalog_to_guvcat', release=release, V_TARG=V_TARG, no_remote=no_remote)
            guv_path = get_Catalog('mos_target_guvcat', release=release, V_TARG=V_TARG, no_remote=no_remote)
            tp = duckdb.execute("""
                    SELECT c.catalogid, c.best, g.fuv_mag, g.nuv_mag
                    FROM read_parquet(?) c
                    JOIN read_parquet(?) g
                        ON c.target_id = g.objid
                    WHERE c.catalogid IN ?
                """, [catalog_to_guv_path, guv_path, catalogids]).fetchdf().to_dict('records')
        results = Table(names = ('icatalogid','fuv','nuv'), dtype=(int,float,float))
        for t in tp:
            if t['best'] is False: continue
            for key in t.keys():
                if t[key] is None:
                    if key in ['fuv_mag','nuv_mag']:
                        t[key] = np.nan
            try:
                t['catalogid']
            except:
                t['catalogid'] = t['catalog']
            results.add_row((t['catalogid'],float(t['fuv_mag']),float(t['nuv_mag'])))
        if len(results) > 0:
            GUV = True
            search_table = join(search_table,results,keys='icatalogid', join_type='left')
        
    if not fast: 
        if db:
            tp = CatalogToAllWise.select(CatalogToAllWise.catalogid, CatalogToAllWise.best, AllWise.w1mpro, AllWise.w2mpro,
                                            AllWise.w3mpro, AllWise.w4mpro)\
                            .join(AllWise).switch(CatalogToAllWise)\
                            .where(CatalogToAllWise.catalogid.in_(catalogids)).dicts()
        else:
            splog.warning('Getting Allwise mags from SDSS-V MOS Targeting Product')

            catalog_to_allwise_path = get_Catalog('mos_target_catalog_to_allwise', release=release, V_TARG=V_TARG, no_remote=no_remote)
            allwise_path = get_Catalog('mos_target_allwise', release=release, V_TARG=V_TARG, no_remote=no_remote)
            tp = duckdb.execute("""
                        SELECT c.catalogid, c.best, w.w1mpro, w.w2mpro, w.w3mpro, w.w4mpro
                        FROM read_parquet(?) c
                        JOIN read_parquet(?) w
                            ON c.target_id = w.cntr
                        WHERE c.catalogid IN ?
                    """, [catalog_to_allwise_path, allwise_path, catalogids]).fetchdf().to_dict('records')
        results = Table(names = ('icatalogid','w1mpro','w2mpro','w3mpro','w4mpro'),
                        dtype=(int,float,float,float,float))
        for t in tp:
            if t['best'] is False: continue
            for key in t.keys():
                if t[key] is None:
                    if key in ['w1mpro','w2mpro','w3mpro','w4mpro']:
                        t[key] = np.nan
            try:
                t['catalogid']
            except:
                t['catalogid'] = t['catalog']
            results.add_row((t['catalogid'],float(t['w1mpro']),float(t['w2mpro']),
                                            float(t['w3mpro']),float(t['w4mpro'])))
        
    if len(results) > 0:
        allwise = True
        search_table = join(search_table,results,keys='icatalogid', join_type='left')
  
    if allwise is True:
        mag = search_table['WISE_MAG']
        mag[:,0] = search_table['w1mpro'].data.filled(fill_value=np.nan)
        mag[:,1] = search_table['w2mpro'].data.filled(fill_value=np.nan)
        mag[:,2] = search_table['w3mpro'].data.filled(fill_value=np.nan)
        mag[:,3] = search_table['w4mpro'].data.filled(fill_value=np.nan)
        search_table['WISE_MAG'] = mag
        search_table.remove_columns(['w1mpro','w2mpro','w3mpro','w4mpro'])

    if twomass is True:
        mag = search_table['TWOMASS_MAG']
        mag[:,0] = search_table['j2mass'].data.filled(fill_value=np.nan)
        mag[:,1] = search_table['h2mass'].data.filled(fill_value=np.nan)
        mag[:,2] = search_table['k2mass'].data.filled(fill_value=np.nan)
        search_table['TWOMASS_MAG'] = mag
        search_table.remove_columns(['j2mass','h2mass','k2mass'])

    if GUV is True:
        mag = search_table['GUVCAT_MAG']
        mag[:,0] = search_table['fuv'].data.filled(fill_value=np.nan)
        mag[:,1] = search_table['nuv'].data.filled(fill_value=np.nan)
        search_table['GUVCAT_MAG'] = mag
        search_table.remove_columns(['fuv','nuv'])
    
    return(search_table)

def get_FieldCadence(designID, rs_plan, db=True,release='sdsswork', V_TARG='*', no_remote=False):
    splog.info("Obtaining Field Cadence")
    if db:
        field = Field.select().join(d2f).join(Design).switch(Field)\
                        .join(Version).switch(Field).where(Design.design_id == designID)\
                        .where(Version.plan==rs_plan)
    else:
        splog.warning('Obtaining Field Cadence from SDSS-V MOS Targeting Product')
        field_path = get_Catalog('mos_target_field', release=release, V_TARG=V_TARG, no_remote=no_remote)
        d2f_path = get_Catalog('mos_target_design_to_field', release=release, V_TARG=V_TARG, no_remote=no_remote)
        design_path = get_Catalog('mos_target_design', release=release, V_TARG=V_TARG, no_remote=no_remote)
        version_path = get_Catalog('mos_target_targetdb_version', release=release, V_TARG=V_TARG, no_remote=no_remote)


        field = duckdb.execute("""
            SELECT f.*
            FROM read_parquet(?) f
            JOIN read_parquet(?) d2f ON f.pk = d2f.field_pk
            JOIN read_parquet(?) d   ON d2f.design_id = d.design_id
            JOIN read_parquet(?) v   ON f.version_pk = v.pk
            WHERE d.design_id = ?
            AND v.plan = ?
        """, [
            field_path,
            d2f_path,
            design_path,
            version_path,
            designID,
            rs_plan
        ]).fetchdf()


    if len(field) > 0:
            t = field[0]
            obsmode = t.cadence.obsmode_pk
            if obsmode is not None:
                obsmode = field[0].cadence.obsmode_pk[0]
            else:
                obsmode = ''
            splog.info(f'Fieldid: {t.field_id}'+'\n'+
                  f'    Version_pk:    {t.version.pk}'+'\n'+
                  f'    RS_tag:        {t.version.tag}'+'\n'+
                  f'    RS_plan:       {t.version.plan}'+'\n'+
                  f'    Field Cadence: {t.cadence.label}'+'\n'+
                  f'    ObsMode:       {obsmode}'
                 )
    elif (str(designID).strip() != '-999') & (str(rs_plan).strip().upper() != 'NA'):
        splog.info(f'Warning: No matching Field found for DesignID ({designID}) and RS_plan ({rs_plan})')
    else:
        splog.info(f'Warning: Invalid DesignID ({designID}) or RS_plan ({rs_plan})')    
    design = Design.select().where(Design.design_id == designID)
    design = design.dicts()
    if len(design) > 0:
        designmode = design[0]['design_mode']
    else:
        designmode = None
    if designmode is None:
        designmode = ''
        if str(designID).strip() != '-999':
            splog.info(f'Warning: No Design Mode found for DesignID ({designID})')
    if len(field) > 0:
        obsmode = field[0].cadence.obsmode_pk
        if obsmode is not None:
            obsmode = field[0].cadence.obsmode_pk[0]
        else:
            obsmode = ''
        return(field[0].cadence.label, obsmode, designmode)
    return('','','')

def target_tab_correction(search_table, db = True, release='sdsswork', V_TARG='*', no_remote=False):
    carton_to_target_pk = search_table['carton_to_target_pk'].data.tolist()

    if db is True:
        splog.info('Checking RevisedMagnitude Table')    
        tp = RevisedMagnitude.select().where(RevisedMagnitude.carton_to_target_pk.in_(carton_to_target_pk)).dicts()
    else:
        splog.info('Checking RevisedMagnitude Table from SDSS-V MOS Targeting Product')
        parquet_file = get_Catalog('mos_target_revised_magnitude', release=release, V_TARG=V_TARG, no_remote=no_remote)

        tp = duckdb.execute("""
                        SELECT carton_to_target_pk as carton_to_target, 
                            g, r, i, z, j, h, k, gaia_g, bp, rp, optical_prov
                        FROM read_parquet(?)
                        WHERE carton_to_target_pk IN ?
                """, [parquet_file, carton_to_target_pk]).fetchdf().to_dict('records')

    results = Table(names = ('carton_to_target_pk','mag_g','mag_r','mag_i','mag_z','mag_j','mag_h','mag_k',
                                'gaia_g','gaia_bp','gaia_rp','optical_prov_rev','v05_rev_mag'),
                    dtype = (int, float, float, float, float, float, float, float, float, float, float, object, bool))
    for t in tp:
        for key in t.keys():
            if t[key] is None:
                if key in ['g','r','i','z','j','h','k','gaia_g','bp','rp']:
                    t[key] = np.nan
                elif key in ['optical_prov']:
                    t[key] = ''
        results.add_row((t['carton_to_target'],float(t['g']),float(t['r']),float(t['i']),float(t['z']),
                            float(t['j']),float(t['h']),float(t['k']),float(t['gaia_g']),float(t['bp']),float(t['rp']),
                            t['optical_prov'], True))
    
    if len(results) > 0:
        splog.info('Updating Magnitudes from RevisedMagnitudes')
        search_table = join(search_table,results,keys='carton_to_target_pk', join_type='left')
        
        
        
        mag = search_table['mag'].data
        corrected = np.where(search_table['v05_rev_mag'].data == True)[0]
        splog.info(f'Updating {len(corrected)} rows')
        mag[corrected,1] = search_table['mag_g'].data[corrected]
        mag[corrected,2] = search_table['mag_r'].data[corrected]
        mag[corrected,3] = search_table['mag_i'].data[corrected]
        mag[corrected,4] = search_table['mag_z'].data[corrected]
        search_table['mag'] = mag

        magt = search_table['bp_mag']
        magt[corrected] = search_table['gaia_bp'].data[corrected]
        search_table['bp_mag'] = magt

        magt = search_table['rp_mag']
        magt[corrected] = search_table['gaia_rp'].data[corrected]
        search_table['rp_mag'] = magt

        magt = search_table['gaia_g_mag']
        magt[corrected] = search_table['gaia_g'].data[corrected]
        search_table['gaia_g_mag'] = magt

        magt = search_table['h_mag']
        magt[corrected] = search_table['mag_h'].data[corrected]
        search_table['h_mag'] = magt

        magt = search_table['optical_prov']
        magt[corrected] = search_table['optical_prov_rev'].data[corrected]
        search_table['optical_prov'] = magt

    return(search_table)


def get_SDSSID(search_table, db=True, release='sdsswork', V_TARG='*', no_remote=False):
    splog.info('Getting SDSS_ID')
    catalogids = np.unique(search_table['icatalogid'].data).tolist()
    if db is True:        
        try:
            tp = SDSS_ID_flat.select(SDSS_ID_flat.catalogid, SDSS_ID_flat.sdss_id)\
                             .where(SDSS_ID_flat.catalogid.in_(catalogids)).dicts()
        except:
            splog._log.exception('Error getting SDSS_ID, trying again....')
            time.sleep(60)
            tp = SDSS_ID_flat.select(SDSS_ID_flat.catalogid, SDSS_ID_flat.sdss_id)\
                             .where(SDSS_ID_flat.catalogid.in_(catalogids)).dicts()
            
    else:
        splog.warning('Getting SDSS_IDs from SDSS-V MOS Targeting Product')
        parquet_file = get_Catalog('mos_target_sdss_id_flat', release=release, V_TARG=V_TARG, no_remote=no_remote)

        tp = duckdb.execute("""
            SELECT catalogid, sdss_id
            FROM read_parquet(?)
            WHERE catalogid IN ?
        """, [parquet_file, catalogids]).fetchdf().to_dict('records')

    results = Table(names=('icatalogid','SDSS_ID'), dtype=(int,int))
    for t in tp:
        results.add_row((t['catalogid'],t['sdss_id']))

    if len(results) == 0:
        splog.info('Warning: No SDSS_ID matches found - Setting all SDSS_ID to -999')
        search_table['SDSS_ID'] = -999
        return(search_table)
    results.sort(['SDSS_ID'])
    results = unique(results, keys='icatalogid', keep='first')
    if len(results) > 0:
        search_table = join(search_table, results, keys='icatalogid',join_type='left')
    else:
        splog.info('Warning: No SDSS_ID matches found - Setting all SDSS_ID to -999')
        search_table['SDSS_ID'] = -999
    try:
        search_table['SDSS_ID'] = search_table['SDSS_ID'].filled(-999)
    except:
        pass

    sci_idx = search_table['category'] == 'science'
    if np.any(sci_idx):
        n_no_match = np.sum(search_table['SDSS_ID'][sci_idx] == -999)
        if n_no_match > 0:
            splog.info('Warning: SDSS_IDs not found for {} science targets'.format(n_no_match))

    return(search_table)


def get_AltCatids(search_table, db=True, release='sdsswork', V_TARG='*', no_remote=False):
    splog.info('Getting All Catalogids for SDSS_IDs')
    sdssids = np.unique(search_table['SDSS_ID'].data).tolist()
    if db is True:
        tp = SDSS_ID_stacked.select()\
                            .where(SDSS_ID_stacked.sdss_id.in_(sdssids)).dicts()
        
    else:
        splog.warning('Getting All Catalogids for SDSS_IDs from SDSS-V MOS Targeting Product')
        parquet_file = get_Catalog('mos_target_sdss_id_stacked', release=release, V_TARG=V_TARG, no_remote=no_remote)

        tp = duckdb.execute("""
            SELECT sdss_id, catalogid21, catalogid25, catalogid31
            FROM read_parquet(?)
            WHERE sdss_id IN ?
        """, [parquet_file, sdssids]).fetchdf().to_dict('records')

    results = Table(names=('SDSS_ID','CATALOGID_V0','CATALOGID_V0P5','CATALOGID_V1'),
                    dtype=(int, int, int, int))
    for t in tp:
        for col in t:
            try:
                if np.isnan(t[col]):
                    t[col] = -999
            except:
                if t[col] is None:
                    t[col] = -999
        results.add_row((int(t['sdss_id']),int(t['catalogid21']),int(t['catalogid25']),int(t['catalogid31'])))
    if len(results) > 0:
        search_table = join(search_table, results, keys='SDSS_ID', join_type='left')
    else:
        search_table['CATALOGID_V0']   = -999
        search_table['CATALOGID_V0P5'] = -999
        search_table['CATALOGID_V1']   = -999
    for col in ['CATALOGID_V0','CATALOGID_V0P5','CATALOGID_V1']:
        try:
            search_table[col] = search_table[col].filled(-999)
        except:
            pass
    return(search_table)

def get_targetflags(search_table, data, db=True, release='sdsswork', V_TARG='*', no_remote=False):
    warnings.filterwarnings("default", module="sdss_semaphore")

    sdssids = np.unique(search_table['SDSS_ID'].data).tolist()

    try:
        sem_opts = dict(verbose = True, sdssc2bv = os.getenv('SDSSC2BV',None))
        TargetingFlags(**sem_opts)
    except:
        sem_opts = {}
    if db is True:

        splog.info('Getting Targeting flags')

        try:
            tp = SDSS_ID_flat.select(SDSS_ID_flat.sdss_id,CartonToTarget.carton_pk)\
                             .join(Target, on=(SDSS_ID_flat.catalogid == Target.catalogid))\
                             .join(CartonToTarget, on=(Target.pk == CartonToTarget.target_pk))\
                             .where(SDSS_ID_flat.sdss_id.in_(sdssids)).tuples()
        except:
            splog._log.exception('Error getting Targeting Flags, trying again....')
            time.sleep(60)
            tp = SDSS_ID_flat.select(SDSS_ID_flat.sdss_id,CartonToTarget.carton_pk)\
                             .join(Target, on=(SDSS_ID_flat.catalogid == Target.catalogid))\
                             .join(CartonToTarget, on=(Target.pk == CartonToTarget.target_pk))\
                             .where(SDSS_ID_flat.sdss_id.in_(sdssids)).tuples()
    else:
        splog.info('Getting Target flags from SDSS-V MOS Target Product')
        sdss_id_flat_path = get_Catalog('mos_target_sdss_id_flat', release=release, V_TARG=V_TARG, no_remote=no_remote)
        target_path = get_Catalog('mos_target_target', release=release, V_TARG=V_TARG, no_remote=no_remote)
        carton_to_target_path = get_Catalog('mos_target_carton_to_target', release=release, V_TARG=V_TARG, no_remote=no_remote)

        tp = duckdb.execute("""
            SELECT
                s.sdss_id,
                c.carton_pk
            FROM read_parquet(?) s
            JOIN read_parquet(?) t
                ON s.catalogid = t.catalogid
            JOIN read_parquet(?) c
                ON t.target_pk = c.target_pk
            WHERE s.sdss_id IN ?
        """, [sdss_id_flat_path, target_path, carton_to_target_path, sdssids]).fetchall()


    if len(tp) == 0:
        splog.info('No Matching Targets')
        try:
            SDSSC2BV = str(TargetingFlags(**sem_opts).version)
        except:
            SDSSC2BV = '1'
        search_table['SDSS5_TARGET_FLAGS'] = Column(name = 'SDSS5_TARGET_FLAGS',
                                                    dtype = 'uint8', shape=(1,),
                                                    length=len(search_table)).astype(object)
        search_table['SDSSC2BV'] = Column(SDSSC2BV, name = 'SDSSC2BV', dtype = object)

        data['SDSS5_TARGET_FLAGS'] = Column(name = 'SDSS5_TARGET_FLAGS',
                                            dtype = "uint8",shape=(1,),
                                            length=len(data)).astype(object)#, shape = (,F))
        data['SDSSC2BV'] = Column(name = 'SDSSC2BV', dtype = object)
        return(search_table, data)

    manual_counts = {}
    flags_dict = {}
    pks_dict = {}
    for sdss_id, carton_pk in tp:
        try:
            flags_dict[sdss_id]
            pks_dict[sdss_id]
        except KeyError:
            flags_dict[sdss_id] = TargetingFlags(**sem_opts)
            pks_dict[sdss_id] = []

        try:
            pks_dict[sdss_id].append(carton_pk)
            flags_dict[sdss_id].set_bit_by_carton_pk(0, carton_pk) # 0 since this is the only object
            manual_counts.setdefault(carton_pk, set())
            manual_counts[carton_pk].add(sdss_id)
        except Exception as e:
            pass
    # Now we will create two columns:
    # - one for all our source identifiers
    # - one for all our targeting flags

    sdss_ids = list(flags_dict.keys())
    flags =TargetingFlags(list(flags_dict.values()),**sem_opts)
    
    # A sanity check.
    for carton_pk, count in flags.count_by_attribute("carton_pk", skip_empty=True).items():
        assert count == len(manual_counts[carton_pk])
                
    N, F = flags.array.shape
    results = Table() #names=('icatalogid','SDSS5_TARGET_FLAGS'), dtype = (int,"{F}B"))
    results.add_column(sdss_ids, name = 'SDSS_ID')
    results.add_column(flags.array, name = 'SDSS5_TARGET_FLAGS')
    
    try:
        SDSSC2BV = str(TargetingFlags(**sem_opts).version)
    except:
        SDSSC2BV = '1'
    
    results['SDSSC2BV'] = Column(SDSSC2BV, name = 'SDSSC2BV', dtype = object)
    
    if data is not None:
        data['SDSS5_TARGET_FLAGS'] = Column(name = 'SDSS5_TARGET_FLAGS', dtype = f"{F}B")
        data['SDSSC2BV'] = Column(name = 'SDSSC2BV', dtype = object)
    search_table = join(search_table, results, keys='SDSS_ID',join_type='left')
    STF = search_table['SDSS5_TARGET_FLAGS']
    sdssids = search_table['SDSS_ID'].data
    STF[np.where(sdssids == -999)[0]] = np.zeros(F, dtype='uint8')
    search_table['SDSS5_TARGET_FLAGS'] = STF
    return(search_table, data)

def get_CartonInfo(search_table, db=True, release='sdsswork', V_TARG='*', no_remote=False):
    carton_to_target_pk = search_table['carton_to_target_pk'].data.tolist()
    if db is True:
        tp = CartonToTarget.select(CartonToTarget.pk,Carton.program, Carton.carton, Version.plan, Mapper.label).join(Carton).join(Version).\
                            switch(Carton).join(Mapper).where(CartonToTarget.pk.in_(carton_to_target_pk)).dicts()
        
    else:
        splog.warning('Getting Carton and Mapper info from SDSS-V MOS Targeting Product')

        carton_to_target_path = get_Catalog('mos_target_carton_to_target', release=release, V_TARG=V_TARG, no_remote=no_remote)
        carton_path = get_Catalog('mos_target_carton', release=release, V_TARG=V_TARG, no_remote=no_remote)
        version_path = get_Catalog('mos_target_targetdb_version', release=release, V_TARG=V_TARG, no_remote=no_remote)
        mapper_path = get_Catalog('mos_target_mapper', release=release, V_TARG=V_TARG, no_remote=no_remote)

        tp = duckdb.execute("""
                SELECT ctt.carton_to_target_pk as pk, c.program, c.carton, v.plan, m.label
                FROM read_parquet(?) ctt
                JOIN read_parquet(?) c
                    ON ctt.carton_pk = c.carton_pk
                JOIN read_parquet(?) v
                    ON c.version_pk = v.pk
                JOIN read_parquet(?) m
                    ON c.mapper_pk = m.pk
                WHERE ctt.carton_to_target_pk IN ?
            """, [
                carton_to_target_path,
                carton_path,
                version_path,
                mapper_path,
                carton_to_target_pk
            ]).fetchdf().to_dict('records')

        
    results = Table(names = ('carton_to_target_pk', 'program_db', 'carton', 'CatVersion', 'mapper'),
                    dtype = (int, object, object, object, object))
    for t in tp:
        for key in t.keys():
            if t[key] is None:
                if key in ['program','carton','plan','label']:
                    t[key] = ''
        results.add_row((t['pk'],t['program'],t['carton'],t['plan'],t['label']))

    for c2t in np.array(carton_to_target_pk):
        if c2t in results['carton_to_target_pk'].data:
            carton_to_target_pk.remove(c2t)

    if db:
        tp = CartonToTarget.select(CartonToTarget.pk,Carton.program, Carton.carton, Version.plan).join(Carton).join(Version).\
                            switch(Carton).where(CartonToTarget.pk.in_(carton_to_target_pk)).dicts()
        
    else:
        tp = duckdb.execute("""
                SELECT ctt.carton_to_target_pk as pk, c.program, c.carton, v.plan
                FROM read_parquet(?) ctt
                JOIN read_parquet(?) c
                    ON ctt.carton_pk = c.carton_pk
                JOIN read_parquet(?) v
                    ON c.version_pk = v.pk
                WHERE ctt.carton_to_target_pk IN ?
            """, [
                carton_to_target_path,
                carton_path,
                version_path,
                carton_to_target_pk
            ]).fetchdf().to_dict('records')
    for t in tp:
        for key in t.keys():
            if t[key] is None:
                if key in ['program','carton','plan','label']:
                    t[key] = ''
        results.add_row((t['pk'],t['program'],t['carton'],t['plan'],''))
    if len(results) > 0:
        search_table = join(search_table,results,keys='carton_to_target_pk', join_type='left')

    return(search_table)

