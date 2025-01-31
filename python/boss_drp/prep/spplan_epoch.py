#!/usr/bin/env python3
import boss_drp
from boss_drp.utils.splog import splog
from boss_drp.prep.spplan import write_plan
from boss_drp.prep.GetconfSummary import get_confSummary
from boss_drp.field import field_to_string, Field
from boss_drp.utils import load_env, get_dirs


try:
    from sdssdb.peewee.sdss5db import opsdb, targetdb
    opsdb.database.set_profile(load_env('DATABASE_PROFILE', default='pipelines'))
    import sdssdb
    SDSSDBVersion=sdssdb.__version__
except:
    if load_env('DATABASE_PROFILE', default='pipelines').lower() in ['pipelines','operations']:
        splog.log('ERROR: No SDSSDB access')
        exit()
    else:
        splog.log('WARNING: No SDSSDB access')
    SDSSDBVersion='N/A'

from sdss_access.path import Path
from sdss_access import Access

import numpy as np
import argparse
from os import getenv, environ, makedirs
from glob import glob
import os.path as ptt
from pydl.pydlutils import yanny
from astropy.table import Table, vstack, Column
from astropy.io import fits
from astropy.time import Time
from sys import argv
from collections import OrderedDict
from subprocess import Popen, PIPE, getoutput
import time



idlspec2dVersion = boss_drp.__version__
#idlutilsVersion = getoutput("idlutils_version")


def expsByEpoch(fieldPk):
    """
        Returns dictionary of epoch indices containing lists of exposure_pks
    """

    Design = targetdb.Design
    Conf = opsdb.Configuration
    Exp = opsdb.Exposure
    Field = targetdb.Field
    Cadence = targetdb.Cadence
    d2f = targetdb.DesignToField
    DesignToStatus = opsdb.DesignToStatus
    CompletionStatus = opsdb.CompletionStatus
    Version = targetdb.Version
    Obs = targetdb.Observatory


    exposures = Design.select(Design.design_id, d2f.exposure.alias("exp_num"),Exp.pk)\
                      .join(Conf).join(Exp).switch(Design).join(d2f)\
                      .where(d2f.field_pk == fieldPk)\
                      .switch(Design).join(DesignToStatus)\
                      .dicts()

    field = Field.select(Cadence.nexp,Field.field_id)\
                 .join(Cadence).where(Field.pk == fieldPk).dicts()[0]
    expCount = [np.sum(field["nexp"][:i+1]) for i in range(len(field["nexp"]))]

    previous_idx = 0

    exps_epoch = dict()
    for epoch, idx in enumerate(expCount):
        exps_epoch[epoch] = [e["pk"] for e in exposures if
                             (e["exp_num"] < idx and e["exp_num"] >= previous_idx)]
        des_in_epoch = [e["design_id"] for e in exposures if
                        (e["exp_num"] < idx and e["exp_num"] >= previous_idx)]
        previous_idx = idx

    return exps_epoch

def getDesignStatus(allexps):
 
    designs = np.unique(allexps['design'].data)
    designs = designs.tolist()
    
    Design = targetdb.Design
    Conf = opsdb.Configuration
    Exp = opsdb.Exposure
    Field = targetdb.Field
    Cadence = targetdb.Cadence
    d2f = targetdb.DesignToField
    DesignToStatus = opsdb.DesignToStatus
    CompletionStatus = opsdb.CompletionStatus
    Version = targetdb.Version
    Obs = targetdb.Observatory
    
    tstatus = Design.select(Design.design_id, DesignToStatus.completion_status_pk, DesignToStatus.manual, DesignToStatus.mjd)\
                    .join(DesignToStatus).where(Design.design_id.in_(designs)).dicts()
    
    allexps['manual'] = 0
    allexps['Status'] = 1
    allexps.add_column(Column(np.NaN, dtype=object, name = 'designmjd'))
    
    manual = allexps['manual'].data
    status = allexps['Status'].data
    designmjd = allexps['designmjd'].data

    for t in tstatus:
        idx = np.where(allexps['design'].data.astype(int) == t['design_id'])[0]
        manual[idx] = 1
        status[idx] = t['status']
        designmjd[idx] = t['mjd']
    allexps['manual'] = manual
    allexps['Status'] = status
    allexps['designmjd'] = designmjd
    return(allexps)



def get_FieldMeta(allexps, obs, plates=False, release = 'sdsswork'):
    """ 
        Obtains the field metadata associated with an exposure
    """
    
    Design = targetdb.Design
    Conf = opsdb.Configuration
    Exp = opsdb.Exposure
    Field = targetdb.Field
    Cadence = targetdb.Cadence
    d2f = targetdb.DesignToField
    DesignToStatus = opsdb.DesignToStatus
    CompletionStatus = opsdb.CompletionStatus
    Version = targetdb.Version
    Obs = targetdb.Observatory
    
    allexps.add_column(obs,name ='obs')
    rs_plans = np.full(len(allexps),'', dtype=object)

    if not plates:
        allexps.add_column(-1,name ='design')
        designid = allexps['design'].data
        
        splog.log('Getting DesignIDs')
        configID = allexps['confid'].data.astype(int).astype(str).tolist()
        conf = Conf.select(Conf.configuration_id, Conf.design_id).where(Conf.configuration_id.in_(configID)).dicts()
        for cfid in conf:
            idx = np.where(allexps['confid'].data.astype(int) == int(cfid['configuration_id']))[0]
            designid[idx] = cfid['design']
            
            comfSumm = get_confSummary(cfid['configuration_id'], obs=obs, release=release)
            
            rs_plan = comfSumm.meta['robostrategy_run']
            rs_plans[idx] = rs_plan

        allexps['design'] = designid
        rs_plans = rs_plans.astype(str)
        allexps.add_column(rs_plans,name ='rs_plan')
        
        output = None
        splog.log('Getting Field Meta Data')
        for rs_plan in np.unique(rs_plans):
            idx = np.where(rs_plans == rs_plan)[0]
            sub_desid = designid[idx]
            field = Field.select(Field.pk, Design.design_id, Field.cadence.label.alias('fcadence'),Field.cadence.max_length.alias('max_length'),
                                 Field.observatory.label.alias('obs'), Field.version.plan.alias('rs_plan'))\
                         .join(d2f).join(Design).where(Design.design_id.in_(np.unique(designid).tolist())).switch(Field).join(Version)\
                         .where(Version.plan==rs_plan).switch(Field).join(Cadence).switch(Field).join(Obs).dicts()
            for fid in field:
                idx2 = np.where((rs_plans == rs_plan) & (designid == fid['design_id']))[0]
                if len(idx2) == 0: continue
                output1 = Table(allexps[idx2].copy())

                output1.add_column(fid['pk'], name='field_pk')
                output1['rs_plan']=fid['rs_plan']
                output1['obs']=fid['obs']
                output1.add_column(Column(fid['fcadence'], dtype=object, name='field_cadence'))
                output1.add_column([fid['max_length']], name='max_length')

                if output is None:
                    output = output1
                else:
                    
                    if output1['max_length'].shape[1] > output['max_length'].shape[1]:
                        pad = np.zeros([output['max_length'].shape[0], output1['max_length'].shape[1] - output['max_length'].shape[1]])
                        output['max_length'] = np.hstack([output['max_length'], pad])
                    if output1['max_length'].shape[1] < output['max_length'].shape[1]:
                        pad = np.zeros([output1['max_length'].shape[0], output['max_length'].shape[1] - output1['max_length'].shape[1]])
                        output1['max_length'] = np.hstack([output1['max_length'], pad])
                    output = vstack([output,output1])
        
        
        try:
            output['field_cadence']=output['field_cadence'].astype(str)
        except:
            print(output)
            splog.log("No Matching DB Field Meta Data")
            return(None)
            #output['field_cadence']=output['field_cadence'].astype(str)
        output['design']=output['design'].astype(str)

        allexps = output
    else:
        allexps['rs_plan'] = 'plates'
        allexps['field_pk'] = -999
        allexps['field_cadence'] = 'plates'
        allexps['max_length'] = 1000
        allexps['obs'] = 'APO'
        allexps['design'] = ''

    return(allexps)


def getExpID(expPK):
    """
        Gets exposure number from exposure primary key
    """
    Exp = opsdb.Exposure

    exp = Exp.select(Exp.exposure_no, Exp.pk, Exp.start_time).where(Exp.pk == expPK).dicts()[0]
    mjd = Time(str(exp['start_time']), format='iso').to_value('mjd')
    return(exp['exposure_no'], exp['start_time'], mjd)


def write_spPlancomb(allexps, fpk, field, topdir=None, run2d=None, clobber=False, plates=False, fpk_flag=None,
                    abandoned=False, min_epoch_len=0, daily=False):
    """
        Writes spPlancomb file for each epoch in a Table of exposures

    """
    allexps['sdrname'] = allexps['name'].data
    shape = (allexps['name'].shape[1],)
    allexps.remove_column('name')
    allexps.add_column(Column([['']*shape[0]], dtype='S30', name = 'name'))
    for i, row in enumerate(allexps):
        names = np.asarray(allexps[i]['sdrname'].data).astype(object)
        for j, name in enumerate(names):
            names[j] = name.decode().replace('sdR','spFrame').replace('.gz','').replace('.fit','.fits').encode('utf-8')
        allexps[i]['name'] = names.tolist()
    allexps.remove_column('sdrname')
    for epoch in list(dict.fromkeys(allexps['epoch_combine'].data)):
        epochexps = allexps[np.where(allexps['epoch_combine'] == epoch)[0]]
        EpochID=epoch
        if np.max(epochexps['mjd'].data)-np.min(epochexps['mjd'].data) < min_epoch_len:
            splog.log('Skipping epoch with len '+str(np.max(epochexps['mjd'].data)-np.min(epochexps['mjd'].data)))
            continue
        else:
            N_MJD = np.max(epochexps['mjd'].data)-np.min(epochexps['mjd'].data)+1
        if len(epochexps) == 0:
            continue
        mjd = max(epochexps['mjd'])
        if (EpochID == mjd) and (daily is False):
            EpochID = -1
            if not abandoned:
                splog.log('Skipping Abandoned Epoch')
                continue
        epochexps['epoch_combine'] = mjd
        
        epochexps = epochexps[['confid','fieldid','mjd','mapname','flavor','exptime','name','epoch_combine',
                               'planfile2d','field_pk','field_cadence','rs_plan','obs']]
        fc = Field(topdir, run2d, field, epoch = not daily)
        outdir = fc.dir()
        makedirs(outdir, exist_ok=True)
        if not daily:
            if fpk_flag is not None:
                planfile = ptt.join(outdir, 'spPlancombepoch-'+field_to_string(field)+'_'+str(fpk_flag)+'-'+str(mjd)+'.par')
            else:
                planfile = ptt.join(outdir, 'spPlancombepoch-'+field_to_string(field)+'-'+str(mjd)+'.par')
        else:
            if fpk_flag is not None:
                planfile = ptt.join(outdir, 'spPlancomb-'+field_to_string(field)+'_'+str(fpk_flag)+'-'+str(mjd)+'.par')
            else:
                planfile = ptt.join(outdir, 'spPlancomb-'+field_to_string(field)+'-'+str(mjd)+'.par')

        plan2ds =       ' '.join(np.unique(np.array(epochexps["planfile2d"])).tolist())
        field_pk =      ' '.join(np.unique(np.array(epochexps["field_pk"]).astype(str)).tolist())
        field_cadence = ' '.join(np.unique(np.array(epochexps["field_cadence"])).tolist())
        rs_plan =       ' '.join(np.unique(np.array(epochexps["rs_plan"])).tolist())
        obs =           ' '.join(np.unique(np.array(epochexps["obs"])).tolist())
        for col in ['planfile2d', 'rs_plan', 'field_cadence', 'obs', 'design', 'max_length', 'epoch_length',
                    'expid','field_pk', 'manual', 'Status','designmjd','start_time','start_mjd']:
            if col in epochexps.colnames:
                epochexps.remove_column(col)

        meta=OrderedDict({'fieldid':          str(field)             +"   # Field number",
                                    'MJD':              str(mjd)               +"   # Modified Julian Date for most recent observation",
                                    'OBS':              obs                    +"   # Observatory",
                                    'DITHER':           'F'                    +"   # Is the Field Dithered (T: True, F: False)",
                                    'FieldCadence':     field_cadence          +"   # Field Cadence",
                                    'FieldPk':          field_pk               +"   # Field Primary Key",
                                    'EpochID':          str(EpochID)           +"   # Completed Plan Epoch ID (or MJD of incomplete epoch or daily coadd)",
                                    'N_MJD':            str(N_MJD)             +"   # Number of MJDs included in the Epoch",
                                    'planfile2d':       plan2ds                +"   # Plan file(s) for Daily 2D spectral reductions",
                                    'planfilecomb':     ptt.basename(planfile) +"   # Plan file for Combine (this file)",
                                    'idlspec2dVersion': idlspec2dVersion       +"   # Version of idlspec2d when building plan file",
                                    'sdssdb_Version':   SDSSDBVersion          +"   # Version of sdssdb when building this plan file",
                                    'RS_Version':       rs_plan                +"   # Robostrategy Version for this field",
                                    })
        write_plan(planfile, epochexps, meta=meta, clobber=clobber, override_manual=False)
    
    allexps_out = allexps.copy()

    for col in ['max_length','designmjd','start_time']:
        if col in allexps_out.colnames:
                allexps_out.remove_column(col)
    if True: #not daily:
        fc = Field(topdir, run2d, field)
        if fpk_flag is not None:
            expfile = ptt.join(fc.dir(),'SciExp-'+field_to_string(field)+'_'+str(fpk_flag)+'.par')
        else:
            expfile = ptt.join(fc.dir(),'SciExp-'+field_to_string(field)+'.par')

        if ptt.exists(expfile):
#            if clobber is False:
#                splog.info('WARNING: Will not over-write SciExp file: ' + ptt.basename(expfile))
#                return
#            else:
            splog.info('WARNING: Over-writing SciExp file: ' + ptt.basename(expfile))
        else: splog.info('Writing SciExp file '+ ptt.basename(expfile))
        allexps_out.convert_unicode_to_bytestring()
        yanny.write_table_yanny(allexps_out, expfile,tablename='SPEXP', overwrite=True)
    return


    
def get_exp_spx(topdir, run2d, field, plates=False, lco=False, release = 'sdsswork'):
    """
        Builds a Table of exposures from the spPlan2d files in a field directory, and calls get_FieldMeta to add the field metadata
    """
    splog.log('Building table of exposures from spPlan2d files')
    allexp = []
    fc = Field(topdir, run2d, field)
    for plan in glob(ptt.join(fc.dir(), 'spPlan2d*.par')):
        yplan = yanny.read_table_yanny(plan, 'SPEXP')
        if not plates:
            obs = yplan.meta['OBS'].lower()
            if lco:
                if obs == 'apo':
                    splog.log(f'Skipping {obs} Field')
                    return None
            else:
                if obs == 'lco':
                    splog.log(f'Skipping {obs} Field')
                    return None
        else:
            obs = 'apo'
            if obs == 'lco':
                splog.log(f'Skipping {obs} Field')
                return None
        yplan['planfile2d']=yplan.meta['planfile2d']
        yplan.meta = OrderedDict()
        yplan = yplan[np.where(yplan['flavor']=='science')[0]]
        allexp.append(yplan)
    allexp = vstack(allexp)
    allexp = get_FieldMeta(allexp, obs, plates=plates, release=release)
    if allexp is not None:
        names = allexp['name']
        expids = [int(ptt.splitext(x)[0].split('-')[-1]) for x in np.asarray(names[:,0],dtype=str).tolist()]
        allexp['expid'] = expids
    return(allexp)



def fps_field_epoch(field, topdir=None, run2d=None, clobber=False, lco = False, abandoned=False,
                    started=False, min_epoch_len=0, release = 'sdsswork', mjd=None,
                    mjdstart = None, mjdend = None):
    """
        Separates a Table of fps exposures into epoches
    """
    allexps = get_exp_spx(topdir, run2d, field, lco=lco, release=release)
    if allexps is None:
        return
    allexps = filter_mjd(allexps, mjd=mjd, mjdstart=mjdstart, mjdend=mjdend)
    
    if allexps is None:
        return
    if len(allexps) == 0:
        return
    splog.log('Getting Design Status')
    allexps = getDesignStatus(allexps)
    
    allexps.sort(['expid'])

    if not started:
        allexps = allexps[allexps['Status'] == 3]
    if (allexps is None) or (len(allexps) == 0):
        return
    fpks = np.flip(np.sort(np.unique(allexps['field_pk'].data)))
    fpk = max(fpks)
    splog.log('Getting Exposures in Epochs')
    epochs = expsByEpoch(fpk)
    allexps['epoch_combine'] = -1
    allexps['epoch_length'] = 0.0
    allexps.add_column(Column('',dtype=object,name='start_time'))
    allexps.add_column(0.0,name='start_mjd')
    for epoch in epochs.keys():
        # get all exposures included in epoch by designid defintion
        for expPk in epochs[epoch]:
            expid, start_time, start_mjd = getExpID(expPk)
            allexps['epoch_combine'][np.where(allexps['expid'] == expid)[0]] = epoch
            allexps['start_time'][np.where(allexps['expid'] == expid)[0]] = start_time
            allexps['start_mjd'][np.where(allexps['expid'] == expid)[0]] = start_mjd
    ml = allexps[np.where( allexps['epoch_combine'] != -1)[0]]['max_length'].data
    ec = allexps[np.where( allexps['epoch_combine'] != -1)[0]]['epoch_combine'].data
    el = allexps['epoch_length']

    el[np.where( allexps['epoch_combine'] != -1)[0]] =  [x[ec[i]] for i, x in enumerate(ml)]
    el[np.where(el.data == 0.0)[0]] = 0.5

    allexps[np.where(allexps['epoch_combine'] == -1)[0]]['epoch_combine'] = allexps[np.where(allexps['epoch_combine'] == -1)[0]]['mjd']
    
    for epoch in np.unique(allexps['epoch_combine'].data):
        inepoch = np.where(allexps['epoch_combine'] == epoch)[0]
        epochexps = allexps[inepoch]
        max_length = max(epochexps['epoch_length'])
        mjd = max(epochexps['mjd'])
        out = np.where(epochexps['mjd'] < mjd-max_length)[0]
        ec = allexps['epoch_combine']
        mj = allexps['mjd']
        ec[inepoch[out]] = -1
    
    
    mjds = allexps[np.where(allexps['epoch_combine'] == -1)[0]]['mjd'].data
    splog.log('Building partial epochs')
    for i, mjd in enumerate(np.flip(np.unique(np.sort(allexps['mjd'])))):
        idx = np.where(allexps['mjd'].data == mjd)[0]
        max_length = (allexps['epoch_length'].data)[idx[-1]]
        idx_set = np.where((float(mjd) - np.asarray(allexps['start_mjd'].data).astype(float) <= max_length) & (allexps['epoch_combine'].data == -1))[0]
        ec[idx_set] = mjd

    write_spPlancomb(allexps, fpk, field, topdir=topdir, run2d=run2d, clobber=clobber, abandoned=abandoned, min_epoch_len=min_epoch_len)
    return
    
def plate_field_epoch(field, topdir=None, run2d=None, clobber=False,mjd=None, mjdstart=None, mjdend=None,
                     abandoned=False, min_epoch_len=0, daily=False, release = 'sdsswork'):
    """
        Separates a Table of plate exposures into epoches
    """
    obs = 'APO'
    FieldCadence = 'plate'
    pk = -999
    plan = 'plate'
    fpk = Table({'obs':[obs], 'FieldCadence': [FieldCadence],'pk':[pk], 'plan':[plan]})
    allexps = get_exp_spx(topdir, run2d, field, plates=True, release=release)
    allexps = filter_mjd(allexps, mjd=mjd, mjdstart=mjdstart, mjdend=mjdend)
    allexps['epoch_combine'] = -1
    
    
    epochs = allexps['epoch_combine']

    dmjd = 1000
    if int(field) in [15000,15001,15002,15038,15070,15071,15171,15172,15173,15252,15253]:
        dmjd = 3 # RM Plates

    if daily is True:
        allexps['epoch_combine'] = allexps['mjd']
        allexps['epoch_length'] = 0
        allexps['max_length'] = dmjd
        write_spPlancomb(allexps, fpk, field, topdir=topdir, run2d=run2d, clobber=clobber, plates=True, daily=True)
        
    else:
        i = 0
        allexps['max_length'] = dmjd
        for map in np.unique(allexps['mapname'].data):
            idx = np.where(allexps['mapname'].data == map)[0]
            sallexps = allexps[idx]
            for j, mjd in enumerate(np.flip(np.sort(np.unique(sallexps['mjd'].data)))):
                mallexps = sallexps[np.where(np.asarray(sallexps['mjd'].data).astype(int) == int(mjd))[0]]
                if (sum(mallexps['epoch_combine'].data == -1) == 0): continue
                idx1 = np.where((int(mjd) - np.asarray(sallexps['mjd'].data).astype(int) < dmjd) & (epochs[idx] == -1))[0]
                sallexps['epoch_combine'][idx1] = i
                epochs[idx[idx1]] = i
                if len(idx1) > 0: i = i+1
            epochs_raw = np.asarray(sallexps['epoch_combine'].copy().data)
            eps = np.sort(np.unique(epochs_raw))
            for j, ep in enumerate(np.flip(np.sort(np.unique(epochs_raw)))):
                epochs[idx[np.where(epochs_raw == ep)[0]]] = eps[j]
        epochs[idx[np.where(epochs[idx] == -1)[0]]] = sallexps[np.where(epochs[idx] == -1)[0]]['mjd'].data
        write_spPlancomb(allexps, fpk, field, topdir=topdir, run2d=run2d, clobber=clobber, plates=True,
                         abandoned=abandoned,  min_epoch_len=min_epoch_len)
    return

def fps_field_daily(field, topdir=None, run2d=None, clobber=False,
                    lco=False,mjd=None, mjdstart=None, mjdend=None,
                    release = 'sdsswork'):
    """
        Sperate a Table of fps exposures into daily Tables
    """
    allexps = get_exp_spx(topdir, run2d, field, release = release)
    allexps['epoch_combine'] = allexps['mjd']
    allexps['epoch_length'] = 0
    fpks = np.unique(allexps['field_pk'].data)
    fpk = ' '.join(fpks.astype(str))
    allexps = filter_mjd(allexps, mjd=mjd, mjdstart=mjdstart, mjdend=mjdend)
    allexps.sort(['expid'])
    write_spPlancomb(allexps, fpk, field, topdir=topdir, run2d=run2d, clobber=clobber, daily=True)
    return

def filter_mjd(allexps, mjd=None, mjdstart=None, mjdend=None):
    
    if mjd is not None: allexps = allexps[np.where(allexps['mjd'].data == int(mjd))[0]]
    if mjdstart is not None: allexps = allexps[np.where(allexps['mjd'].data >= int(mjdstart))[0]]
    if mjdend is not None: allexps = allexps[np.where(allexps['mjd'].data <= int(mjdend))[0]]
    return(allexps)


def spplan_epoch(topdir=None, run2d=None, fieldid=None, fieldstart=None, fieldend=None, abandoned=False,
                clobber=False, daily=False, lco=False, mjd=None, mjdstart=None, mjdend=None, started=False,
                min_epoch_len=0, release = 'sdsswork', **kwds):
    """
        Generates list of fields and calls the correct field epoch function for each field
    """
    fc = Field(topdir, run2d, '*')
    fieldlist = get_dirs(ptt.dirname(fc.dir()), field = True,
                         match=fieldid, start=fieldstart, end=fieldend)
    splog.info('Number of Field Directories = '+ str(len(fieldlist))+'\n')
    for i, field in enumerate(fieldlist):
        if i > 0:
            splog.info('------------------------------------')
        splog.info('Field Directory: '+str(field))
            
        if int(field) < 16000:
            if lco is True:
                continue
            plate_field_epoch(field, topdir=topdir, run2d=run2d, clobber=clobber, abandoned=abandoned,
                              min_epoch_len=min_epoch_len, daily=daily, release=release, mjd=mjd,
                              mjdstart=mjdstart, mjdend=mjdend)
        else:
            if daily is False:
                if SDSSDBVersion == 'N/A':
                    splog.info(f'Skipping {field} due to no SDSSDB access')
                fps_field_epoch(field, topdir=topdir, run2d=run2d, clobber=clobber, lco=lco, abandoned=abandoned,
                                started=started, min_epoch_len=min_epoch_len, release=release, mjd=mjd,
                                mjdstart=mjdstart, mjdend=mjdend)
            else:
                fps_field_daily(field, topdir=topdir, run2d=run2d, clobber=clobber, release=release, mjd=mjd,
                                mjdstart=mjdstart, mjdend=mjdend)
    
def spplancombin(topdir=None, run2d=None, run1d=None, mjd=None, mjdstart=None, mjdend=None,
                fieldid=None, fieldstart=None, fieldend=None, clobber=False, abandoned=False,
                minexp=1, daily=False, lco=False, apo=False, logfile=None, started=False,
                min_epoch_len = 0, release = 'sdsswork', **extra_kwds):
    """
        Builds the spplancomb plan files
    """

    if logfile is not None:
        splog.open(logfile=logfile, logprint=False)
        splog.log('Log file '+logfile+' opened '+ time.ctime())

    if topdir is None: topdir = getenv('BOSS_SPECTRO_REDUX')
    splog.info('Setting TOPDIR='+ topdir)
    if run2d is None: run2d = getenv('RUN2D')
    splog.info('Setting RUN2D='+ run2d)
    if run1d is None: run1d = getenv('RUN1D')
    splog.info('Setting RUN1D='+ run1d)


    if lco is True:
        environ['OBSERVATORY'] = 'LCO'
        opsdb.database.connect()  # This will recreate the base model class and reload all the model classes
    else:
        environ['OBSERVATORY'] = 'APO'
        opsdb.database.connect()  # This will recreate the base model class and reload all the model classes

    spplan_epoch(topdir=topdir, run2d=run2d, clobber=clobber, daily=daily,
                mjd=mjd, mjdstart=mjdstart, mjdend=mjdend, lco=lco, abandoned=abandoned,
                fieldid=fieldid, fieldstart=fieldstart, fieldend=fieldend, started=started,
                min_epoch_len=min_epoch_len, release = release)
                
    if logfile is not None:
        splog.log('Successful completion of spplan_epoch at '+ time.ctime())
        splog.close()


