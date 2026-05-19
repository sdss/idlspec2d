#!/usr/bin/env python3
from boss_drp.prep.readfibermaps.db_tools import get_targetflags, get_AltCatids
from boss_drp.utils import load_env
from boss_drp.utils.splog import splog
from boss_drp.summary import summary_names, Summary_names
from boss_drp.utils.parquet.write import write_parquet
from boss_drp.utils.parquet.stream import stream_writer
from boss_drp.utils.parquet.schema import Schema
spAll_schema = Schema()
spLite_schema = Schema()


import os
os.environ['SDSSC2BV'] = "3"

from sdss_semaphore.targeting import TargetingFlags

try:
    from sdssdb.peewee.sdss5db.targetdb import database
    test = database.set_profile(load_env('DATABASE_PROFILE', default='pipelines'))
    # from sdssdb.peewee.sdss5db.targetdb import CartonToTarget, Carton, Version, Mapper, Target
    nodb=False
except:
    if load_env('DATABASE_PROFILE', default='pipelines').lower() in ['pipelines','operations']:
        splog.log('ERROR: No SDSSDB access')
        exit()
    else:
        splog.log('ERROR: No SDSSDB access')
    nodb= True

from astropy.table import Table# , Column, join
#from astropy.io import fits
from os import getenv, remove
import os.path as ptt
#from glob import glob
import numpy as np
import time
import shutil
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import pyarrow.compute as pc
import pyarrow as pa

import pandas as pd
from pathlib import Path


def get_Targeting_file(run2d, boss_spectro_redux=getenv('BOSS_SPECTRO_REDUX')):
    return str(Path(summary_names.spAllfile).parent / f'spTargeting-{run2d}.parquet')


class TargetFlagsUpdater:
    def __init__(self, search_parquet=None, search_table=None,
                 key="SDSS_ID", flag_col="SDSS5_TARGET_FLAGS", meta_col="SDSSC2BV"):
        self.key = key
        self.flag_col = flag_col
        self.meta_col = meta_col
        self.sdssc2bv = None

        if search_parquet is not None:
            ids, flags, self.sdssc2bv = self._load_from_parquet(search_parquet)

        elif search_table is not None:
            ids, flags, self.sdssc2bv = self._load_from_table(search_table)

        else:
            raise ValueError("Provide either search_parquet or search_table")

        # Sort once so lookup can use searchsorted
        order = np.argsort(ids)
        self.ids = np.asarray(ids)[order]
        self.flags = np.asarray(flags, dtype=np.uint8)[order]

        if self.flags.ndim == 1:
            self.flags = self.flags[:, None]

    def _load_from_parquet(self, path):
        # Read only the needed columns
        t = pq.read_table(path, columns=[self.key, self.flag_col])

        ids = t[self.key].to_numpy(zero_copy_only=False)

        # Works for list / fixed_size_list columns
        flags = np.asarray(t[self.flag_col].to_pylist(), dtype=np.uint8)

        # Read metadata if present
        sdssc2bv = None
        pf = pq.ParquetFile(path)
        md = pf.metadata.metadata if pf.metadata is not None else None
        if md is not None and self.meta_col.encode() in md:
            sdssc2bv = md[self.meta_col.encode()].decode()

        return ids, flags, sdssc2bv

    def _load_from_table(self, table):
        if isinstance(table, Table):
            ids = np.asarray(table[self.key])
            col = table[self.flag_col]
            data = col.data if hasattr(col, "data") else col
            if np.ma.isMaskedArray(data):
                flags = np.asarray(np.ma.filled(data, 0), dtype=np.uint8)
            else:
                flags = np.asarray(data, dtype=np.uint8)
            sdssc2bv = table.meta.get(self.meta_col, None)

        elif isinstance(table, pd.DataFrame):
            ids = table[self.key].to_numpy()
            flags = np.asarray(table[self.flag_col].to_list(), dtype=np.uint8)
            sdssc2bv = table.attrs.get(self.meta_col, None)

        else:
            raise TypeError(f"Unsupported type for search_table: {type(table)}")

        return ids, flags, sdssc2bv

    def _match_ids(self, ids):
        idx = np.searchsorted(self.ids, ids)
        valid = (ids > 0) & (idx < len(self.ids))

        pos = np.flatnonzero(valid)
        matched = np.zeros_like(valid, dtype=bool)
        matched[pos] = (self.ids[idx[pos]] == ids[pos])

        return idx, valid & matched

    def _set_ids(self, ids, updated):
        idx, valid = self._match_ids(ids)

        for out_i, src_i in zip(np.flatnonzero(valid), idx[valid]):
            updated[out_i] = self.flags[src_i].tolist()   
        return updated     

    def set(self, arr, nrows):
        # dict batch from stream_writer
        if isinstance(arr, dict):
            ids = np.asarray(arr[self.key])
            current = arr[self.flag_col]

            # Convert current values to a mutable Python list-of-lists
            if hasattr(current, "tolist"):
                updated = current.tolist()
            else:
                updated = list(current)

            updated = self._set_ids(ids, updated)

            arr[self.flag_col] = updated
            return arr

        # Arrow Table / RecordBatch-like input
        if hasattr(arr, "column_names") and self.key in arr.column_names:
            ids = np.asarray(arr[self.key])
            current = arr[self.flag_col]
            updated = current.to_pylist() if hasattr(current, "to_pylist") else list(current)

            updated = self._set_ids(ids, updated)

            arr[self.flag_col] = pa.array(updated)
            return arr

        # NumPy structured array / recarray
        if hasattr(arr, "dtype") and arr.dtype.names and self.key in arr.dtype.names:
            ids = np.asarray(arr[self.key])
            current = np.asarray(arr[self.flag_col])
            updated = current.tolist()

            updated = self._set_ids(ids, updated)

            arr[self.flag_col] = np.asarray(updated, dtype=current.dtype)
            return arr

        raise TypeError(f"Unsupported input type for arr: {type(arr)}")
    
    def modify(self, table):
        if self.sdssc2bv is not None:
            meta = dict(table.schema.metadata or {})
            meta[self.meta_col.encode()] = str(self.sdssc2bv).encode()
            table = table.replace_schema_metadata(meta)
        return table

def update_Targeting_flags(run2d, boss_spectro_redux, schema = None, clobber=False, build_only=False,
                           nobackup=False, release='sdsswork', no_remote=False, V_TARG='*'):

    summary_names.set(boss_spectro_redux,run2d)
    spall_file = summary_names.spAllfile_parquet
    Targeting_file = get_Targeting_file(run2d, boss_spectro_redux=boss_spectro_redux)
    if not ptt.exists(Targeting_file):
        clobber = True
    if clobber:
        splog.info('Determining Updated Flags')
        splog.info(f'Reading {spall_file} to get SDSS_IDs')
        spall = pq.read_table(spall_file, columns = ['SDSS_ID'])
        sdssids = spall['SDSS_ID']
        splog.info(f'Number Rows in loaded File: {len(sdssids)}')

        sdssids = pc.drop_null(sdssids)

        mask = pc.and_kleene(
            pc.is_valid(sdssids),
            pc.and_kleene(
                pc.not_equal(sdssids, -999),
                pc.not_equal(sdssids, 0),
            ),
        )
        sdssids = pc.unique(pc.filter(sdssids, mask)).to_numpy()
        splog.info(f'Number Unique SDSS_IDs: {len(sdssids)}')
        search_table = Table([sdssids], names = ['SDSS_ID'])
        search_table, _junk= get_targetflags(search_table, None, db = (not nodb), 
                                             release=release, no_remote=no_remote, V_TARG=V_TARG)
        for col in ['SDSS5_TARGET_FLAGS','SDSSC2BV']:
            search_table[col].fill_value = 0
            
            try:
                data = search_table[col].data.filled()
            except:
                data = search_table[col].data
            search_table[col] = data.astype('uint8')

        # search_table = getAltCatids(search_table, db = (not nodb), 
        #                             release=release, no_remote=no_remote, V_TARG=V_TARG)

        meta = {'Date':time.ctime(),'RUN2D':run2d}
        splog.info(f'Writing {Targeting_file}')

        write_parquet(search_table, Path(Targeting_file), None,
                schema=None, vo=True, metadata=meta,
                column_null=None)

    else:
        if ptt.exists(Targeting_file):
            splog.info('Reading: '+Targeting_file)
            search_table = pd.read_parquet(Targeting_file)

        else:
            raise Exception(Targeting_file+' is missing')


    if build_only:
        exit()
    updater = TargetFlagsUpdater(search_table=search_table)

    epochs={'daily':{},
            #'epoch':dict(epoch = True),
            #'allepoch':dict(custom='allepoch',allsky=True)
            }

    for epoch in epochs.keys():
        summ = Summary_names()
        summ.set(boss_spectro_redux, run2d, **epochs[epoch])

        spAll_schema.yanny_file = summary_names.datamodel
        spAll_schema.datamodel_arrow_schema('EXT1')
        spAll_schema.datamodel_header_metadata('HDR0')

        spLite_schema.yanny_file = summary_names.datamodel
        spLite_schema.datamodel_arrow_schema('EXTLITE')
        spLite_schema.datamodel_header_metadata('HDR0')
        spLite_schema.get_mapping()

        for ff, schema in [(summ.spAllfile_parquet, spAll_schema), (summ.spAlllitefile_parquet, spLite_schema)]:
                              
            if ptt.exists(ff):
                temp_ff = Path(ff)
                temp_ff = temp_ff.parent / f'tmp_{temp_ff.name}'
                temp_ff = str(temp_ff)
                if not nobackup:
                    if ptt.exists(ff+'.bak'):
                        splog.info(f'Removing {ptt.basename(ff)}.bak')
                        remove(ff+'.bak')
                    try:
                        splog.info(f'Saving a backup of {ff} to {ff}.bak')
                        shutil.copy(ff,ff+'.bak')
                    except OSError as exc:
                        raise OSError(f"Failed to save backup to destination {ff}") from exc


                dataset = ds.dataset(ff, format="parquet", partitioning="hive")

                temp_parquet = Path(summary_names.MJD_dir).parent / Path(ff).name
                if temp_parquet.resolve() == Path(ff).resolve():
                    temp_parquet = temp_parquet.with_suffix(".temp.parquet")

                meta = {'SDSSC2BV':updater.sdssc2bv, 'Date':time.ctime(),'RUN2D':run2d}

                stream_writer(dataset, None, temp_parquet, ff,
                            schema, updater.set, hdr=meta,
                            print_func = print)



if __name__ == "__main__":
    update_Targeting_flags(os.getenv('RUN2D'), os.getenv('BOSS_SPECTRO_REDUX'), schema = None, clobber=True,
                           nobackup=True, release='sdsswork', no_remote=False, V_TARG='*', build_only=False)