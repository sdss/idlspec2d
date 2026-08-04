import re

from boss_drp.summary import Summary_names, summary_names, fieldlist_name
from boss_drp.post.fieldmerge import build_custom_fieldlist
from boss_drp.field import fieldgroup
from boss_drp.utils.splog import splog
from boss_drp.utils import jdate
from boss_drp.utils.parquet.to_fits import parquet_to_fits
from boss_drp.utils.parquet.write import write_parquet
from boss_drp.utils.parquet.table import build_table_from_fits
# from boss_drp.utils.parquet.VOparquet import convert_parquet_to_voparquet
from boss_drp.utils.parquet.utils import fill_null_safe
from boss_drp.utils.parquet.schema import Schema
from boss_drp.utils.parquet.stream import stream_writer
from boss_drp.post.fieldmerge_tools.spAll2lite import spAll_toLite
from boss_drp.post import get_Targeting_file, TargetFlagsUpdater
from boss_drp.post import plot_sky_targets, plot_sky_locations

import json
import hashlib
import time
from datetime import date, timedelta
from pathlib import Path
import numpy as np
import os
from functools import partial
import textwrap
import shutil
import inspect
from parse import parse
import json
import re

import pyarrow.parquet as pq
import pyarrow.dataset as ds
import pyarrow.compute as pc


# import pandas as pd
from astropy.io import fits
from astropy.time import Time



def clean_wrap(text, prefix='', pad=0, width = 100):
    total_prefix_len = len(prefix) #+ pad

    return prefix + textwrap.fill(
        text,
        width=width - pad # - total_prefix_len,
        #subsequent_indent=' ' * total_prefix_len
    )

# ----------------------------
# Hashing utility
# ----------------------------
def file_md5(path, chunk_size=8_388_608):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()

def compute_hashes(fits_files):
    return {str(f): file_md5(f) for f in fits_files}

spAll_schema = Schema()
spLine_schema = Schema()
spLite_schema = Schema()
spAll_freeze_schema = Schema()


class SpecPrimary:
    def __init__(self):
        self.counts = {}
        self.best = {}
        self.row_counter = 0

    def compute(self, dataset, schema):
        # ============================================================
        # compute NSPECOBS + best row index
        # ============================================================
        self.counts = {}
        self.best = {}

        splog.info("Computing NSPECOBS and primary rows...")
        row_counter = 0

        for batch in dataset.to_batches():

            sdssid = fill_null_safe(batch["SDSS_ID"], "SDSS_ID", schema.column_meta)
            mjd = fill_null_safe(batch["MJD"], "MJD", schema.column_meta)

            sn = batch["SN_MEDIAN"]
            sn_arr = np.stack(sn.to_numpy(zero_copy_only=False), axis=0)
            lengths = pc.list_value_length(sn).to_numpy()
            if not np.all(lengths == lengths[0]):
                raise ValueError("Inconsistent list lengths in SN_MEDIAN")

            jfilt = 2 if lengths[0] != 1 else 0
            sn_r = sn_arr[:, jfilt]

            score = sn_r  

            for i in range(len(sdssid)):                
                sid = sdssid[i]
                try:
                    sid = sid.as_py()
                except:
                    pass

                if sid <= 0:
                    row_counter += 1
                    continue

                self.counts[sid] = self.counts.get(sid, 0) + 1
                s = score[i]
                m = mjd[i]

                current = self.best.get(sid)
                if current is None or (s > current[0]) or (s == current[0] and m > current[1]):
                    self.best[sid] = (s, m, row_counter)

                row_counter += 1

        splog.info(f"Tracked {len(self.counts):,} unique SDSS_ID")

    def set(self, arr, nrows):
        sdssid = arr["SDSS_ID"]
        try: sdssid = sdssid.to_numpy(zero_copy_only=False)
        except: pass
        batch_rows = np.arange(self.row_counter, self.row_counter + nrows)

        specprimary = np.zeros(nrows, dtype=np.int16)
        nspecobs = np.zeros(nrows, dtype=np.int16)

        valid_mask = sdssid > 0
        if np.any(valid_mask):
            valid_ids = sdssid[valid_mask]
            nspecobs[valid_mask] = [self.counts[s] for s in valid_ids]
            best_rows = np.array([self.best[s][2] for s in valid_ids])
            specprimary[valid_mask] = (
                batch_rows[valid_mask] == best_rows
            ).astype(np.int16)

        specprimary[~valid_mask] = -999
        nspecobs[~valid_mask] = -999

        arr["SPECPRIMARY"] = specprimary
        arr["NSPECOBS"] = nspecobs

        return arr

specprimary = SpecPrimary()

# ---------------------------
# Main function
# ---------------------------
def stack_all_parquet(
        partition_dir, partition_file_name_model, _args,
        output_parquet, schema, hdr, 
        batch_size=200_000,
        compute_specprimary=True,
        lite = False, freeze_output = False,
        frozen_partition_dir = None, freeze_args={},
        freeze_schema = None,
        bkup = False
):

    temp_parquet = Path(summary_names.MJD_dir).parent / Path(output_parquet).name
    if temp_parquet.resolve() == Path(output_parquet).resolve():
        temp_parquet = temp_parquet.with_suffix(".temp.parquet")
    parquet_files = list(partition_dir.rglob(str(str(partition_file_name_model).format(**_args))))

    if frozen_partition_dir is not None:
        if len(freeze_args) == 0:
            freeze_args = _args  
        frozen_files = list(frozen_partition_dir.rglob(str(str(partition_file_name_model.name).format(**freeze_args))))
        parquet_files.extend(frozen_files)


    splog.info(f"Parquet files found for stacking: {len(parquet_files)}")
    dataset = ds.dataset(parquet_files, format="parquet", partitioning="hive")

    total_rows = dataset.count_rows()
    for key in hdr.keys():
        if key == 'SDSSC2BV':
            try:
                hdr[key] = dataset.meta['SDSSC2BV']
            except:
                continue
        if key == 'Nobj':
            hdr[key] = total_rows
        if key == 'Nlines':
            hdr[key] = len(dataset.names)


    splog.info(f"Total rows: {total_rows:,}")

    # ============================================================
    # compute NSPECOBS + best row index
    # ============================================================
    counts = {}
    best = {}

    if compute_specprimary:
        specprimary.compute(dataset, schema)

    # ============================================================
    # Write ONE enriched parquet with large row groups
    # ============================================================
    



    if lite:
        columns = {name: ds.field(name) for name in dataset.schema.names}
        for col in spLite_schema.mapping['spAllColumn'].tolist():
            columns.pop(col)
        for col in spLite_schema.mapping['spLiteColumn'].tolist():
            columns[f'{col}'] = columns.pop(f'{col}_LITE')
        for col in list(columns.keys()):
            if col not in spLite_schema.column_meta.keys():
                columns.pop(col)
    else:
        columns = {name: ds.field(name) for name in dataset.schema.names}

        for col in list(columns.keys()):
            if col not in schema.column_meta.keys():
                columns.pop(col)
    

    splog.info(f"Writing parquet to ... {output_parquet}")

    

                      

    # if freeze_output:
    #     freeze_output = frozen_partition_dir / str(partition_file_name_model.name).format(**freeze_args)

    if bkup:
        if Path(output_parquet).exists():
            summary_names.bk.mkdir()
            # if not Path(summary_names.bk.spAllfile_parquet).exists():
            shutil.copy2(summary_names.spAllfile_parquet, summary_names.bk.spAllfile_parquet)
            shutil.copy2(summary_names.spAlllitefile_parquet, summary_names.bk.spAlllitefile_parquet)
            shutil.copy2(summary_names.splinefile_parquet, summary_names.bk.splinefile_parquet)

    stream_writer(dataset, columns, temp_parquet, output_parquet,
                  schema, specprimary.set, batch_size=batch_size, hdr=hdr,
                  print_func = splog.info)#, freeze_output = freeze_output)

    if freeze_output:
        if freeze_schema is None:
            freeze_schema = schema
        #IF not spLine
        #TODO this is not working properly because the frozen file is missing the _LITE columns...
        freeze_output = frozen_partition_dir / str(partition_file_name_model.name).format(**freeze_args)
        columns_freeze = {col: ds.field(col) for col in dataset.schema.names}
        stream_writer(dataset, columns_freeze, temp_parquet, freeze_output,
                      freeze_schema, None, batch_size=batch_size, hdr=hdr, 
                      print_func = splog.info, freeze_output=freeze_output)

# ----------------------------
# Rebuild check
# ----------------------------
def parquet_needs_rebuild(parquet_path, current_hashes, current_sdssc2bv, prefix = ''):
    if not parquet_path.exists():
        p = parquet_path

        try:
            rel = p.relative_to(os.getenv('BOSS_SPECTRO_SCRATCH'))
            display = Path("$BOSS_SPECTRO_SCRATCH") / rel
        except ValueError:
            try:
                rel = p.relative_to(os.getenv('BOSS_SPECTRO_REDUX'))
                display = Path("$BOSS_SPECTRO_REDUX") / rel
            except ValueError:
                display = p  # fallback
        splog.error(f"{prefix}Missing {display}")
        return True

    pf = pq.ParquetFile(parquet_path)
    metadata = pf.schema_arrow.metadata

    if metadata is None or b"fits_hashes" not in metadata or b"SDSSC2BV" not in metadata:
        splog.error(f'{prefix}invalid meta for {parquet_path}')
        return True

    stored_hashes = json.loads(
        metadata[b"fits_hashes"].decode()
    )

    stored_sdssc2bv = int(metadata[b"SDSSC2BV"].decode())
    if current_sdssc2bv is None:
        current_sdssc2bv = stored_sdssc2bv

    parquet_path.touch()
    return (stored_hashes != current_hashes) or (stored_sdssc2bv != current_sdssc2bv)

def get_frozen_mjd(path, parquet_name=None):
    if not path.exists():
        return []
    with open(path, 'r') as f:
        mapping = json.load(f)   
        return mapping[parquet_name]['MJDs']
    return[]

def extract_mjd(p, run2d):
    pattern = re.compile(rf'(?<={re.escape(run2d)}-)(\d+)(?=_)')
    m = pattern.search(p.name)
    return int(m.group(1)) if m else None

def get_latest_frozen_mjd(partition_file_name_model, _args, frozen_partition_dir, run2d):
    frozen_files = list(frozen_partition_dir.rglob(str(partition_file_name_model).format(**_args)))
    frozen_mjds = [extract_mjd(f, run2d) for f in frozen_files]
    frozen_mjds = [mjd for mjd in frozen_mjds if mjd is not None]
    if not frozen_mjds:
        return frozen_partition_dir / "null.parquet", None
    _args['mjd'] = max(frozen_mjds)
    return frozen_partition_dir / str(partition_file_name_model).format(**_args), max(frozen_mjds)

def compose(*funcs):
    def inner(x):
        for f in funcs:
            x = f(x)
        return x
    return inner

def build_target_summary(indir,run2d,epoch=False, allsky=False, custom=None, force_rebuild=False, keep_active=False,
                        datamodel=None, line_datamodel=None, logfile = None, outroot = None, dev=False,
                        to_fits=False, mjdstart = None, mjdend = None, clobber=False, MJD_dir=None,
                        update_target_flags = False,  freeze_output= False, bkup = False,
                        *args, **kwrds):
    start = time.time()


    if logfile is None:
        if outroot is not None:
            logfile = Path(outroot).with_suffix('.log')
        else:
            if epoch is True:
                logfile = Path(indir) / run2d / 'summary' / 'epoch' / f'spAll_epoch-{run2d}.log'
            elif custom is None:
                logfile = Path(indir) / run2d / 'summary' / 'daily' / f'spAll-{run2d}.log'
            else:
                logfile = (
                    Path(indir)
                    / run2d
                    / 'summary'
                    / fieldgroup(custom, custom=True)
                    / f'spAll_{custom}-{run2d}.log'
                )

        if dev:
            logfile = Path(str(logfile).replace('spAll', 'spAll_dev'))

    logfile = Path(logfile)
    if logfile.parent != Path():
        logfile.parent.mkdir(parents=True, exist_ok=True)

    splog.open(logfile=logfile, backup=False)
    splog.log(f'Log file {logfile} opened '+ time.ctime())


    if datamodel is not None:
        summary_names.datamodel = datamodel
    if line_datamodel is not None:
        summary_names.line_datamodel = line_datamodel

    spAll_schema.yanny_file = summary_names.datamodel
    spAll_schema.datamodel_arrow_schema('EXT1')
    spAll_schema.datamodel_header_metadata('HDR0')

    spLite_schema.yanny_file = summary_names.datamodel
    spLite_schema.datamodel_arrow_schema('EXTLITE')
    spLite_schema.datamodel_header_metadata('HDR0')
    spLite_schema.get_mapping()

    spLine_schema.yanny_file = summary_names.line_datamodel
    spLine_schema.datamodel_arrow_schema('EXT1')
    spLine_schema.datamodel_header_metadata('HDR0')

    spAll_freeze_schema.yanny_file = summary_names.datamodel
    spAll_freeze_schema.datamodel_arrow_schema('EXT1')
    spAll_freeze_schema.datamodel_header_metadata('HDR0')
    for row in spLite_schema.mapping:
        spAll_freeze_schema.column_meta[row['spLiteColumn']+'_LITE'] = spLite_schema.column_meta[row['spLiteColumn']]

    summary_names.set(indir, run2d, tmpext='.tmp', epoch=epoch, allsky=allsky, 
                      custom=custom, outroot=outroot, dev=dev, MJD_dir=MJD_dir)


    partition_dir = Path(summary_names.MJD_dir)
    partition_dir.mkdir(parents=True, exist_ok=True)

    frozen_partition_dir = Path(summary_names.outdir) / 'mjd' 

    ## Load Fieldlist
    fieldlist_name.build(indir, run2d, epoch=epoch, custom_name=custom)
    fieldlist_file = fieldlist_name.name
    fieldlist_parquet = fieldlist_name.parquet
    if allsky is False:
        fmjd = None
        if Path(fieldlist_parquet).exists():
            splog.info(f'Reading {fieldlist_parquet}')
            try:
                fmjd = pq.read_table(fieldlist_parquet).to_pandas().to_records(index=False)
            except:
                fmjd = None
                time.sleep(90)
                try:
                    fmjd = pq.read_table(fieldlist_parquet).to_pandas().to_records(index=False)
                except:
                    fmjd = None
                    pass
            
        if fmjd is None:
            splog.info(f'Failure loading {fieldlist_parquet}, trying fits fallback')
            if Path(fieldlist_file).exists():
                splog.info(f'Reading {fieldlist_file}')
                try:
                    fmjd = fits.getdata(fieldlist_file)
                except:
                    time.sleep(90)
                    try:
                        fmjd = fits.getdata(fieldlist_file)
                    except:
                        pass
        if fmjd is None:
            splog.error(f'No Valid FieldList file exists at {fieldlist_file} or {fieldlist_parquet}')
            exit()
    else:
        fmjd = build_custom_fieldlist(indir, custom, run2d, run2d)
    
    if mjdstart:
        fmjd = fmjd[fmjd['MJD'] >= mjdstart]
    if mjdend:
        fmjd = fmjd[fmjd['MJD'] <= mjdend]

    if not mjdstart:
        patterns = [summary_names.daily_spAll_parquet.format(run2d=run2d, mjd='*',obs='*'),
                    summary_names.daily_spline_parquet.format(run2d=run2d, mjd='*',obs='*')]
            
    else:
        patterns = []
        for mjd in range(mjdstart, int(float(Time( str(date.today())).jd)-2400000.5)+1):
            patterns.extend([summary_names.daily_spAll_parquet.format(run2d=run2d, mjd=f'{mjd}',obs='*'),
                                summary_names.daily_spline_parquet.format(run2d=run2d, mjd=f'{mjd}',obs='*')])

    if clobber:
        splog.info('Removing old daily spAll parquet files')
        for p in patterns:
            for f in partition_dir.rglob(p):
                if f.is_file():
                    f.unlink()

    f_mjds = {}
    #Find frozen MJDs from frozen_partition_dir and remove any non-frozen parquet files for those MJDs in partition_dir to avoid confusion
    for obs_f in ['APO', 'LCO']:
        f_mjds[obs_f] = []

        f_args = dict(run2d=run2d, mjd='*',obs='manual')
        for sdir, t in [('spAll',summary_names.daily_spAll_parquet), ('spLine',summary_names.daily_spline_parquet)]:
            f_parc, mjd_f = get_latest_frozen_mjd(t, f_args,
                                            frozen_partition_dir/ '{obs}'.format(**f_args) / sdir, run2d)
            f_mjds[obs_f].extend(get_frozen_mjd(f_parc.parent / f_parc.with_suffix('.json'), f_parc.name))
            
        if frozen_partition_dir.exists():
            if not partition_dir.samefile(frozen_partition_dir):
                for sdir, t in [('spAll',summary_names.daily_spAll_parquet), ('spLine',summary_names.daily_spline_parquet)]:
                    f_parc_alt, mjd_f_alt = get_latest_frozen_mjd(t, f_args, 
                                                        partition_dir / '{obs}'.format(**f_args) / sdir, run2d)
                    f_mjds[obs_f].extend(get_frozen_mjd(f_parc_alt.parent / f_parc_alt.with_suffix('.json'), f_parc_alt.name))

        if len(f_mjds[obs_f]) > 0:
            f_mjds_arr = np.asarray(f_mjds[obs_f]).astype(int)
            f_mjds[obs_f] = set(f_mjds_arr)
            splog.info(f'Found {len(f_mjds_arr)} frozen MJDs for {obs_f}: {sorted(f_mjds[obs_f])}')
            for p in patterns:
                for f in partition_dir.rglob(p):
                    mjd = None
                    if f.is_file():
                        result = parse(summary_names.daily_spAll_parquet, f.name)
                        if result is None:
                            result = parse(summary_names.daily_spline_parquet, f.name)
                            if result is None:
                                continue
                        if result['obs'] == 'manual':
                            continue
                        
                        if int(result['mjd']) in f_mjds_arr:
                            f.unlink()
    

    i = 0
    meta_spall = None
    meta_spLine = None
    rebuilt = False

    nmissing = 0
    nfound = 0
    nfrozen = 0
    current_sdssc2bv = None

    if update_target_flags:
        # Update Targeting flags if requested, using the spTargeting file for the run2d which should have the most up-to-date flags.
        #  This will be merged in as part of the spAll_toLite conversion
        sptarget= get_Targeting_file(run2d, boss_spectro_redux=indir)
        splog.info(f'Loading spTargeting file ({sptarget}) for updated flags')
        updater = TargetFlagsUpdater(search_parquet = sptarget)
        current_sdssc2bv = updater.sdssc2bv


    for mjd, obs in sorted(set(zip(fmjd['MJD'], fmjd['OBSERVATORY']))):
        # 
        fmjd_summ = Summary_names()
        spAll_fits_files = []
        spline_fits_files = []
        fields = []
        missing = []
        if mjd in f_mjds[obs]:
            splog.info(f"[{mjd}:{obs}] Found frozen MJD, skipping to next.")
            nfrozen += len(fmjd[(fmjd['MJD'] == mjd) & (fmjd['OBSERVATORY'] == obs)]['FIELD'])
            continue
        for field in fmjd[(fmjd['MJD'] == mjd) & (fmjd['OBSERVATORY'] == obs)]['FIELD']:
            if (custom is not None) and (obs.lower() in ['apo', 'lco']):
                field = f"{custom}_{obs.lower()}"
            elif custom is not None:
                field = f"{custom}"
            else:
                field = int(field)
            fmjd_summ.set(indir, run2d, field=field,mjd=str(mjd), epoch=epoch, allsky=allsky, 
                          custom=custom, outroot=outroot, dev=dev, obs=obs)
    
            found = False
            if os.path.exists(fmjd_summ.spAllfile):
                spAll_fits_files.append(fmjd_summ.spAllfile)
                found = True
            if os.path.exists(fmjd_summ.splinefile):
                spline_fits_files.append(fmjd_summ.splinefile)
                found = True

            if not found:
                missing.append(field)
                nmissing += 1
            else:
                fields.append(field)
                nfound += 1

            i+=1

        padLength = len(inspect.currentframe().f_code.co_name) + 2 
        if (len(spAll_fits_files) + len(spline_fits_files)) == 0:
            splog.info(f"[{mjd}:{obs}] No files found, skipping.")
            splog.info(clean_wrap(f"{', '.join(map(str, missing))}", pad = padLength, 
                                  prefix = f"[{mjd}:{obs}] Missing Fields: "))
                                  
            continue
        _args = dict(run2d=run2d, mjd=f'{mjd}',obs=f'{obs}')
        spAll_parquet_path = partition_dir /f'{obs}'/'spAll'/ summary_names.daily_spAll_parquet.format(**_args)
        spAll_parquet_path.parent.mkdir(parents=True, exist_ok=True)
        spline_parquet_path = partition_dir /f'{obs}'/'spLine'/ summary_names.daily_spline_parquet.format(**_args)
        spline_parquet_path.parent.mkdir(parents=True, exist_ok=True)

        spAll_fits_hashes = compute_hashes(spAll_fits_files)
        spline_fits_hashes = compute_hashes(spline_fits_files)
        if (parquet_needs_rebuild(spAll_parquet_path, spAll_fits_hashes, current_sdssc2bv, prefix=f'[{mjd}:{obs}] ') or
            parquet_needs_rebuild(spline_parquet_path, spline_fits_hashes, current_sdssc2bv, prefix=f'[{mjd}:{obs}] ')):
            if spAll_parquet_path.is_file():
                splog.info(f"[{mjd}:{obs}] Rebuilding parquet...")
            else:
                splog.info(f"[{mjd}:{obs}] Building parquet...")
            fields = map(str, fields)
            splog.info(clean_wrap(f"{', '.join(fields)}", pad = padLength, #= len("build_target_summary: "),
                                  prefix = f"[{mjd}:{obs}] Found Fields: "))
            if len(missing) > 0:
                missing = map(str, missing)
                splog.info(clean_wrap(f"{', '.join(missing)}", pad = padLength, #= len("build_target_summary: "),
                                    prefix = f"[{mjd}:{obs}] Missing Fields: "))

            table, hdr = build_table_from_fits(spAll_fits_files)
            meta_spall = {'SDSSC2BV':hdr['SDSSC2BV'], 'Date':time.ctime(),'RUN2D':run2d, 
                          "fits_hashes": json.dumps(spAll_fits_hashes)}
            

            modifier = partial(spAll_toLite, spLite_schema=spLite_schema)

            if update_target_flags:
                splog.info(f"[{mjd}:{obs}] Updating targeting flags from {sptarget}")
                modifier = compose(
                    partial(spAll_toLite, spLite_schema=spLite_schema),
                    updater.modify, updater.set
                )

            write_parquet(table, spAll_parquet_path, modifier,  
                          spAll_schema, metadata=meta_spall, vo=False)

            table, hdr = build_table_from_fits(spline_fits_files)
            meta_spLine = {'Date':time.ctime(),'RUN2D':run2d, 
                           'NLines':hdr['DIMS0'], 'Nobj':hdr['DIMS1'],
                           "fits_hashes": json.dumps(spline_fits_hashes)}

            write_parquet(table, spline_parquet_path, None,
                          spLine_schema, metadata=meta_spLine, vo=False)

            splog.info(f"[{mjd}:{obs}] Done.")
            rebuilt = True
        else:
            splog.info(f"[{mjd}:{obs}] Up to date. No rebuild needed.")


        # This is just left in here incase a "frozen" parquet is put in scratch, but the default location should be
        # to include these in the BOSS_SPECTRO_REDUX equivalent 
        f_args = dict(run2d=run2d, mjd=f'{mjd}',obs='manual')
        try:

            if not partition_dir.samefile(frozen_partition_dir):
                frozen_spAll = partition_dir/'{obs}'.format(**f_args)/'spAll'/summary_names.daily_spAll_parquet.format(**f_args)
                if frozen_spAll.exists():
                    frozen_spAll.touch()
                frozen_spLine = partition_dir/'{obs}'.format(**f_args)/'spLine'/summary_names.daily_spline_parquet.format(**f_args)
                if frozen_spLine.exists():
                    frozen_spLine.touch()
        except FileNotFoundError:
            pass
    splog.info(f'Found {nfound} Fields ({nmissing} missing and {nfrozen} frozen) and added to intermediate daily parquet files')

    if keep_active:
        if (i == 0) or (mjdstart is not None):
            splog.debug('Updating last read time of all intermediate parquet files (to keep active in scratch)')
            for p in partition_dir.rglob("*"):
                if p.is_file():
                    p.touch()

    if meta_spall is None:
        meta_spall = {'Date':time.ctime(), 'RUN2D':run2d, 'SDSSC2BV':''}
    elif 'fits_hashes' in meta_spall:
        _ = meta_spall.pop('fits_hashes',None)

    if meta_spLine is None:
        meta_spLine = {'Date':time.ctime(),'RUN2D':run2d, 
                       'NLines':'', 'Nobj':''}
    elif 'fits_hashes' in meta_spall:
        _ = meta_spall.pop('fits_hashes',None)

    # Only rebuild global Parquet if something changed
    if ((not Path(summary_names.spAllfile_parquet).exists()) or 
        (not Path(summary_names.spAlllitefile_parquet).exists()) or 
        (not Path(summary_names.splinefile_parquet).exists())):
        splog.info('Global Parquet files do not exist, rebuilding...')
        rebuilt = True   
    if rebuilt or force_rebuild:
        _args = dict(run2d='*', mjd='*',obs='*')
        t_frozen_partition_dir = frozen_partition_dir / '{obs}'.format(**f_args) if frozen_partition_dir is not None else None
        stack_all_parquet(
            partition_dir,  Path("*") / "spAll" / summary_names.daily_spAll_parquet,
            _args, summary_names.spAllfile_parquet, spAll_schema,
            meta_spall, frozen_partition_dir=t_frozen_partition_dir / 'spAll' if frozen_partition_dir is not None else None,
            freeze_output=freeze_output, freeze_args = f_args, bkup=bkup, freeze_schema = spAll_freeze_schema,
        )
        
        stack_all_parquet(
            partition_dir,  Path("*") / "spAll" / summary_names.daily_spAll_parquet,
            _args, summary_names.spAlllitefile_parquet, spLite_schema,
            meta_spall, frozen_partition_dir=t_frozen_partition_dir / 'spAll' if frozen_partition_dir is not None else None,
            compute_specprimary=True, lite=True, freeze_args = f_args,
        )

        stack_all_parquet(
            partition_dir,  Path("*") / "spLine" / summary_names.daily_spline_parquet,
            _args, summary_names.splinefile_parquet, spLine_schema,
            meta_spLine, frozen_partition_dir=t_frozen_partition_dir / 'spLine' if frozen_partition_dir is not None else None,
            compute_specprimary=False, freeze_args = f_args,freeze_output = freeze_output
        )

    if to_fits:
        splog.info('Converting Parquet Summary Files to Fits')
        splog.info('Converting spAll')
        col_widths = parquet_to_fits(summary_names.spAllfile_parquet, summary_names.spAllfile,
                        hdr_cards=list(spAll_schema.primary_hdr.keys()),
                        column_desc={ name: meta["description"] for name, meta in spAll_schema.column_meta.items()},
                        column_null={ name: val for name, meta, in spAll_schema.column_meta.items()
                                     if (val := meta.get("null")) not in ("", None)})


        splog.info('Converting spAll-lite')
        _ = parquet_to_fits(summary_names.spAlllitefile_parquet, summary_names.spAlllitefile,
                        hdr_cards=list(spLite_schema.primary_hdr.keys()), **col_widths,
                        column_desc={ name: meta["description"] for name, meta in spLite_schema.column_meta.items()},
                        column_null={ name: val for name, meta, in spLite_schema.column_meta.items() 
                                     if (val := meta.get("null")) not in ("", None)})

        
        splog.info('Converting spAllLine')
        _ = parquet_to_fits(summary_names.splinefile_parquet, summary_names.splinefile,
                            hdr_cards=list(spLine_schema.primary_hdr.keys()),
                            column_desc={ name: meta["description"] for name, meta in spLine_schema.column_meta.items()},
                            column_null={ name: val for name, meta, in spLine_schema.column_meta.items() 
                                         if (val := meta.get("null")) not in ("", None)})

    if custom is None:
        plot_sky_locations()
        plot_sky_targets(nobs=True)

    splog.info(f'Elapsed Time: {str(timedelta(seconds=time.time()-start))}')
    splog.info('Successful completion of build_spall at '+ time.ctime())

if __name__ == "__main__":
    todaymjd = int(float(Time( str(date.today())).jd)-2400000.5)
    ndays =  None #10
    if ndays is not None:
        mjdstart = todaymjd - ndays
    else:
        mjdstart = None
    
    mjdend = None

    build_target_summary(os.getenv('BOSS_SPECTRO_REDUX'),os.getenv('RUN2D'),
                         epoch=False, allsky=False, custom=None, to_fits=False, 
                         force_rebuild=True,datamodel=None, line_datamodel=None, 
                         mjdstart = mjdstart, mjdend = mjdend, clobber= True, 
                         update_target_flags=True, freeze_output = False, bkup = True) 


