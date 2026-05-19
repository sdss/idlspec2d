from boss_drp.post.build_target_summary import build_table_from_fits, compute_hashes
from boss_drp.post.fieldmerge_tools.spAll2lite import spAll_toLite
from functools import partial
from boss_drp.utils.parquet.schema import Schema
from boss_drp.summary import summary_names
from boss_drp.utils.parquet.write import write_parquet

import os
from pathlib import Path
import time
import json



spAll_schema = Schema()
spLite_schema = Schema()

spAll_schema.yanny_file = summary_names.datamodel
spAll_schema.datamodel_arrow_schema('EXT1')
spAll_schema.datamodel_header_metadata('HDR0')
spLite_schema.yanny_file = summary_names.datamodel
spLite_schema.datamodel_arrow_schema('EXTLITE')
spLite_schema.datamodel_header_metadata('HDR0')
spLite_schema.get_mapping()


def master_toparquet():

    run2d = 'master'
    obs = 'manual'
    spAll_fits = Path(os.getenv('BOSS_SPECTRO_REDUX'))/'summary'/'daily','spAll-master.fits.gz'
    spAll_fits_hashes = compute_hashes([spAll_fits])
    table, hdr = build_table_from_fits(spAll_fits)
    meta_spall = {'SDSSC2BV':hdr['SDSSC2BV'], 'Date':time.ctime(),'RUN2D':run2d, 
                            "fits_hashes": json.dumps(spAll_fits_hashes)}
    mjd = max(table['MJD'])
    modifier = partial(spAll_toLite, spLite_schema=spLite_schema)

    _args = dict(run2d=run2d, mjd=f'{mjd}',obs=f'{obs}')
    frozen_partition_dir = Path(summary_names.outdir) / 'MJD' 

    spAll_parquet_path = frozen_partition_dir /f'{obs}'/'spAll'/ summary_names.daily_spAll_parquet.format(**_args)


    write_parquet(table, spAll_parquet_path, modifier,  
                    spAll_schema, metadata=meta_spall, vo=False)
    
