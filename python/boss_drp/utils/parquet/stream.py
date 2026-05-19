from .utils import add_astropy_string_metadata, fill_null_safe, add_astropy_string_metadata
from .VOparquet import convert_parquet_to_voparquet

from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import json


from typing import Callable, Optional
import shutil

#from tqdm.auto import tqdm
from boss_drp.utils.tqdm import progress as tqdm
def stream_writer(dataset: pa.dataset, columns,
                  temp_parquet: Path,
                  vo_parquet: Path, 
                  schema,   
                  aux: Optional[Callable[[dict],dict]],               
                  batch_size=200_000,
                  hdr= None, print_func = print, freeze_output = None
                  ):

    # ============================================================
    # Write ONE enriched parquet with large row groups
    # ============================================================


    writer = None
    buffer_tables = []
    buffered_rows = 0
    hdr = hdr or {}
    
    scanner = dataset.scanner(columns=columns) 

    pbar = tqdm(dataset, total= dataset.count_rows(), desc="Streaming to Parquet")

    for batch in scanner.to_batches():

        nrows = batch.num_rows
        arr = {}

        # Fill nulls
        if schema is not None:
            for col in batch.schema.names:
                arr[col] = fill_null_safe(batch[col], col, schema.column_meta)
        else:
            for col in batch.schema.names:
                arr[col] = batch[col]

        if aux is not None:
            arr = aux(arr, nrows)

        table = pa.Table.from_pydict(arr)
        buffer_tables.append(table)
        buffered_rows += nrows

        # Flush when reaching desired row group size
        if buffered_rows >= batch_size:
            combined = pa.concat_tables(buffer_tables)
            encoded_meta = {k.encode(): str(v).encode() for k, v in hdr.items()}
            hdr_schema = combined.schema.with_metadata(encoded_meta)

            if writer is None:
                writer = pq.ParquetWriter(
                    temp_parquet,
                    hdr_schema,
                    compression="zstd",
                    use_dictionary=True,
                )
            combined = combined.cast(hdr_schema)
            combined = add_astropy_string_metadata(combined)
            writer.write_table(combined)
            buffer_tables = []
            buffered_rows = 0
        
        pbar.update(batch.num_rows)
    pbar.close()

    # Final flush
    if buffer_tables:
        combined = pa.concat_tables(buffer_tables)
        encoded_meta = {k.encode(): str(v).encode() for k, v in hdr.items()}
        hdr_schema = combined.schema.with_metadata(encoded_meta)

        if writer is None:
            writer = pq.ParquetWriter(
                        temp_parquet,
                        hdr_schema,
                        compression="zstd",
                        use_dictionary=True,
                    )    
        combined = combined.cast(hdr_schema)
        combined = add_astropy_string_metadata(combined)
        writer.write_table(combined)

    if writer:
        writer.close()


    if freeze_output:
        freeze_output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(temp_parquet, freeze_output)
        full_table = pq.read_table(temp_parquet)

        col_mjd = full_table["MJD"]
        col_obs = full_table["OBS"]

        # Normalize OBS (lower + optional trim)
        obs_norm = pc.utf8_lower(pc.utf8_trim_whitespace(col_obs))

        min_val = pc.min(col_mjd).as_py()
        max_val = pc.max(col_mjd).as_py()
        mjds = pc.unique(col_mjd).to_pylist()
        mapping = {freeze_output.name: {'MJD': (min_val, max_val), 'MJDs': mjds}}

        # Get unique OBS values (drop nulls if needed)
        unique_obs = pc.unique(obs_norm)

        for obs_val in unique_obs.to_pylist():
            if obs_val is None:
                continue  # skip nulls

            mask = pc.equal(obs_norm, obs_val)
            obs_mjds = pc.unique(pc.filter(col_mjd, mask)).to_pylist()

            key = f"{obs_val.upper()}_MJDs"
            mapping[freeze_output.name][key] = obs_mjds


        with open(freeze_output.parent / freeze_output.with_suffix('.json'), "w") as f:
            json.dump(mapping, f)
        print_func(f"Frozen Parquet written to {freeze_output}")
        try:
            temp_parquet.unlink()
        except Exception:
            pass
        print_func("Done.")
        return

    print_func(f"Parquet written to {temp_parquet}")

    print_func("Done.")
    # ============================================================
    # Second streamed pass: add VOParquet footer metadata
    # ============================================================
    print_func("Converting to VOParquet...")
    convert_parquet_to_voparquet(
        temp_parquet,
        vo_parquet,
        column_meta=schema.column_meta,
        hdr=hdr,
    )

    try:
        temp_parquet.unlink()
    except Exception:
        pass

    print_func(f"VOParquet written to {vo_parquet}")
    print_func("Done.")


