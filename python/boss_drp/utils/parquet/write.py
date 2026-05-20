
from .pandas import table_to_pandas_safe, pandas_dtype_map, pandas_to_votable_dtype
from .utils import add_astropy_string_metadata

import io
import pandas as pd
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq


from typing import Callable, Optional
from astropy.io.votable import writeto
from astropy.io.votable.tree import VOTableFile, Resource, TableElement as voTable, Field, Info

def _as_bytes_map(d):
    out = {}
    for k, v in d.items():
        kb = k if isinstance(k, bytes) else str(k).encode("utf-8")
        vb = v if isinstance(v, bytes) else str(v).encode("utf-8")
        out[kb] = vb
    return out

def write_parquet(table, parquet_path,
                  modifier: Optional[Callable[[pa.Table], pa.Table]],
                  schema=None, vo=True, metadata=None,
                  column_null=None):

    metadata = metadata or {}
    column_null = column_null or {}

    pandas_df = table_to_pandas_safe(table)

    if schema is not None:
        column_null = {
            name: val for name, meta in schema.column_meta.items()
            if (val := meta.get("null")) not in ("", None)
        }
        for col in column_null:
            dtype_str = str(pandas_df[col].dtype)
            if dtype_str in pandas_dtype_map:
                pandas_df[col] = pandas_df[col].astype(pandas_dtype_map[dtype_str])
                pandas_df[col] = pandas_df[col].replace(int(column_null[col]), pd.NA)

        arrow_table = pa.Table.from_pandas(
            pandas_df,
            schema=schema.schema_def,
            preserve_index=False,
        )
    else:
        for col in column_null:
            pandas_df[col] = pandas_df[col].astype("Int64")
            pandas_df[col] = pandas_df[col].replace(column_null[col], pd.NA)
        pl_df = pl.from_pandas(pandas_df)
        arrow_table = pl_df.to_arrow()

    # Non-VO metadata you already support
    meta = {}
    if schema is not None:
        for card in schema.primary_hdr:
            if card in metadata:
                meta[card] = metadata[card]
        if 'fits_hashes' in metadata:
            meta['fits_hashes'] = metadata['fits_hashes']

    if not vo:
        encoded_meta = _as_bytes_map(meta)
        existing_meta = arrow_table.schema.metadata or {}
        combined_meta = {**existing_meta, **encoded_meta}

        if modifier is not None:
            arrow_table = modifier(arrow_table)


        arrow_table = arrow_table.replace_schema_metadata(combined_meta)

        pq.write_table(
            arrow_table,
            parquet_path,
            compression="zstd",
            row_group_size=500_000,
        )
        return

    # Build the VO VOTable header
    votable = VOTableFile(version="1.4")
    resource = Resource()
    votable.resources.append(resource)

    vo_table = voTable(votable)
    vo_table.name = parquet_path.name
    resource.tables.append(vo_table)

    for key, value in meta.items():
        vo_table.infos.append(Info(name=key, value=value))

    for col in pandas_df.columns:
        votype, arraysize = pandas_to_votable_dtype(pandas_df[col])
        field_kwargs = {"name": col, "datatype": votype}
        if arraysize is not None:
            field_kwargs["arraysize"] = arraysize

        field = Field(votable, **field_kwargs)
        if schema is not None:
            col_meta = schema.column_meta.get(col, {}) if isinstance(schema.column_meta, dict) else {}
            if "description" in col_meta:
                field.description = col_meta["description"]
            if "unit" in col_meta:
                field.unit = col_meta["unit"]
            if "ucd" in col_meta:
                field.ucd = col_meta["ucd"]

        vo_table.fields.append(field)

    # Serialize the data-less VOTable into bytes
    buf = io.BytesIO()
    writeto(votable, buf)
    votable_xml = buf.getvalue()

    # Standard VOParquet footer metadata
    vo_footer_meta = {
        b"IVOA.VOTable-Parquet.version": b"1.0",
        b"IVOA.VOTable-Parquet.content": votable_xml,
    }

    if modifier is not None:
        arrow_table = modifier(arrow_table)

    arrow_table = add_astropy_string_metadata(arrow_table)

    existing_meta = arrow_table.schema.metadata or {}
    combined_meta = {**existing_meta, **vo_footer_meta}
    arrow_table = arrow_table.replace_schema_metadata(combined_meta)

    pq.write_table(
        arrow_table,
        parquet_path,
        compression="zstd",
        row_group_size=500_000,
    )

