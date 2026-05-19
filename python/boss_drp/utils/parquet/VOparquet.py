
import io
from pathlib import Path
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from astropy.io.votable import writeto
from vo_parquet.vo_parquet_table import VOParquetTable
from astropy.io.votable.tree import VOTableFile, Resource, TableElement as voTable, Field, Param, Info


def _arrow_field_to_votable(pa_field: pa.Field, col_meta: dict | None = None):
    """
    Infer VOTable datatype + arraysize from a pyarrow Field,
    with optional overrides from column metadata.
    """
    col_meta = col_meta or {}
    t = pa_field.type
    arraysize = None

    # Unwrap list-like columns
    if pa.types.is_fixed_size_list(t):
        arraysize = str(t.list_size)
        t = t.value_type
    elif pa.types.is_list(t) or pa.types.is_large_list(t):
        arraysize = "*"
        t = t.value_type

    # Dictionary-encoded columns
    if pa.types.is_dictionary(t):
        t = t.value_type

    # Explicit overrides win
    if col_meta.get("votable_datatype"):
        datatype = col_meta["votable_datatype"]
    elif pa.types.is_boolean(t):
        datatype = "boolean"
    elif pa.types.is_integer(t):
        datatype = "long"
    elif pa.types.is_floating(t) or pa.types.is_decimal(t):
        datatype = "double"
    else:
        # strings, binary, timestamps, dates, anything else
        datatype = "char"

    if col_meta.get("arraysize") is not None:
        arraysize = col_meta["arraysize"]
    elif datatype == "char" and arraysize is None:
        arraysize = "*"

    return datatype, arraysize



def _build_votable_footer(schema: pa.Schema, column_meta: dict, hdr: dict, table_name: str):
    """
    Build the VOTable XML blob that goes into the VOParquet footer.
    """
    votable = VOTableFile(version="1.4")
    resource = Resource()
    votable.resources.append(resource)

    vo_table = voTable(votable)
    vo_table.name = table_name

    if hdr.get("description"):
        vo_table.description = str(hdr["description"])

    # Store free-form header metadata as PARAMs
    for k, v in hdr.items():
        if k == "description":
            continue
        vo_table.params.append(
            Param(
                votable,
                name=str(k),
                datatype="char",
                arraysize="*",
                value=str(v),
            )
        )

    # Store per-column metadata as FIELDs
    for pa_field in schema:
        meta = column_meta.get(pa_field.name, {}) if isinstance(column_meta, dict) else {}
        datatype, arraysize = _arrow_field_to_votable(pa_field, meta)

        kwargs = {
            "name": pa_field.name,
            "datatype": datatype,
        }
        if arraysize is not None:
            kwargs["arraysize"] = arraysize

        field = Field(votable, **kwargs)

        if meta.get("description"):
            field.description = meta["description"]
        if meta.get("unit"):
            field.unit = meta["unit"]
        if meta.get("ucd"):
            field.ucd = meta["ucd"]

        vo_table.fields.append(field)

    resource.tables.append(vo_table)

    buf = io.BytesIO()
    writeto(votable, buf)
    return buf.getvalue()


def convert_parquet_to_voparquet(src_parquet, dst_parquet, column_meta, hdr=None):
    """
    Rewrite an existing plain Parquet file into a VOParquet-compatible file.
    This is a second streamed pass, not an in-place edit.
    """
    hdr = hdr or {}

    src_parquet = Path(src_parquet)
    dst_parquet = Path(dst_parquet)

    pf = pq.ParquetFile(src_parquet)
    base_schema = pf.schema_arrow

    # Keep any existing footer KV metadata and add the VOParquet blob
    existing_meta = dict(base_schema.metadata or {})
    vo_blob = _build_votable_footer(
        schema=base_schema,
        column_meta=column_meta,
        hdr=hdr,
        table_name=dst_parquet.name,
    )

    existing_meta = dict(base_schema.metadata or {})
    existing_meta[b"IVOA.VOTable-Parquet.version"] = b"1.0"
    existing_meta[b"IVOA.VOTable-Parquet.content"] = vo_blob

    for k, v in hdr.items():
        if k == "description":
            continue
        existing_meta[str(k).encode("utf-8")] = str(v).encode("utf-8")

    out_schema = base_schema.with_metadata(existing_meta)

    # existing_meta[b"IVOA.VOTable-Parquet.content"] = vo_blob
    # for k, v in hdr.items():
    #     if k == "description":
    #         continue
    #     existing_meta[str(k).encode()] = str(v).encode()

    # out_schema = base_schema.with_metadata(existing_meta)

    writer = pq.ParquetWriter(
        dst_parquet,
        out_schema,
        compression="zstd",
        use_dictionary=True,
    )

    try:
        dataset = ds.dataset(src_parquet, format="parquet")
        for batch in dataset.to_batches():
            writer.write_table(pa.Table.from_batches([batch], schema=base_schema))
    finally:
        writer.close()

    return dst_parquet

