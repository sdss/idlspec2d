#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Sean Morrison (smorris0@illinois.edu)
# @Date: 2026-03-17

from __future__ import annotations


from boss_drp.utils.splog import splog

import fitsio
import pyarrow as pa
import pyarrow.parquet as pq
# from tqdm.auto import tqdm
from boss_drp.utils.tqdm import progress as tqdm
import numpy as np
from typing import Tuple, Dict, Type


def arrow_scalar_dtype(t: pa.DataType) -> Type[np.generic]:
    if pa.types.is_int8(t): return np.int8
    if pa.types.is_int16(t): return np.int16
    if pa.types.is_int32(t): return np.int32
    if pa.types.is_int64(t): return np.int64
    if pa.types.is_uint8(t): return np.uint8
    if pa.types.is_uint16(t): return np.uint16
    if pa.types.is_uint32(t): return np.uint32
    if pa.types.is_uint64(t): return np.uint64
    if pa.types.is_float32(t): return np.float32
    if pa.types.is_float64(t): return np.float64
    if pa.types.is_boolean(t): return np.bool_
    raise TypeError(f"Unsupported scalar type: {t}")

def infer_widths(pf: pq.ParquetFile, batch_size: int = 250_000) -> Tuple[Dict, Dict]:
    splog.info("Calculating String + Array Column Lengths")

    string_cols = [
        f.name for f in pf.schema_arrow
        if pa.types.is_string(f.type) or pa.types.is_large_string(f.type)
    ]

    list_cols = [
        f.name for f in pf.schema_arrow
        if pa.types.is_list(f.type) or pa.types.is_large_list(f.type)
    ]

    all_cols = string_cols + list_cols

    string_widths = {}
    list_widths = {}

    total_rows = pf.metadata.num_rows
    pbar = tqdm(pf, total=total_rows, desc="Inferring column widths")

    for batch in pf.iter_batches(batch_size=batch_size, columns=all_cols):
        schema = batch.schema

        for i, field in enumerate(schema):
            name = field.name
            col = batch.column(i)
            t = field.type

            # ---- STRING ----
            if name in string_cols:
                arr = col.to_pylist()
                mx = max(
                    (0 if v is None else len(str(v).encode("utf-8")))
                    for v in arr
                )
                string_widths[name] = max(string_widths.get(name, 1), mx)

            # ---- LIST ----
            elif name in list_cols:
                values = col.to_pylist()
                for v in values:
                    if v is None:
                        continue
                    list_widths[name] = max(list_widths.get(name, 0), len(v))

        pbar.update(batch.num_rows)

    pbar.close()

    return string_widths, list_widths

def build_dtype(schema: pa.Schema,
                string_widths: dict[str, int],
                list_widths: dict[str, int]) -> np.dtype:
    fields = []

    for field in schema:
        t = field.type

        if pa.types.is_string(t) or pa.types.is_large_string(t):
            w = string_widths.get(field.name, 1)
            fields.append((field.name, f"S{w}"))
            continue

        if pa.types.is_fixed_size_list(t):
            base = arrow_scalar_dtype(t.value_type)
            fields.append((field.name, (base, (t.list_size,))))
            continue

        if pa.types.is_list(t) or pa.types.is_large_list(t):
            base = arrow_scalar_dtype(t.value_type)
            w = list_widths.get(field.name)
            if w is None:
                raise ValueError(f"Could not infer width for list column {field.name}")
            fields.append((field.name, (base, (w,))))
            continue

        fields.append((field.name, arrow_scalar_dtype(t)))

    return np.dtype(fields)

def batch_to_recarray(batch: pa.RecordBatch, dtype: np.dtype) -> np.recarray:
    out = np.empty(batch.num_rows, dtype=dtype).view(np.recarray)

    for i, field in enumerate(batch.schema):
        name = field.name
        col = batch.column(i)
        t = field.type

        if pa.types.is_string(t) or pa.types.is_large_string(t):
            w = dtype.fields[name][0].itemsize
            arr = np.empty(batch.num_rows, dtype=f"S{w}")
            py = col.to_pylist()
            for j, v in enumerate(py):
                arr[j] = b"" if v is None else str(v).encode("utf-8")[:w]
            out[name] = arr
            continue

        if pa.types.is_fixed_size_list(t):
            py = col.to_pylist()
            out[name] = np.asarray(py, dtype=arrow_scalar_dtype(t.value_type))
            continue

        if pa.types.is_fixed_size_list(t) or pa.types.is_list(t) or pa.types.is_large_list(t):
            # Variable-length lists are awkward for FITS binary tables.
            # This will work only if every row has the same length.
            # py = col.to_pylist()
            # out[name] = np.asarray(py, dtype=arrow_scalar_dtype(t.value_type))
            # continue
            base = arrow_scalar_dtype(t.value_type)
            width = dtype.fields[name][0].shape[0]

            py = col.to_pylist()
            out_arr = np.zeros((batch.num_rows, width), dtype=base)

            for j, v in enumerate(py):
                if v is None:
                    continue

                v_arr = np.asarray(v, dtype=base)
                n = len(v_arr)

                if n > width:
                    # truncate (should not happen if inference was correct)
                    v_arr = v_arr[:width]
                    n = width

                if name.upper() == "SDSS5_TARGET_FLAGS":
                    # pad on right (your special case)
                    out_arr[j, :n] = v_arr
                else:
                    # pad on left
                    out_arr[j, width - n:] = v_arr

            out[name] = out_arr
            continue
        arr = col.to_numpy(zero_copy_only=False)
        out[name] = arr

    return out


def parquet_to_fits(
        parquet_path: str,
        fits_path: str,
        batch_size: int = 250_000,
        hdr_cards: list = [],
        string_widths: Dict | None = None,
        list_widths: Dict | None = None,
        column_desc: Dict | None = None,
        column_null: Dict | None = None,
        ) -> Tuple[Dict, Dict]:
    

    pf = pq.ParquetFile(parquet_path)
    if (list_widths is None) or (string_widths is None):
        string_widths, list_widths = infer_widths(pf, batch_size=batch_size)
    dtype = build_dtype(pf.schema_arrow, string_widths, list_widths)

    total_rows = pf.metadata.num_rows
    pbar = tqdm(total=total_rows, desc="Writing FITS")

    primaryhdr = {} #TODO: get from parquet metadata
    for card in hdr_cards:
        try:
            primaryhdr[card] = pf.metadata[card] 
        except:
            primaryhdr[card] = ''

    with fitsio.FITS(fits_path, "rw", clobber=True) as f:
        f.write(None, header=primaryhdr)
        first = True
        for batch in pf.iter_batches(batch_size=batch_size):
            rec = batch_to_recarray(batch, dtype)
            if first:
                f.write(rec)
                first = False
            else:
                f[-1].append(rec)

            pbar.update(batch.num_rows)
        pbar.close()

        if column_desc is not None:
            hdu = f[-1]
            hdr = hdu.read_header()

            

            records = []
            for rec in hdr.records():
                records.append(dict(rec))

                name = rec["name"]
                if name.startswith("TTYPE"):
                    colname = rec["value"]

                    # keep comment on TTYPE
                    desc = column_desc.get(colname)
                    if desc is not None:
                        records[-1]["comment"] = desc
                    else:
                        records[-1]["comment"] = " "
                    # insert matching TNULL immediately after this TTYPE
                    cnull = column_null.get(colname)
                    if cnull is not None:
                        records.append({
                            "name": name.replace("TTYPE", "TNULL"),
                            "value": int(cnull),
                            "comment": f'Null value for {colname}'
                        })
                if name.startswith("TFORM"):
                    records[-1]["comment"] = ""

            hdu.write_keys(fitsio.FITSHDR(records), clean=False)

    return dict(string_widths=string_widths, list_widths=list_widths)

