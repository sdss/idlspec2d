import pyarrow as pa
import numpy as np

def add_astropy_string_metadata(table: pa.Table) -> pa.Table:
    md = dict(table.schema.metadata or {})

    for name in table.column_names:
        col = table[name].combine_chunks()
        t = col.type

        if pa.types.is_string(t) or pa.types.is_large_string(t):
            vals = [x.as_py() for x in col if x.as_py() is not None]
            max_len = max((len(v) for v in vals), default=0)
            md[f"table::len::{name}".encode("utf-8")] = str(max_len).encode("utf-8")

        elif pa.types.is_binary(t) or pa.types.is_large_binary(t):
            vals = [x.as_py() for x in col if x.as_py() is not None]
            max_len = max((len(v) for v in vals), default=0)
            md[f"table::len::{name}".encode("utf-8")] = str(max_len).encode("utf-8")

    return table.replace_schema_metadata(md)

def fill_null_safe(column: pa.Array, column_name: str = None, column_meta: dict = None):
    """
    Fill nulls in Arrow column using the sentinel from column_meta if available,
    otherwise fall back to reasonable defaults.
    Returns NumPy array suitable for FITS.
    """
    nrows = len(column)
    typ = column.type

    # Determine sentinel
    if column_meta and column_name and "fill_value" in column_meta[column_name]:
        sentinel = column_meta[column_name]["fill_value"]
    else:
        # Fallback
        sentinel = None

    # --------------------------
    # Fixed-size list
    # --------------------------
    if pa.types.is_fixed_size_list(typ):
        size = typ.list_size
        base_type = typ.value_type

        if sentinel is None:
            if pa.types.is_integer(base_type):
                sentinel_val = 0
            elif pa.types.is_floating(base_type):
                sentinel_val = np.nan
            elif pa.types.is_boolean(base_type):
                sentinel_val = False
            else:
                sentinel_val = ""

        else:
            sentinel_val = sentinel
        fill_row = np.full(size, sentinel_val)
        fill_array = np.tile(fill_row, (nrows, 1))
        arrow_sentinel = pa.array(fill_array.tolist(), type=typ)

        filled_outer = column.fill_null(arrow_sentinel)
        vals = filled_outer.values
        
        if vals.null_count > 0:
            vals = vals.fill_null(sentinel_val)
            
        arr = vals.to_numpy().reshape(-1, size)
        if sentinel_val == 0:
            flat = pa.array(arr.reshape(-1), type=pa.int64())
        else:
            flat = pa.array(arr.reshape(-1))
        arrow_arr = pa.FixedSizeListArray.from_arrays(flat, list_size=size)
        return arrow_arr.to_numpy(zero_copy_only=False)


    # --------------------------
    # Variable-length list
    # --------------------------
    elif pa.types.is_list(typ):
        if sentinel is None:
            sentinel_val = []
        else:
            sentinel_val = sentinel

        arrow_sentinel = pa.array([sentinel_val] * nrows, type=typ)
        return column.fill_null(arrow_sentinel).to_numpy(zero_copy_only=False)

    # --------------------------
    # Scalars (int, float, bool, string)
    # --------------------------
    else:
        if sentinel is None:
            if pa.types.is_integer(typ):
                sentinel_val = 0
            elif pa.types.is_floating(typ):
                sentinel_val = np.nan
            elif pa.types.is_boolean(typ):
                sentinel_val = False
            else:
                sentinel_val = ""
        else:
            sentinel_val = sentinel

        arrow_sentinel = pa.scalar(sentinel_val, type=typ)
        return column.fill_null(arrow_sentinel).to_numpy(zero_copy_only=False)
