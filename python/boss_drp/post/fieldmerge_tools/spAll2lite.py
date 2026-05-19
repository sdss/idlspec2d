import pyarrow as pa
import pyarrow.compute as pc
import numpy as np

def list_transform_fallback(arr, func):
    py = arr.to_pylist()
    out = []
    for sub in py:
        if sub is None:
            out.append(None)
        else:
            out.append([func(x) if x not in ("", None) else None for x in sub])
    return pa.array(out)

def list_min(arr):
    py = arr.to_pylist()
    out = []
    for sub in py:
        if not sub:
            out.append(None)
        else:
            vals = [v for v in sub if v is not None]
            out.append(min(vals) if vals else None)
    return pa.array(out, type=pa.int64())

def list_mean(arr):
    py = arr.to_pylist()
    out = []
    for sub in py:
        if not sub:
            out.append(None)
        else:
            vals = [v for v in sub if v is not None]
            out.append(sum(vals)/len(vals) if vals else None)
    return pa.array(out, type=pa.float64())

def spAll_toLite(spAll, spLite_schema):
    # preserve original field objects (they contain metadata)
    original_fields = list(spAll.schema)    # list of pa.Field (metadata included)
    # original_names = [f.name for f in original_fields]
    original_schema_meta = spAll.schema.metadata
    original_arrays = [spAll.column(i) for i in range(spAll.num_columns)]

    new_arrays = []
    new_names = []
    new_fields = []

    # ---------- 1. ID columns ----------
    for col in spLite_schema.id_cols:
        str_col = pc.cast(spAll[col], pa.string())
        lists = pc.split_pattern(str_col, ' ')

        # Split strings and treat empty as null
        try:
            lists = pc.list_transform(lists, lambda x: pc.if_else(x == '', None, x))
        except:
            lists = list_transform_fallback(lists, lambda x: x)

        if col == 'CARTON_TO_TARGET_PK':
            try:
                arr = pc.list_get(lists, 0)
            except:
                arr = pc.list_element(lists,0)
            arr = pc.cast(arr, pa.int64())
        else:
            try:
                int_lists = pc.list_transform(lists, lambda x: pc.cast(x, pa.int64()))
            except:
                int_lists = list_transform_fallback(lists, int)
            try:
                arr = pc.list_min(int_lists)
            except:
                arr = list_min(int_lists)

        # Replace nulls
        mapping_row = spLite_schema.mapping[spLite_schema.mapping['spAllColumn'] == col]
        null = mapping_row['nullval'][0]
        arr = pc.if_else(pc.is_null(arr), pa.scalar(int(null), type=arr.type), arr)

        new_arrays.append(arr)
        new_name = mapping_row['spLiteColumn'][0]
        new_names.append(new_name+'_LITE')

        meta = spLite_schema.column_meta[col]
        new_fields.append(pa.field(new_name+'_LITE', arr.type, metadata=meta))

    # ---------- 2. Numeric columns ----------
    for col in spLite_schema.numeric_cols:
        str_col = pc.cast(spAll[col], pa.string())

        lists = pc.split_pattern(str_col, ' ')
        try:
            lists = pc.list_transform(lists, lambda x: pc.if_else(x == '', None, x))
            float_lists = pc.list_transform(lists, lambda x: pc.cast(x, pa.float64()))
        except:
            lists = list_transform_fallback(lists, lambda x: x)
            float_lists = list_transform_fallback(lists, float)     
        try:
            arr = pc.list_mean(float_lists)
        except:
            arr = list_mean(float_lists)

        # Replace nulls with np.nan
        mapping_row = spLite_schema.mapping[spLite_schema.mapping['spAllColumn'] == col]
        null = mapping_row['nullval'][0]
        null = np.nan if null == 'nan' else float(null)

        arr = pc.if_else(pc.is_null(arr), pa.scalar(null, type=arr.type), arr)

        new_arrays.append(arr)
        new_name = mapping_row['spLiteColumn'][0]
        new_names.append(new_name.replace('_LIST', '')+'_LITE')

        meta = spLite_schema.column_meta[col.replace('_LIST', '')]
        new_fields.append(pa.field(new_name.replace('_LIST', '')+'_LITE', arr.type, metadata=meta))


    # ---------- 3. Build a new table with original columns + new columns ----------
    # Get existing arrays (preserve order)
    existing_arrays = [spAll.column(i) for i in range(spAll.num_columns)]


    # Build final schema from original fields + new_fields (preserves metadata)
    final_fields = original_fields + new_fields
    final_schema = pa.schema(final_fields, metadata=original_schema_meta)

    all_arrays = original_arrays + new_arrays

    if len(all_arrays) != len(final_fields):
        raise ValueError("Number of arrays doesn't match number of final schema fields")


    # Recreate table with explicit schema (keeps field metadata)
    spAll = pa.Table.from_arrays(all_arrays,schema=final_schema)

    return(spAll)