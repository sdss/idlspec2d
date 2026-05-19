import numpy as np
import pandas as pd 

pandas_dtype_map = {
    'int8': 'Int8',
    'int16': 'Int16',
    'int32': 'Int32',
    'int64': 'Int64',
    'uint8': 'UInt8',
    'uint16': 'UInt16',
    'uint32': 'UInt32',
    'uint64': 'UInt64',
}


def pandas_to_votable_dtype(series):
    dtype = str(series.dtype)


    def is_array_column(series):
        return series.dtype == "object" and any(
            isinstance(x, (list, np.ndarray)) for x in series.dropna()
        )

    # --- array columns ---
    if is_array_column(series):
        # try to infer element type from first non-null
        first = next((x for x in series if x is not None), None)

        if first is not None:
            arr = np.asarray(first)

            if np.issubdtype(arr.dtype, np.integer):
                return "long", "*"
            elif np.issubdtype(arr.dtype, np.floating):
                return "double", "*"

        return "double", "*"  # fallback

    # --- scalar columns ---
    if dtype.startswith("int") or dtype.startswith("Int"):
        return "long", None
    elif dtype.startswith("float"):
        return "double", None
    elif dtype == "bool":
        return "boolean", None
    elif "datetime" in dtype:
        return "char", "*"
    else:
        return "char", "*"


def table_to_pandas_safe(table):
    """
    Convert Astropy Table to Pandas DataFrame.
    Multidimensional columns become list-of-lists.
    """

    data_dict = {}

    for col in table.colnames:
        data = table[col]

        if len(data.shape) > 1:
            # 2D column → list per row
            data_dict[col] = [row.tolist() for row in data]
        else:
            # 1D scalar column
            data_dict[col] = data.data  # numpy array
        if isinstance(data_dict[col], np.ndarray) and not data_dict[col].dtype.isnative:
            data_dict[col] = data_dict[col].byteswap().newbyteorder()

    return pd.DataFrame(data_dict)