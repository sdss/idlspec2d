import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits

def normalize_array_columns_for_vstack(tables):
    """
    Make all array columns have the same shape across tables by padding.
    """
    colnames = tables[0].colnames

    for col in colnames:
        # Only worry about 2D array columns
        if len(tables[0][col].shape) < 2:
            continue

        # Find max second dimension
        max_len = max(t[col].shape[1] for t in tables)

        for i, t in enumerate(tables):
            if t[col].shape[1] != max_len:
                pad_width = max_len - t[col].shape[1]
                if col.upper() == "SDSS5_TARGET_FLAGS":
                    t[col] = np.pad(t[col], [(0,0),(0,pad_width)], mode="constant", constant_values=0)
                else:
                    t[col] = np.pad(t[col], [(0,0),(pad_width,0)], mode="constant", constant_values=0)
            t[col] = np.array([row.tolist() for row in t[col]], dtype=object)
    return(tables)

def build_table_from_fits(fits_files):
    tables = []

    for f in fits_files:
        t = Table.read(f)  # adjust HDU if needed
        t.meta = {} # clearing metadata to avoid vstack issues
        tables.append(t)
    tables = normalize_array_columns_for_vstack(tables)

    full_table = vstack(tables)

    # Normalize string columns (avoid fixed-width unicode issues)
    for col in full_table.colnames:
        if full_table[col].dtype.kind in ["U", "S"]:
            full_table[col] = full_table[col].astype(str)

    hdr = fits.getheader(fits_files[-1])
    return full_table, hdr