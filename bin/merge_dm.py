from splog import Splog
#splog = Splog()

########################################
from pydl.pydlutils.yanny import read_table_yanny, yanny
from astropy.table import Table, vstack, join, Column, MaskedColumn
import numpy as np
from astropy.io import fits


def tableToModel(table, dm_ext, name, splog, old=False, drop_cols=None, verbose=False):
    if drop_cols is not None:
        drop_cols = np.atleast_1d(drop_cols)
    #dm_table = Table()
    cols = table.colnames
    for col in cols:

        if col.upper() not in dm_ext['Column'].data:
            if drop_cols is not None:
                if col in drop_cols:
                    table.remove_column(col)
                    continue
            if verbose:
                splog.log(col+' missing from datamodel for '+name)
            table.remove_column(col)
            continue
        if 'K' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            if 'uK' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
                dtype = np.uint64
            else:
                dtype = int
        elif 'J' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            dtype = np.int32
        elif 'I' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            dtype = np.int16
        elif 'D'  in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            dtype = float
        elif 'E'  in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            dtype = np.float32
        elif 'L' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            dtype = bool
            if table[col].dtype.type == np.str_:
                test = np.zeros(len(table[col].data))
                test[np.where(table[col].data == 'F')] = 0
                test[np.where(table[col].data == 'T')] = 1
                test=test.astype(bool)
            else:
                test = table[col].astype(bool).data
        elif 'B' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
            dtype = 'uint8'
        else:
            dtype = object
        
        shape = ''.join(c for c in dm_ext[dm_ext['Column'] == col.upper()]['type'][0] if c.isdigit())

        if old:
            data = table[col].data
            if dtype == object:
                if sum(data == np.asarray(['']*len(data))) != 0:
                    continue
        data = table[col].data
        table.remove_column(col)
        if shape == '':
            if dtype == bool:
                try:
                    table.add_column(Column(test, name = col.upper()))
                except:
                    splog.log(col)
                    table.add_column(Column(test, name = col.upper()))
            elif dtype == 'uint8':
                try:
                    data[data.astype(object) == ''] = 0
                    table.add_column(Column(data.astype(dtype), name = col.upper()))
                except:
                    splog.log(col)
                    table.add_column(Column(data.astype(dtype), name = col.upper()))
            else:
                try:
                    table.add_column(Column(data.astype(dtype), name = col.upper()))
                except:
                    splog.log(col)
                    table.add_column(Column(data.astype(dtype), name = col.upper()))
        else:
            table.add_column(Column(data.astype(dtype), name = col.upper(), shape=(shape,)))
        data = None
#        table.remove_column(col)
    return(table)

def merge_dm(table=None, ext = 'Primary', name = None, hdr = None, dm ='spfibermap_dm.par',
             old_tab = None, splog=Splog(), drop_cols = None, verbose=False):

    model=read_table_yanny(dm,'MODEL')
    model.convert_bytestring_to_unicode()
    dm_model = model[model['Name'] == ext][0]
    
    if dm_model['ext'] != 'None' :
        dm_ext = read_table_yanny(dm, dm_model['ext'])
        dm_ext.convert_bytestring_to_unicode()
        if table is not None:
            table = tableToModel(table, dm_ext, name, splog, old=False, drop_cols=drop_cols, verbose=verbose)
###############################
            for col in table.colnames:
                if (table[col].dtype == int) or (table[col].dtype == np.int16) or (table[col].dtype == np.int32):
                    coldat = table[col].data
                    try:
                        fill = int(dm_ext[dm_ext['Column'] == col.upper()]['null'][0])
                    except:
                        if verbose:
                            splog.log('WARNING: datamodel ('+dm+') is missing fill value for '+col)
                        fill = -999
                        try:
                            dm_ext[dm_ext['Column'] == col.upper()]['null'][0] = -999
                        except:
                            pass
                    coldat[np.where(coldat == 999999)[0]] = fill
                    coldat = np.ma.masked_values(coldat, fill)
                    table[col] = MaskedColumn(coldat,fill_value = fill)
###############################
            if old_tab is not None:
                old_tab = tableToModel(old_tab, dm_ext, name, splog, old=True, drop_cols=drop_cols, verbose=verbose)

                for col in old_tab.colnames:
                    if (old_tab[col].dtype == int) or (old_tab[col].dtype == np.int16) or (old_tab[col].dtype == np.int32):
                        coldat = old_tab[col].data
                        fill = int(dm_ext[dm_ext['Column'] == col.upper()]['null'][0])
                        coldat[np.where(coldat == 999999)[0]] = fill
                        coldat = np.ma.masked_values(coldat, fill)
                        old_tab[col] = MaskedColumn(coldat,fill_value = fill)

                    if 'B' in dm_ext[dm_ext['Column'] == col.upper()]['type'][0]:
                        if old_tab[col].shape[1] > table[col].shape[1]:
                            coldat = table[col].data
                            pad = old_tab[col].shape[1] - table[col].shape[1]
                            table[col] = np.pad(coldat, [(0,0),(pad,0)], mode = 'constant', constant_values= 0)
                        elif old_tab[col].shape[1] < table[col].shape[1]:
                            coldat = table[col].data
                            pad = table[col].shape[1] - old_tab[col].shape[1]
                            old_tab[col] = np.pad(coldat, [(0,0),(pad,0)], mode = 'constant', constant_values= 0)
                    if (old_tab[col].dtype.kind in ['U','S']):
                        try:
                            if table[col].dtype == object:
                                old_tab[col] = old_tab[col].astype(object)
                        except:
                            pass
                table = vstack([old_tab, table])
###############################
            for col in table.colnames:
                if (table[col].dtype == object) or ('|S' in str(table[col].dtype)):
                    try:
                        maxlen = len(max(table[col].data, key=len))
                        if maxlen == 0:
                            maxlen = 20
                    except:
                        maxlen = 20
                    table[col] = table[col].astype('<S'+str(maxlen))


        dm_ext = read_table_yanny(dm, dm_model['ext'])
        dm_ext.convert_bytestring_to_unicode()
        cols = []
        for row in dm_ext:
            try:
                null = int(row['null']) if row['null'] != '' else None
            except:
                null = row['null']
            if table is None:
                data = None
            elif row['Column'] not in table.colnames:
                data = None
            else:
                data = table[row['Column']].data
                shape = table[row['Column']].dtype
                table.remove_column(row['Column'])
            if row['type'] != 'A':
                if row['type'] != 'B':
                    print(row['type'],row['Column'])
                    if row['type'] == 'uK':
                        bzero=fits.BinTableHDU(Table([data]),uint = True).columns[0].bzero
                        cols.append(fits.Column(name = row['Column'], format = row['type'][1:], null=null, array= data, bzero=bzero))
                    else:
                        cols.append(fits.Column(name = row['Column'], format = row['type'], null=null, array= data, ))
                else:
                    if data is not None:
                        N, F = data.shape
                    else:
                        F = 0
                    cols.append(fits.Column(name = row['Column'], format=f"{F}B", dim=f"({F})", array=data, ))
            else:
                if data is None:
                    shape = '10'
                else:
                    shape = str(shape).replace('|S','').replace('<U','').replace('<S','')

                cols.append(fits.Column(name = row['Column'], format = shape+'A', array= data))
            shape = None
            data = None
        hdu = fits.BinTableHDU.from_columns(cols, name = name, uint = True)
        
        for card in hdu.header.cards:
            match = np.where(str(card[1]) == dm_ext['Column'].data)[0]
            if len(match) > 0:
                comment = dm_ext[match[0]]['description']
                hdu.header.set(card[0],card[1],comment)
    else:
        hdu = fits.PrimaryHDU()
    if dm_model['hdr'] != 'None':
        dm_hdr = read_table_yanny(dm, dm_model['hdr'])
        dm_hdr.convert_bytestring_to_unicode()
        for row in dm_hdr:
            val = '' if row['card'] not in hdr else hdr[row['card']]
            hdu.header.set(row['card'], val, row['description'])
    return(hdu)
