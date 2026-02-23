#!/usr/bin/env python3
import numpy as np
import warnings
import re
import pandas as pd

coaddids = pd.DataFrame([
                        {'name':'daily','id':'00'}, #boss daily coadds
                        {'name':'epoch','id':'01'}, #boss field-epoch coadds
                        {'name':'allepoch','id':'02'},  #boss allepoch coadds
                        {'name':'spiders','id':'02'},#DR18 boss allepoch coadd alias name
                        {'name':'efeds','id':'02'},#DR18 boss allepoch coadd alias name
                        {'name':'allepoch_apo','id':'03'},  #boss allepoch (apo) coadds
                        {'name':'allepoch_lco','id':'04'},  #boss allepoch (lco) coadds
                        {'name':'allvisit','id':'10'},  #apogee allvist
                        {'name':'allstar','id':'11'}, #apogee allstar
                        {'name':'test','id':'99'}
                        ])
#
#coaddids = {'daily':'00', #boss daily coadds
#            'epoch':'01', #boss field-epoch coadds
#            'allepoch':'02', #boss allepoch coadds
#            'spiders':'02', #DR18 boss allepoch coadd alias name
#            'allvisit':'10', #apogee allvist
#            'allstar':'11', #apogee allstar
#            'test':'99'
#            }
#
#coaddids_inv = {v: k for k, v in coaddids.items()}


class UndefinedCoadd(Exception):
    """Exception raise for coadds not in coaddids"""
    def __init__(self, coadd, message = None):
        if message is None:
            message=f"Coadd name '{coadd}' not in {','.join(coaddids.name)}"
        self.message = message
        self.coadd = coadd
        super().__init__(message)

class MissingFiberIDs(Exception):
    """Exception raise for MissingFiberIDs not supplied to calculate legacy specobjids"""
    def __init__(self, message = None):
        if message is None:
            message=f"Fiber IDs are not supplied and are required to calculate legacy specobjids"
        self.message = message
        super().__init__(message)

lsh = lambda x,s: np.uint64(x)*np.uint64(2**s)

def encode_legacy(fieldid, mjd, fiberid, tag):
    try:
        tag = str(tag.strip())
        if '_' in tag:
            n,m,p = tag.split('_')
            n = int(n[1:])
            m = int(m)
            p = int(p)
            run2d = (n-5)*10000 + m*100 + p
        elif tag.isnumeric():
            run2d = int(tag)
        else:
            print(f"WARNING: Unable to parse RERUN from {tag} for CAS-style SPECOBJID; Using 0 instead")
            run2d = 0
    except Exception as e:
        print(f"WARNING: Unable to parse RERUN from {tag} for CAS-style SPECOBJID; Using 0 instead")
        run2d = 0
    specobjid = np.uint64(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        specobjid |= lsh(fieldid,50) | lsh(fiberid,38) | lsh(int(mjd)-50000,24) | lsh(run2d,10)
    return(specobjid)

def get_coadd(coadd):
    coadd = str(coadd).lower().strip()
    if not coadd.isnumeric():
        try:
            coadd = coaddids[coaddids.name == coadd].id.iloc[0] #coaddids[coadd]
        except:
            raise(UndefinedCoadd(coadd))
    elif len(str(coadd)) > 2:
        raise(UndefinedCoadd(coadd, message = f'Coadd ID ({coadd}) must be less then 100'))
    return coadd

def encode_tag(tag):
    if ('v' in tag) and (len(tag.split('_')) == 3) :
        n,m,p = tag.strip().split('_')
        n = re.sub(r'\D', '', n)
        p = p.split('-')[0]
        p = re.sub(r'\D', '', p)
        tag = n.zfill(2)+m.zfill(2)+p.zfill(2)
    elif ('v' in tag) and (len(tag.split('.')) == 3) :
        n,m,p = tag.strip().split('.')
        n = re.sub(r'\D', '', n)
        p = p.split('-')[0]
        p = re.sub(r'\D', '', p)
        tag = n.zfill(2)+m.zfill(2)+p.zfill(2)
    elif len(tag.split('.')) == 2:
        n,m = tag.strip().split('.')
        tag = n.zfill(2)+m.zfill(2)+'0'.zfill(2)
    elif len(tag.split('.')) == 3:
        n,m,p = tag.strip().split('.')
        p = p.split('-')[0]
        p = re.sub(r'\D', '', p)
        tag = n.zfill(2)+m.zfill(2)+p.zfill(2)
    elif 'dr' in tag.lower():
        tag = tag.replace('DR','')
        tag = tag.zfill(6)
    elif tag.isnumeric():
        tag = tag.zfill(6)
    else:
        print(f"WARNING: Unable to parse Tag from {tag[0]} for SDSS-V SPECOBJID; Using 000000 instead")
        tag = '000000'
    return tag
    
def encode_V(sdssid, fieldid, mjd, coadd, tag):

    if isinstance(coadd, (list, np.ndarray)):
        if len(set(coadd)) > 1:
            coadd = np.array([get_coadd(x) for x in coadd])
        else:
            coadd = get_coadd(coadd[0])
    else:
        coadd = get_coadd(coadd)
  
    if True:
    #try:
        if isinstance(tag, (list, np.ndarray)):
            if len(set(tag)) == 1:
                tag = encode_tag(tag[0])
            else:
                tag = np.array([get_coadd(x) for x in tag])
        else:
            tag = encode_tag(tag)
    #except:
    #    print(f"WARNING: Unable to parse Tag from {tag[0]} for SDSS-V SPECOBJID; Using 000000 instead")
    #    tag = '000000'

    if type(coadd) is str:
        coadd = np.atleast_1d(np.char.add(coadd.zfill(2), tag))
        if len(coadd) == 1:
            coadd = np.asarray([coadd[0]]*len(sdssid))

    else:
        coadd = np.char.add(np.char.zfill(coadd,2), tag)

    mask = sdssid != ''
    maxlen = np.char.str_len(sdssid).max() + 20
    result = np.full_like(sdssid, '', dtype=f'<U{maxlen}')
    
    result[mask] = np.char.add( np.char.add(np.char.add(sdssid[mask], fieldid[mask]),
                                                        mjd[mask]),coadd[mask])
    return(result)
    
def encode(sdssid, fieldid, mjd, coadd, tag, fiberid=None, allnew = False):
    """ Defines the SDSS-V version of SpecobjIDs

    This function accepts a single sdss_id or list of sdss_ids
    with their corresponding fieldid, mjd, coadd version, and
    pipeline tag. If the sdss_id is a single value, the fieldid
    and mjd parameters are assumed to be a single value, and if
    they are not, then the first value of the array is used. If
    the sdss_id is an array, then the fieldid and mjd are assumed
    to be a single set of values (in which case the values are
    used for all sdss_id) or an array of the same length as the
    sdss_id.

    Parameters:
        sdssid, fieldid, mjd (float or string or array of floats or strings):
            the sdss_id and its corresponding field and mjds
        coadd (string or int):
            the name of the coadd or the numeric id of the coadd
            as defined by coaddids
        tag (string):
            the pipeline tag used
        fiberid (optional, float or string or array of floats or strings)
            the fiber ids used for legacy style specobjids
        allnew (optional, boolean)
            calculates SDSS-V style specobjids for all targets

    Returns:
        SpecobjIDs (string or array of string matching length of sdssid):
            the SpecobjIDs
    Example:
        >>> encode(1,16000,59999,'daily','v6_1_3')
        >>> '100160005999900060103'
        >>> encode([1,2],16000,59999,'daily','v6_1_3')
        >>> ['100160005999900060103', '200160005999900060103']
    """
  
    arrsdssid = False
    if isinstance(sdssid, list):
        arrsdssid = 'list'
    elif isinstance(sdssid, np.ndarray):
        arrsdssid = 'arr'
    
    sdssid = np.atleast_1d(sdssid).astype(str)
    try:
        fieldid = np.atleast_1d(fieldid)
        if not isinstance(fieldid, np.ma.MaskedArray):
            fieldid = np.ma.masked_array(fieldid)
        fieldid = np.char.zfill(fieldid.filled(0).astype(str),7)
    except:
        fieldid = np.char.zfill(np.atleast_1d(fieldid).astype(str),7)
    if len(fieldid) == 1:
        fieldid = np.asarray([fieldid[0]]*len(sdssid))
  
    mjd = np.atleast_1d(mjd).astype(str)
    if len(mjd) == 1:
        mjd = np.asarray([mjd[0]]*len(sdssid))

    if fiberid is not None:
        fiberid = np.atleast_1d(fiberid)
        if len(fiberid) == 1:
            fiberid = np.asarray([fiberid[0]]*len(sdssid))

    legacy = False
    tagt = np.atleast_1d(tag)
    if ('v' in tagt[0]) and (len(tagt[0].split('_')) == 3):
        n,m,p = tagt[0].strip().split('_')
        p = p.split('-')[0]
        p = re.sub(r'\D', '', p)
        n = re.sub(r'\D', '', n)
        if int(n) < 6:
            legacy = True
        elif int(p) <= 4 and int(n) == 6 and int(m) == 0:
            legacy = True
    elif ('v' in tagt[0]) and (len(tagt[0].split('.')) == 3):
        n,m,p = tagt[0].strip().split('.')
        p = p.split('-')[0]
        p = re.sub(r'\D', '', p)
        n = re.sub(r'\D', '', n)
        if int(n) < 6:
            legacy = True
        elif int(p) <= 4 and int(n) == 6 and int(m) == 0:
            legacy = True
    elif tagt[0].isnumeric():
        if int(tagt) in [103,104,26]:
            legacy = True
    if legacy and not allnew:
        specojbid = []
        
        if fiberid is not None:
            fiberid = np.atleast_1d(fiberid)
            if len(fiberid) == 1:
                fiberid = np.asarray([fiberid[0]]*len(sdssid))
        else:
            raise MissingFiberIDs()
            
        for pid,mjd,fid in zip(fieldid, mjd, fiberid):
            specojbid.extend(encode_legacy(fieldid, mjd, fiberid, tag))
        specojbid = np.asarray(specojbid).astype(str)
    else:
        specojbid= encode_V(sdssid, fieldid, mjd, coadd, tag)
    
    bsid = np.where(sdssid == '-999')[0]
    if len(bsid) >0:
        specojbid[bsid] = ''
    
    if not arrsdssid:
        return(specojbid[0])
    elif arrsdssid == 'list':
        return(specojbid.tolist())
    return(specojbid)



def decode_V(specobjid):
    """
    Decode SDSS-V (DR19+) SpecobjID
    
    Parameters:
        specobjid (float or string)
    
    Returns:
        unwrap (dictionary)
            dictionary of the attributes used to build specobjid
    """
    unwrap = {}
    unwrap['specobjid'] = str(specobjid)
    if specobjid == '':
        return unwrap
    specobjid = str(specobjid)
    try:
        unwrap['fieldid'] = int(specobjid[-20:-13])
    except:
        unwrap['fieldid'] = specobjid[-20:-13]
    unwrap['mjd'] = int(specobjid[-13:-8])
    unwrap['sdss_id'] = int(specobjid[:-20])

    try:
        unwrap['coadd'] = coaddids[coaddids.id == specobjid[-8:-6]].name.iloc[0]
        #coaddids_inv[specobjid[-8:-6]]
    except:
        unwrap['coadd'] = specobjid[-8:-6]

    if int(specobjid[-6:]) in [104,103,26]:
        unwrap['tag'] = str(int(specobjid[-6:]))
        unwrap['instrument'] = 'sdss'
    elif int(specobjid[-6:]) in [17]:
        unwrap['tag'] = 'DR17'
        unwrap['instrument'] = 'apogee'
    elif unwrap['coadd'] in ['allstar','allvisit']:
        unwrap['tag'] = '{0:d}.{1:d}.{2:d}'.format(int(specobjid[-6:-4]),
                                                   int(specobjid[-4:-2]),
                                                   int(specobjid[-2:]))
        unwrap['instrument'] = 'apogee'
    elif int(specobjid[-6:]) == 0:
        unwrap['tag'] = 'unknown'
        unwrap['instrument'] = 'unknown'
    elif int(specobjid[-6:-4]) == 0:
        unwrap['tag'] = str(int(specobjid[-6:]))
        unwrap['instrument'] = 'unknown'
    else:
        unwrap['tag'] = 'v{0:d}_{1:d}_{2:d}'.format(int(specobjid[-6:-4]),
                                                    int(specobjid[-4:-2]),
                                                    int(specobjid[-2:]))
        unwrap['instrument'] = 'boss'
    return(unwrap)
  
def decode_legacy(specobjid):
    """
    Decode pre-SDSS-V (and SDSS-V DR18) SpecobjID
    
    Parameters:
        specobjid (float or string)
    
    Returns:
        unwrap (dictionary)
            dictionary of the attributes used to build specobjid
    """
    unwrap = {}
    unwrap['specobjid'] = specobjid
    specobjid = np.atleast_1d(specobjid).astype(np.uint64)
    unwrap['fieldid'] = np.bitwise_and(specobjid >> 50, 2**17 - 1)[0]
    unwrap['fiber'] = np.bitwise_and(specobjid >> 38, 2**12 - 1)[0]
    unwrap['mjd'] = int(np.bitwise_and(specobjid >> 24, 2**14 - 1)[0] + 50000)
    unwrap['tag'] = np.bitwise_and(specobjid >> 10, 2**14 - 1)
    N = ((unwrap['tag'] // 10000) + 5).tolist()
    M = ((unwrap['tag'] % 10000) // 100).tolist()
    P = (unwrap['tag'] % 100).tolist()

    unwrap['tag'] = [str(p)
                    if (n == 5 and m == 0) else 'v{0:d}_{1:d}_{2:d}'.format(n, m, p)
                    for n, m, p in zip(N, M, P) ][0]
    unwrap['coadd'] = 'epoch'
    unwrap['instrument'] = 'boss'
    return(unwrap)


def decode(specobjid):
    """
    Decodes SpecobjIDs (whether the old or new style)

    Parameters:
        specobjid (float or string or array of floats or strings):
            the SpecobjIDs to decode (can be either single value, list, or numpy array)
      
    Returns:
        attrib_uid (dictionary, list of dictionaries, or array of dictionaries)
            dictionary of the attributes used to build SpecobjIDs matches SpecobjIDs dimentions
    """
  
    arrsid = False
    if isinstance(specobjid, list):
        arrsid = 'list'
    elif isinstance(specobjid, np.ndarray):
        arrsid = 'arr'
    specjobid = np.atleast_1d(specobjid).astype(str)
    attrib_sid = [decode_V(x) if (len(x) > 20 or len(x) == 0) else decode_legacy(x) for x in specjobid]

    if not arrsid:
        return(attrib_sid[0])
    elif arrsid == 'arr':
        return(np.asarray(attrib_sid))
    return(attrib_sid)
  

if __name__ == '__main__' :
    from specobjid import encode, decode, coaddids
    #encode(sdssid, fieldid, mjd, coadd, tag, fiberid=None):

    print(encode(1,16000,59999,'daily','v6_1_3',1))
    print(decode(100160005999900060103))
    print('  ')
    print(encode(1,16000,59999,'daily','v6_0_3',1))
    print(decode(18014398952125516800))
    print('  ')
    print(encode(1,889,52663,'epoch','26',1))
    print(decode(1000925336738752512))
    print('  ')
    print(encode(1,8954,57453,'epoch','v5_13_2',1))
    print(decode(10081308165788686336))

