def getcard(hdr, card, default=None, noNaN=False):
    try:
        if hdr.count(card) > 0:
            if type(hdr[card]) is str:
                hdr[card] = hdr[card].strip().replace("'","")
                if noNaN is True:
                    if hdr[card].strip().upper() == 'NAN':
                        hdr[card] = default
                return(hdr[card])
            elif hdr[card] is None:
                return(default)
            else:
                return(hdr[card])
        else:
            return(default)
    except:
        return(default)
