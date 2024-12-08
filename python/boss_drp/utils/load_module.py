from os import getenv



def load_env(key, default=None):
    val = getenv(key)
    if val is None:
        val = default
        if val is None:
            print('ERROR: '+key+' is not set')
            exit()
    return(val)
