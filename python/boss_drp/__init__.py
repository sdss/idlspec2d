#import subprocess
#import os
#from pkg_resources import resource_filename

from sdsstools import get_package_version
#from sdsstools import get_config, Configuration
import os
import numpy as np

__version__ = get_package_version(__file__, 'boss_drp') or 'dev'

#def _config_None2null(config):
#    for key1 in config.keys():
#        if type(config[key1]) is Configuration:
#            config[key1] = _config_None2null(config[key1])
#        elif type(config[key1]) is str:
#            if config[key1].lower() == 'none':
#                config[key1] = None
#            else:
#                pass
#        else:
#            pass
#    return(config)
#    
#config = get_config('boss_drp', allow_user=True,
#                    default_envvars={'BOSS_SPECTRO_REDUX':os.getenv('BOSS_SPECTRO_REDUX'),
#                                     'RUN2D':os.getenv('RUN2D'),
#                                     'RUN1D':os.getenv('RUN1D'),
#                                     'DAILY_DIR':os.getenv('DAILY_DIR'),
#                                     'OBSERVATORY':os.getenv('OBSERVATORY')})
#config = _config_None2null(config)
#
#
#
#def config_to_args(config, args, cmd, slurm=False, load=True, update_config=False):
#    cmd = np.atleast_1d(cmd).tolist()
#    cmdl = np.char.lower(np.atleast_1d(cmd)).tolist()
#    def_config = config.copy()
#    vargs = args if type(args) is dict else vars(args)
#    if load:
#        try:
#            config.load(args.config, use_base=True)
#            config = _config_None2null(config)
#        except:
#            pass
#    for ckey in ['general']+cmd+cmdl:
#        try:
#            step_config = config[ckey]
#        except:
#            try:
#                if '.' in key:
#                    step_config = config[key.split('.')[0]][key.split('.')[1]]
#            except:
#                continue
#        if step_config is None:
#            continue
#        for skey in step_config.keys():
#            if skey not in vargs.keys():
#                continue
#            val = vargs[skey]
#            if val in [None,False]:
#                args.__dict__[skey] = step_config.get(skey)
#
#        if update_config:
#            for skey in vargs.keys():
#                if skey not in step_config.keys():
#                    continue
#                step_config[skey] = vargs[skey]
#        if slurm:
#            try:
#                slurm_config = step_config['slurm']
#            except:
#                continue
#            for skey in slurm_config.keys():
#                if skey in vargs.keys():
#                    val = vargs[skey]
#                    if val in [None,False]:
#                        args.__dict__[skey] = slurm_config.get(skey)
#            
#    if not update_config:
#        return(args)
#    return(args, config)
##config_to_args(config, args, 'spTrace')
##    for ckey in config.keys():
##        if ckey.lower() not in ['general']+cmd:
##            continue
##        step_config = config.get(ckey)
##        for skey in in step_config.keys():
##            if skey in vars(args).keys():
##                val = vars(args)[skey]
##                if val in [None,False]:
##                    args.__dict__[skey] = step_config.get(skey)
##        if slurm:
##            for skey in step_config['slurm'].keys():
##                if skey in vars(args).keys():
##                    val = vars(args)[skey]
##                    if val in [None,False]:
##                        args.__dict__[skey] = step_config['slurm'].get(skey)
#
# #   return(args)



#__version__ = subprocess.getoutput("idlspec2d_version")
#if 'command not found' in __version__:
    #my_env = os.environ.copy()
    #my_env["PATH"] = f"{os.path.abspath(resource_filename('boss_drp','../../bin/'))}:{my_env['PATH']}"
    #my_env['IDLSPEC2D_DIR'] = os.path.abspath(resource_filename('boss_drp','../../'))
    #p = subprocess.Popen("idlspec2d_version", env=my_env, stdout=subprocess.PIPE)
    #out, err = p.communicate()
    #__version__ = out.decode().replace('\n','')
