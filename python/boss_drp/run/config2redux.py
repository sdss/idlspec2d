from jinja2 import Template
from boss_drp.field import field_to_string, Fieldtype, Field
from boss_drp.Config import config
import boss_drp
import copy
import os
import os.path as ptt
import datetime


def idl_flags(key, value=None):
    if value is None:
        return f'/{str(key).strip()}'
    else:
        return f'{str(key).strip()}={str(value).strip()}'

def py_flags(key, value=None):
    if value is None:
        return f'--{str(key).strip()}'
    else:
        return f'--{str(key).strip()} {str(value).strip()}'


def config2redux(plan2d = [], plancombine = '', 
                 obs=None, field=None, mjd=None, epoch=False, 
                 custom=None, allsky=False, custom_single_mjd=None):
    pipe_flags = {}
    if (not epoch) and (not custom):
        daily= True
    else:
        daily= False

    
    cfield = [field_to_string(field)]
    if custom:
        _field = field
        cfield = [f'{custom}']
        if obs:
            cfield.append(obs)
        spfieldfile = f"spFullsky-{'_'.join(cfield)}-{mjd}.fits"
        if custom_single_mjd:
            scriptname = f"redux_{'_'.join(cfield)}-{mjd}_{custom_single_mjd}"
        else:
            scriptname = f"redux_{'_'.join(cfield)}-{mjd}"
        field = '_'.join(cfield)
    elif epoch:
        _field = int(field)
        field = field_to_string(field)
        spfieldfile = f'spField-{field}-{mjd}.fits'
        scriptname = f"redux-{field}-{mjd}"

    else:
        _field = int(field)
        field = field_to_string(field)
        spfieldfile = f'spField-{field}-{mjd}.fits'
        scriptname = f"redux-{field}-{mjd}"

    field_dir = Field(config.pipe.get('general.BOSS_SPECTRO_REDUX'),
                      config.pipe.get('general.RUN2D'), field, 
                      mjd=mjd,custom=(custom is not None), custom_name=custom,
                      epoch = epoch, run1d = config.pipe.get('general.RUN1D'), obs = obs)
    
    
    if not config.pipe['Clobber.clobber_pipe']:
        if ptt.exists(ptt.join(field_dir.dir(),scriptname)):
            return None
    if config.queue.get('no_write'):
        return ptt.join(field_dir.dir(),scriptname)
        
    topdir = config.pipe.get('general.BOSS_SPECTRO_REDUX')
    if topdir == os.getenv('BOSS_SPECTRO_REDUX',''):
        topdir = None

    pipe_flags = dict(date = datetime.datetime.now().strftime("%c"),
                      run2d = config.pipe.get('general.RUN2D'),
                      run1d =config.pipe.get('general.RUN1D'),
                      field = field,
                      mjd   = (mjd if not custom else custom_single_mjd), 
                      planmjd = mjd,
                      obs   = obs,
                      BOSS_SPECTRO_DATA = None,
                      GCAM_DATA = None,
                      topdir = topdir,
                      field_dir = field_dir.dir(),
                      spfieldfile = spfieldfile,
                      daily=daily,
                      epoch=epoch,
                      custom=custom,
                      fibermap = (config.pipe.get('Stage.run_fibermap') and daily),
                      reduce2d = (config.pipe.get('Stage.run_reduce2d') and daily),
                      combine = config.pipe.get('Stage.run_combine'),
                      customcombine = config.pipe.get('Stage.run_combine'),
                      analyze = config.pipe.get('Stage.run_analyze'),
                      XCSAO = config.pipe.get('Stage.run_XCSAO'),
                      fieldlist= (config.pipe.get('Stage.run_fieldlist') and (not custom)),
                      fieldmerge= (config.pipe.get('Stage.run_fieldmerge')),
                      specFiles = (config.pipe.get('Stage.run_specFiles')),
                      spcalib = (config.pipe.get('Stage.run_spcalib') and (not custom)),
                      healpix = (config.pipe.get('Stage.run_healpix') and (not custom) and (not epoch)),
                      spplan2d = plan2d, #Lists
                      spplancomb = ptt.basename(plancombine),
                      idlspec2d_module_change = None,
                      idlutils_module_change = None
    )
    if obs is not None:
        lco = True if obs.lower() == 'lco' else False
    else: 
        lco = False
    if daily:
        pipe_flags['BOSS_SPECTRO_DATA'] = 'BOSS_SPECTRO_DATA_S' if lco else 'BOSS_SPECTRO_DATA_N'
        pipe_flags['GCAM_DATA'] = 'GCAM_DATA_S' if lco else 'GCAM_DATA_N'
        

    if config.pipe.get('general.RUN2D') != config.pipe.get('general.RUN1D'): 
        run1d = config.pipe.get('general.RUN1D')
        pipe_flags['idlspec2d_module_change'] = f'module switch idlspec2d idlspec2d/{run1d.strip()}'
    if config.pipe.get('general.idlutils_1d') is not None:
        idlutils_1d = config.pipe.get('general.idlutils_1d')
        pipe_flags['idlutils_module_change'] = f'module switch idlutils idlutils/{idlutils_1d.strip()}'

            

    if not custom:
        ftype = Fieldtype(fieldid=field, mjd=mjd)
        plates = ftype.plates
        legacy=ftype.legacy
    else:
        plates = False
        legacy = False
    

    extra_lookup = dict(release = config.pipe.get('general.RELEASE'),
                      V_TARG = None if config.pipe.get('general.V_TARG') == '*' else '',
                      remote = config.pipe.get('general.REMOTE'),
                      lco = True if obs == 'lco' else False,
                      plates = plates, legacy=legacy,
                      epoch = epoch, custom=custom, allsky=allsky,
                      custom_field=('_'.join(cfield) if custom is not None else None))
    
    _config = copy.deepcopy(config.pipe)


    if lco and int(mjd) < 60402: #Pre proper implementation of fixed-screen flats for arc2trace
        lco_rm_fps = [23129,23130,23131,23132,23133,23134,23135,23136,23137,23175,
                        23288,23408,23409,23410,31687,31688,112357,112358,112362]
        if _field not in lco_rm_fps:
            _config.reduce = True

                               #conversion, config, extra, config_exclude
    stage = dict(fibermap=      (py_flags,  _config.get('fibermap'),
                                 ['release','remote','V_TARG'],[]),
                 reduce2d=        (idl_flags, _config.get('reduce'),  
                                 ['lco', 'plates', 'legacy'],[]),
                 combine=       (idl_flags, _config.get('combine'), 
                                 ['legacy','plates','epoch'],[]),
                 customcombine= (idl_flags, _config.get('combine.custom_coadd'),
                                 ['epoch'],[]),
                 analyze=       (idl_flags, _config.get('analyze.spec1d'),
                                 {'custom':'custom_field','epoch':'epoch'},[]),
                 XCSAO=         (py_flags,  _config.get('analyze.XCSAO'),
                                 {'custom':'custom_field','epoch':'epoch'},[]),
                 fieldlist=     (py_flags,  _config.get('post.fieldlist'),
                                 ['epoch'],[]),
                 fieldmerge=    (py_flags,  _config.get('post.fieldmerge'),
                                 ['custom','epoch'],['skip_specprimary']),
                 specFiles=     (py_flags,  _config.get('post.spec'),
                                 ['epoch','custom','allsky'],[]),
                 spcalib=       (py_flags,  _config.get('post.spCalib'),
                                 ['epoch'],[])
                )                 

    #CLOBBER flag


    for cmd, (convert, _config, extras, exclude) in stage.items():
        flags = []
        for key, value in _config.items():
            if key in exclude:
                continue
            if value is None:
                continue
            if isinstance(value, bool) or str(value).lower() in ['true', 'false']:
                if str(value).lower() == 'true':
                    flags.append(convert(key))
                continue
            if isinstance(value, dict):
                continue
            if key == 'skip_specprimary':
                if value == 'update':
                    flags.append(convert('update_specprimary'))
                    continue
            flags.append(convert(key, value=value))

        if (pipe_flags[cmd]) and (epoch):
            flags.append(convert('epoch'))

        if isinstance(extras, dict):
            for key1, key2 in extras.items():
                val2 = extra_lookup[key2]
                if val2 is None:
                    continue
                if isinstance(val2, bool) or str(val2).lower() in ['true', 'false']:
                    if str(val2).lower() == 'true':
                        flags.append(convert(key1))
                    continue
                if isinstance(val2, dict):
                    continue
                flags.append(convert(key1,val2))
        else:
            for key1 in extras:
                val2 = extra_lookup[key1]
                if val2 is None:
                    continue
                if isinstance(val2, bool) or str(val2).lower() in ['true', 'false']:
                    if str(val2).lower() == 'true':
                        flags.append(convert(key1))
                    continue
                if isinstance(val2, dict):
                    continue
                flags.append(convert(key1, extra_lookup[key1]))

        if '--' in convert('test'):
            flags = ' '.join(flags) 
        else:
            flags = ', '.join(flags) 
            if len(flags) > 0: flags = flags+', '
        pipe_flags[f'{cmd}flags'] = flags



    template = ptt.join(ptt.dirname(boss_drp.__file__), 'etc','templates','redux.j2')
    print(ptt.join(field_dir.dir(),scriptname))
    with open(ptt.join(field_dir.dir(),scriptname), "w", encoding="utf-8") as output_file:
        with open(template) as template_file:
            j2_template = Template(template_file.read())
            output_file.write(j2_template.render(pipe_flags))
    return ptt.join(field_dir.dir(),scriptname)


