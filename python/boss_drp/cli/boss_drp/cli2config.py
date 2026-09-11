from boss_drp.field.generations import generations
from boss_drp.Config import config, update_key

def cli2config(args, config_par = None, exclude = None, set_gen = True):
    config_par = config_par or {}
    exclude = exclude or []
    config.load(args.queue_config, queue_config_file= args.queue_config_file,
                 config_name=args.config, config_file=args.config_file)
    for name, val in config_par.items():
        config.queue.set(name, val)

    exclude_args = ['config','config_file','queue_config','queue_config_file','show_config',
                    'mjdstart','mjdend','fieldstart','fieldend']#
    exclude_args.extend(exclude)

    new_mjdrange = None
    if args.get('mjdstart') is not None and args.get('mjdend') is not None:
        new_mjdrange = [args.get('mjdstart'), args.get('mjdend')]
    elif args.get('mjdstart') is not None:
        new_mjdrange = [args.get('mjdstart'), None]
    elif args.get('mjdend') is not None:
        new_mjdrange = [None, args.get('mjdend')]

    if new_mjdrange is not None:
        if config.pipe['fmjdselect.mjdrange'] is None:
            update_key(config.pipe, 'mjdrange', [new_mjdrange])
        else:
            update_key(config.pipe, 'mjdrange', config.pipe['fmjdselect.mjdrange'] + [new_mjdrange])

    new_fieldrange = None
    if args.get('fieldstart') is not None and args.get('fieldend') is not None:
        new_fieldrange = [args.get('fieldstart'), args.get('fieldend')]
    elif args.get('fieldstart') is not None:
        new_fieldrange = [args.get('fieldstart'), None]
    elif args.get('fieldend') is not None:
        new_fieldrange = [None, args.get('fieldend')]

    new_fieldrange = None
    if new_fieldrange is not None:
        if config.pipe['fmjdselect.fieldrange'] is None:
            update_key(config.pipe, 'fieldrange', [new_fieldrange])
        else:
            update_key(config.pipe, 'fieldrange', config.pipe['fmjdselect.fieldrange'] + [new_fieldrange])

    if new_fieldrange is not None:
        if config.pipe['fmjdselect.fieldrange'] is None:
            update_key(config.pipe, 'fieldrange', [new_fieldrange])
        else:
            update_key(config.pipe, 'fieldrange', config.pipe['fmjdselect.fieldrange'] + [new_fieldrange])

    for arg, value in args.items():
        if value is not None and arg not in exclude_args:
            if not update_key(config.pipe, arg, value):
                print(f"[DEBUG] Key '{arg}' not found anywhere in config.")

    if not config.pipe['fmjdselect.obs']:
        update_key(config.pipe, 'obs', ['apo','lco'])

    if set_gen:
        generations.set_config(config.pipe['fmjdselect.obs'])
