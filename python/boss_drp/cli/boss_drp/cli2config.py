from boss_drp.field.generations import generations
from boss_drp.Config import config, update_key

def cli2config(args, config_par = None, exclude = None, set_gen = True):
    config_par = config_par or {}
    exclude = exclude or []
    config.load(args.queue_config, queue_config_file= args.queue_config_file,
                 config_name=args.config, config_file=args.config_file)
    for name, val in config_par.items():
        config.queue.set(name, val)

    exclude_args = ['config','config_file','queue_config','queue_config_file','show_config']
    exclude_args.extend(exclude)

    for arg, value in args.items():
        if value is not None and arg not in exclude_args:
            if not update_key(config.pipe, arg, value):
                print(f"[DEBUG] Key '{arg}' not found anywhere in config.")

    if not config.pipe['fmjdselect.obs']:
        update_key(config.pipe, 'obs', ['apo','lco'])

    if set_gen:
        generations.set_config(config.pipe['fmjdselect.obs'])
