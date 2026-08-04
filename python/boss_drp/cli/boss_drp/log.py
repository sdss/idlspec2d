#!/usr/bin/env python3

# from boss_drp.utils.daily_log import (daily_log_email, daily_log_to_file,
#                                       daily_log_index)
from boss_drp.utils.argparse_help import AttrDict, _add_obs

from boss_drp.utils import jdate
from boss_drp.field import Field
from boss_drp import daily_dir
import glob
import os
import json
import os.path as ptt
from collections import OrderedDict

jdate = jdate.astype(int)
import click


@click.command(name='log')
@click.option("--obs", multiple=True, default=("apo", "lco"),
              type=click.Choice(["apo", "lco"], case_sensitive=False),
              help="Observatory for status update")
@click.option("--apo", is_flag=True, expose_value=False, help="Run apo",
              callback=lambda ctx, param, value: _add_obs(ctx, param, value, "apo"))
@click.option("--lco", is_flag=True, expose_value=False, help="Run lco",
              callback=lambda ctx, param, value: _add_obs(ctx, param, value, "lco"))
@click.option("--mjd", multiple=True, type=int, help="Update these MJDs")
@click.option("--mjdstart", type=int, help="Starting MJD")
@click.option("--mjdend", type=int, help="Ending MJD")
@click.option("--epoch", is_flag=True, help="Run for epoch Coadds")
@click.option("--custom", type=str, default=None, help="Name of custom Coadd")
@click.option("--topdir", default=os.getenv("BOSS_SPECTRO_REDUX"), type=str,
              help="Optional override value for the environment variable $BOSS_SPECTRO_REDUX")
@click.option("--run2d", default=os.getenv("RUN2D"), type=str,
              help="Optional override value for the enviro variable $RUN2D")
@click.option("--run1d", default=os.getenv("RUN1D"), type=str,
              help="Optional override value for the enviro variable $RUN1D")
@click.option("--email", is_flag=True, help="Send each mjd status as email")
@click.option("--fast", is_flag=True, help="Skip updating index until end")
@click.option("--refresh", is_flag=True, help="Refresh all the existing Status logs for obs")
@click.option("--refresh_error", is_flag=True, help="Refresh existing Status logs for obs with errors")
@click.option("--refresh_critical", is_flag=True, help="Refresh existing Status logs for obs with critical errors")
@click.option("--force", is_flag=True, help="Refresh Summaries pages")
@click.pass_context
def log(ctx, **kwrds):
    """BOSS Pipeline Status Log"""
    from boss_drp.utils.daily_log import (daily_log_email, daily_log_to_file,
                                      daily_log_index, valid_mjd)

    args = AttrDict(ctx.params)

    if args.run2d is not None:
        if args.run1d is None:
            args.run1d = args.run2d
    
    if args.epoch:
        dir_ = 'epoch'
    elif args.custom is not None:
        dir_ = args.custom
    else:
        dir_ = 'daily'
    
    if args.refresh:
        args.fast = True
        for obs in args.obs:
            mjds = glob.glob(ptt.join(daily_dir, 'logs', 'Status', dir_,
                                      args.run2d, f'?????-{obs.upper()}.html'))
            
            mjds = [int(ptt.basename(x).split('-')[0]) for x in mjds]
            mjds.sort()
            obs = obs.lower()
            for mjd in mjds:
                if not valid_mjd(mjd, args.mjd, args.mjdstart, args.mjdend):
                    continue
                print(args.run2d, mjd, obs)
                daily_log_to_file(obs, mjd, topdir=args.topdir, run2d=args.run2d,
                                  run1d=args.run1d, redux=None, html_log=None,
                                  summary=(not args.fast), epoch=args.epoch,
                                  custom = args.custom)
    elif args.refresh_error or args.refresh_critical:
        args.fast = True
        if args.refresh_error:
            _json_file = ptt.join(daily_dir, 'logs', 'Status', dir_,
                                      args.run2d,'error.json')
        else:
            _json_file = ptt.join(daily_dir, 'logs', 'Status', dir_,
                                      args.run2d,'critical.json')
        if ptt.exists(_json_file):
            try:
                with open(_json_file, 'r') as json_file:
                    errors = json.load(json_file, object_pairs_hook=OrderedDict)
                for obs in args.obs:
                    mjds = []
                    for err in errors:
                        if err['OBS'].lower() == obs.lower():
                            mjds.append(err['MJD'])
                    mjds = sorted(list(set(mjds)))
                    for mjd in mjds:
                        if not valid_mjd(mjd, args.mjd, args.mjdstart, args.mjdend):
                            continue
                        print(args.run2d, mjd, obs)
                        daily_log_to_file(obs, mjd, topdir=args.topdir, run2d=args.run2d,
                                          run1d=args.run1d, redux=None, html_log=None,
                                          summary=(not args.fast), epoch=args.epoch,
                                          custom = args.custom)
            except:
                print(f'Error opening {_json_file}... exiting')
                raise SystemExit(1)
        else:
            print(f'{_json_file} does not exist... exiting')
            raise SystemExit(1)
    else:
        if not args.mjd:
            if args.custom is None:
                if args.mjdend is None and args.mjdstart is None:
                    args.mjd = [jdate-1, jdate]
                elif args.mjdend is None and args.mjdstart is not None:
                    args.mjd = range(args.mjdstart, jdate+1)
                elif args.mjdstart is None and args.mjdend is not None:
                    args.mjd = [args.mjdend]
                else:
                    args.mjd = range(args.mjdstart, args.mjdend+1)
            else:
                fdir = Field(args.topdir, args.run2d, '{custom}_{obs}',
                            custom_name = args.custom, custom = True)
                fd = ptt.join(fdir.dir(), 'redux_{custom}_{obs}-?????')
                redux =      glob.glob(fd.format(custom=args.custom, obs='apo'))
                redux.extend(glob.glob(fd.format(custom=args.custom, obs='lco')))
                args.mjd = [int(x.split('-')[-1]) for x in redux]
                
                fdm = ptt.join(fdir.dir(), 'redux_{custom}_{obs}-?????_?????')
                redux1 =      glob.glob(fdm.format(custom=args.custom, obs='apo'))
                redux1.extend(glob.glob(fdm.format(custom=args.custom, obs='lco')))
                args.mjd.extend(list(set([int(x.split('-')[-1].split('_')[0]) for x in redux1])))
        elif args.custom is not None:
            fdir = Field(args.topdir, args.run2d, '{custom}_{obs}',
                        custom_name = args.custom, custom = True)
            fd = ptt.join(fdir.dir(), 'redux_{custom}_{obs}-?????')
            fdm = ptt.join(fdir.dir(), 'redux_{custom}_{obs}-?????_?????')
        args.mjd = list(dict.fromkeys(args.mjd))
        for mjd in args.mjd:
            if args.custom is None:
                it_obs = args.obs
            else:
                it_obs = [ptt.basename(x).split('-')[-2].split('_')[-1] for x in glob.glob(fd.format(custom=args.custom, obs='???').replace('?????',str(mjd)))]
                it_obs = [x for x in it_obs if x in args.obs]
            
                it_obs_m = [ptt.basename(x).split('-')[-2].split('_')[-1] for x in glob.glob(fdm.format(custom=args.custom, obs='???').replace('?????_',f'{mjd}_'))]
                it_obs.extend(list(set([x for x in it_obs_m if x in args.obs])))
               
            for obs in it_obs:
                obs = obs.lower()
                print(args.run2d, mjd, obs)
                if args.email:
                    if args.epoch:
                        subject = f'Epoch Status: {mjd} {obs}'
                    else:
                        subject = f'Status: {mjd} {obs}'
                    daily_log_email(subject, None, obs, mjd, content=None,
                                email_file = ptt.join(daily_dir, 'etc','emails'),
                                topdir=args.topdir, run2d=args.run2d, run1d=args.run1d,
                                redux = None, epoch=args.epoch,
                                custom=args.custom)
                else:
                    daily_log_to_file(obs, mjd, topdir=args.topdir, run2d=args.run2d,
                                run1d=args.run1d, redux=None, html_log=None,
                                summary=(not args.fast), epoch=args.epoch,
                                custom = args.custom)
    if args.force:
        args.mjd = None
    if args.fast or args.force:
        daily_log_index(ptt.join(daily_dir, 'logs', 'Status', dir_, args.run2d), args.run2d,
                        epoch=args.epoch, custom=args.custom, fast_mjds = args.mjd)
        fmjds =  args.mjd
    else:
        fmjds = None
    daily_log_index(ptt.join(daily_dir, 'logs', 'Status', dir_, args.run2d), args.run2d,
                    epoch=args.epoch, custom=args.custom, flag_noSci=True, fast_mjds = fmjds)
