#!/usr/bin/env python3
import re
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import numpy as np
import glob
import os.path as ptt
import os
from astropy.io import fits

# Function to find the closest timestamp before a given time
def find_closest_before(time, sorted_times):
    closest = None
    for t in sorted_times:
        if t < time:
            closest = t
        else:
            break
    return closest

# Function to parse timestamps from log file lines
def parse_log_file(file_path):
    start_times = []
    end_times = []
    end_exp = []
    write_exp = []
    start = False
    fail_idl = False
    with open(file_path, 'r') as file:
        for line in file:
            # Extract start times (assuming "start" keyword in log line)
            if "Exposure Flavor: science" in line:
                match = re.search(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}', line)
                if match:
                    start_times.append(match.group())
                    start = True
            elif start is True and '% Failed to acquire license.' in line:
                fail_idl = True
            # Extract end times (assuming "end" keyword in log line)
            elif "read_sos.py /data/boss/sos" in line:
                if not fail_idl:
                    mjd = line.split()[5]
                    ccd = line.split()[7].split('-')[2]
                    expid = line.split()[7].split('-')[3].split('.')[0]
                    ff = f'/data/spectro/{mjd}/sdR-{ccd}-{expid}.fit.gz'
                    write_exp.append(ptt.getctime(ff))
                    try:
                        hdr = fits.getheader(ff)
                        end_exp.append(hdr['INTEND'])
                    except:
                        end_exp.append(None)
                    match = re.search(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}', line)
                    if match:
                        end_times.append(match.group())
                        start = False
                else:
                    fail_idl = False

    return start_times, end_times, end_exp, write_exp

# Function to calculate elapsed times
def calculate_elapsed_times_from_logs(start_times, end_times, exp = False, create=False):
    if exp:
        start_times = [datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f') for t in start_times]
    elif create:
        start_times = [datetime.fromtimestamp(t) for t in start_times]
    else:
        start_times = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in start_times]
    end_times = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in end_times]

    elapsed_times = []
    for end in end_times:
        closest_start = find_closest_before(end, start_times)
        if closest_start:
            elapsed_times.append((end - closest_start).total_seconds())
        else:
            elapsed_times.append(None)

    return elapsed_times

def parse_runtime(file_paths, all=False, stamp=False):
    now = datetime.now()
    if now < now.replace(hour=18, minute=0, second=0, microsecond=0):
        current_date = now-timedelta(days=1)
    else:
        current_date = now
    current_date = current_date.strftime('%Y-%m-%d')

    for file_path in file_paths:
        dates = [current_date]

        # Parse log file and calculate elapsed times
        start_times, end_times, end_exps, write_exp = parse_log_file(file_path)
        elapsed_times = calculate_elapsed_times_from_logs(start_times, end_times)
        Full_elapsed_times = calculate_elapsed_times_from_logs(end_exps, end_times, exp=True)
        write_elapsed_times = calculate_elapsed_times_from_logs(write_exp, end_times, create=True)
        if all:
            for f in glob.glob(file_path + '.????-??-??'):
                dates.append(f.split('.')[-1])
                st, et, ee, we = parse_log_file(f)
                elapsed_times.extend(calculate_elapsed_times_from_logs(st, et))
                Full_elapsed_times.extend(calculate_elapsed_times_from_logs(ee, et, exp=True))
                write_elapsed_times.extend(calculate_elapsed_times_from_logs(we, et, create=True))

        if all:
            print('Using:', file_path+'.????-??-??')
        else:
            print('Using:', file_path)
        print('Nexp:', len(elapsed_times))
        print('Mean:', np.nanmean(elapsed_times))
        print('Standard Deviation:', np.nanstd(elapsed_times))
        print('Median:', np.nanmedian(elapsed_times))
        print('-------------------------------------------')
        print('Nexp full:', len(Full_elapsed_times))
        print('Mean full:', np.nanmean(Full_elapsed_times))
        print('Std full:', np.nanstd(Full_elapsed_times))
        print('Med full:', np.nanmedian(Full_elapsed_times))
        print('-------------------------------------------')
        print('Nexp write:', len(write_elapsed_times))
        print('Mean write:', np.nanmean(write_elapsed_times))
        print('Std write:', np.nanstd(write_elapsed_times))
        print('Med Write:', np.nanmedian(write_elapsed_times))

        # Plot the distribution of the elapsed times
        fig,ax = plt.subplots(1,figsize=(10, 6))
        bins = np.arange(min(np.nanmin(elapsed_times),np.nanmin(Full_elapsed_times),np.nanmin(write_elapsed_times))-.5,
                         max(np.nanmax(elapsed_times),np.nanmax(Full_elapsed_times),np.nanmax(write_elapsed_times))+.5,1)
        plt.hist([t for t in elapsed_times if t is not None], bins=bins, edgecolor='black',label='ReduxEnd-ReduxStart',alpha=.5)
        plt.hist([t for t in write_elapsed_times if t is not None], bins=bins, edgecolor='black',label='ReduxEnd-WriteEnd',alpha=.5)
        plt.hist([t for t in Full_elapsed_times if t is not None], bins=bins, edgecolor='black',label='ReduxEnd-ExpEnd',alpha=.5)
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(which='minor', length=4)

        plt.legend()
        dates = np.sort(np.asarray(dates)).tolist()
        plt.title('Distribution of Elapsed Times\n'+', '.join(dates))
        plt.xlabel('Elapsed Time (seconds)')
        plt.ylabel('Frequency')
        plt.grid(True)

        # Save the plot
        output_dir = '/data/boss/sos/logs/'
        cdate = '_'+current_date if stamp else ''

        try:
            plt.savefig(ptt.join(output_dir, f"{ptt.basename(file_path)}{cdate}.png"))
            print('Saving to: '+ptt.join(output_dir, f"{ptt.basename(file_path)}{cdate}.png"))
        except Exception as e:
            print(e)
            plt.savefig(f"{ptt.basename(file_path)}{current_date}.png")
            print('Saving to: '+f"{ptt.basename(file_path)}{current_date}.png")
        print("")
