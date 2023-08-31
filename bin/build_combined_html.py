#!/usr/bin/env python3
from glob import glob
import os.path as ptt
from os import rename
import time
def get_mjd(html):
    filename = ptt.splitext(ptt.basename(html))[0]
    int_part = filename.split('-')[1]
    return(int(int_part))

def build_combine_html():
    directory = '/data/boss/sos/combined/'
    if ptt.exists(ptt.join(directory,'index.html')):
        file_time = ptt.getmtime(ptt.join(directory,'index.html'))
        if time.time()-file_time < 3600*15:
            print('No update required')
            return
    with open(ptt.join(directory,'index.html.tmp'), 'w') as f:
        f.write('<html>\n')
        f.write(' <head>\n')
        f.write('   <title>BOSS Spectro Analysis</title>\n')
        f.write('   <style>BODY {font: normal 12pt Helvetica,Arial}</style>\n')
        f.write('   <link rel="shortcut icon" href="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-192x192.png" type="image/x-icon">\n')
        f.write(' </head>\n')
        f.write(' <body>\n')
        f.write('   <img src="sdss-new-logo.png" style="height: 78; width: auto;">\n')
        f.write('   <h2>BOSS Spectro Analysis</h2>\n')
        f.write('   <hr>\n')
        f.write('   <ui>\n')
        f.write('    <li><a href="logfile-current.html">logfile-current.html</a></li>\n')
        for html in sorted(glob(ptt.join(directory,'logfile-?????.html')),key=get_mjd,reverse=True):
            f.write('    <li><a href="'+ptt.basename(html)+'">'+ptt.basename(html)+'</a></li>\n')

        f.write('   <ui>\n')
        f.write(' </body>\n')
        f.write('</html>\n')
    rename(ptt.join(directory,'index.html.tmp'), ptt.join(directory,'index.html'))
if __name__ == '__main__':
    """
    Build sos/combined/index.html
    """
    build_combine_html()
