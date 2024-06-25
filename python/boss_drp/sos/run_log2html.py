#!/usr/bin/env python3
import boss_drp.utils.putils as putils
import re
import os.path as ptt

def run_soslog2html(lf, mjd, obs):
    cmd = f"sos_log2html, '{lf}', '{lf.replace('.fits','.html')}', /fps, /sdssv_sn2, obs='{obs}'"
    i = 0
    while i < 5:
        if i > 0:
            sleep(2)
            print("Trying again to get idl license")
        print("executing: " + cmd+"\n")
        rv = putils.runCommand(f'idl -e "{cmd}"', echo=False)
        if rv[0] != 0:
            print("Failed to update SOS HTML log")
        if 'Failed to acquire license.' not in rv[1]:
            break

    return
    squote = "\'"
    addstring = '</HEAD><BODY ONLOAD=\"timerID=setTimeout('+squote+'location.reload(true)'+squote+',60000)\">'
    currentfile = ptt.join(ptt.dirname(lf),'logfile-current.html')
    with open(lf.replace('.fits','.html'),'r') as source:
        lines = source.readlines()
    with open(currentfile,'w') as source:
        for line in lines:
            source.write(re.sub('</HEAD>',addstring, re.sub('BOSS Spectro','BOSS Spectro (Current)', line)))
    yesterday = str(int(mjd)-1)
    tomorrow = str(int(mjd)+1)
    with open(ptt.join(ptt.dirname(ptt.dirname(lf)),'combined',f'logfile-{mjd}.html'),'w') as source:
            source.write(re.sub('Yesterday: <A HREF=../'+yesterday+'/','Yesterday: <A HREF=',
                                re.sub('Tomorrow: <A HREF=../'+tomorrow+'/','Tomorrow: <A HREF=',line)))

    with open(ptt.join(ptt.dirname(ptt.dirname(currentfile)),'combined','logfile-current.html'),'w') as source:
            source.write(re.sub('Yesterday: <A HREF=../'+yesterday+'/','Yesterday: <A HREF=',
                                re.sub('Tomorrow: <A HREF=../'+tomorrow+'/','Tomorrow: <A HREF=',
                                        re.sub('</HEAD>',addstring,
                                                re.sub('BOSS Spectro','BOSS Spectro (Current)', line)))))
        
