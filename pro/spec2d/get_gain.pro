function get_gain, mjd, camname, obs
    yanny_read,filepath('opGain.par',root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles'),gaindata,/anonymous
    gains=*gaindata
    gains=gains[where(gains.MJD lt mjd and gains.camname eq STRLOWCASE(camname) $
                    and gains.obs eq STRUPCASE(obs))]
    mjd_gain=max(gains.mjd,i)
    splog, 'Gain For: ', gains[i].OBS, ' ', gains[i].camname, ' ', gains[i].mjd
    splog, gains[i].gain
    splog, 'Note: ', gains[i].Note
    return, gains[i].gain
end
