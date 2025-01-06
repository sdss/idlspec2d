
function get_trace_dir, mjd, topdir=topdir, run2d=run2d
    if not keyword_set(topdir) then $
        topdir = getenv('BOSS_SPECTRO_REDUX')
    if not keyword_set(run2d) then run2d = getenv('RUN2D')
    dir_ = djs_filepath(strtrim(mjd,2), root_dir=topdir, $
                        subdirectory=[run2d,'trace'])
    return, dir_
end
