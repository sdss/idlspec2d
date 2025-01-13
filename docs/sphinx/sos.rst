.. title:: SOS: The BOSS on Mountain Pipeline

Son-of-spectro (SOS): The BOSS on Mountain Pipeline
===================================================
At APO (LCO) the log files are stored in |SOS_user|-hub.apo.nmsu.edu/data/boss/sos/<MJD> (|SOS_user|-hub.lco.cl/data/boss/sos/<MJD>) with the logs file stored in /home/|SOS_user|/boss/sos/logs.

The sos systemctl process is currently started automatically at boot.

The python help for this command can be found with SOS -h.


SOS as a systemctl process
--------------------------
SOS is designed to run as a pair of systemctl processes (1 for red and 1 for blue) at the observatories. They are controlled indepednetly, so that if one crashes it can be restared while leaving the other running. In |SOS_user| they are run by the |SOS_user| users on the |SOS_HOST| machines at each observatory. To start/restart/stop each of these processes, they can be controlled via ::

    systemctl --user start|restart|stop SOS
    systemctl --user start|restart|stop SOS_red

or as a joint command ::

    systemctl --user start|restart|stop SOS SOS_red

Backing out to an Older Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If errors in a recent update prevent SOS from functioning properly, older versions of SOS can be loaded instead. As the SOS version is managed by the idlspec2d modules on the mountains, this can be done by issuing the following commands as the |SOS_user| user on |SOS_HOST| ::

    systemctl --user stop SOS SOS_red
    cd /home/|SOS_user|/software/modulefiles/idlspec2d
    unlink default
    ln -s v6_0_9  default
    systemctl --user start SOS SOS_red

where v6_0_9 can be replaced by the latest functional version installed.

Troubleshooting SOS
^^^^^^^^^^^^^^^^^^^
The Son-of-Spectro (SOS) reduction has proven to be quite robust. However, if it appears to not be working one could check the following:

* Is the system process running? Log into "|SOS_user|@|SOS_HOST|", and issue :code:`systemctl --user status SOS SOS_red` command to check if both processes are still running.  If not then :code:`systemctl --user restart SOS SOS_red` (you can just to SOS or SOS_red if only 1 crashed). Then use the catchup commands to rerun any missed exposures and note in the night logs. If both processes are running, then just try the catchup command without restarting the processes (please note these exposures as well). If the reduction still fails then try emailing or slacking |Contact|

* Are there old "lock" files sitting around? SOS generates "lock" files to prevent multiple processes from changing these files simultaneously. If any have been sitting around for several minutes or more, then something is wrong. You should delete the file ending with ".lock". However, this could result in a corrupted "logfile*.fits".

The three places that could have a lockfile are ::

    /data/boss/sos/{mjd}/*.lock
    /data/boss/sos/combined/*.lock
    /data/boss/sos/{mjd}/trace/{mjd}/*.lock
    /home/sdss5/software/sdsscore/main/*/sdHdrfix/*.lock

It is normal for lock files to be there temporarily (for a few minutes); if code appears hung or the lockfiles are from a few days ago, then something is wrong. The lockfiles are links to the file they are locking, e.g. ::


    -rw-rw-r-- 1 sdss5 sdss5    542 Mar  9 01:17 sdHdrFix-57090.par
    lrwxrwxrwx 1 sdss5 sdss5     18 Mar  9 18:17 sdHdrFix-57090.par.lock -> sdHdrFix-57090.par

.. warning ::
    If you need to remove a lock, it is very important to remove the lock and not the file to which it is pointing. In the example above, remove the file with extension .par.lock, not .par file. i.e. issue :code:`rm $SDSSCORE_DIR/$OBSERVATORY/sdHdrfix/sdHdrFix-57090.par.lock`. DON'T issue :code:`rm $SDSSCORE_DIR/$OBSERVATORY/sdHdrfix/sdHdrFix-57090.par`.


Re-Reducing Data with SOS
-------------------------
In cases where data nees to be re-reduced, this can be done manually (using the |SOS_user| user on |SOS_HOST|) after loading the idlspec2d module. To reduce the blue and red channels independently ::

    SOS -b -c -m MJDXX -e EXPIDXX
    SOS -r -c -m MJDXX -e EXPIDXX

or jointly using ::

    SOS -r -c -m MJDXX -e EXPIDXX

the MJD can be negelected if the data is for the current MJD and the EXPID can be a single exposure ID or a range of exposure IDs.

.. note ::
    If the S/N in one of the BOSS cameras does not appear in Kronos, you can run this command to re-reduce an exposure using this command. You can check the log files at /home/|SOS_user|/boss/sos/logs/sos_log-b1-error or (sos_log-r1-error). tail -n 150 sos_log-b1-error will show the last 150 lines of the error long and should show enough of the log to determine if the current exposure crashed due to  "Failed to acquire license."


Frequently Used Commands
------------------------
.. note::
    All commands are executed from |SOS_user| on |SOS_HOST|. All commands assume that idlspec2d is set up  :code:`module load idlspec2d`.

Start SOS Processes After Observing Has Started For the Night:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use this command if the SOS process was not started until after observing has started for the night or an issue is identified later in the night. This procedure will process all images already present for the current MJD (or when utilizing the --mjd flag for the specified MJD). This can be combined with a --exp to reprocess a single exposure or range of exposures. Unlike the older version of sos, this command is designed to be run in parallel to the systemd processes allowing for reprocessing of the earlier data while new data is still being taken and processed. ::

    module load idlspec2d
    SOS -b -c
    SOS -r -c

these logs have _catchup added to the file names

Other flag
""""""""""

In addition to the catchup mode,  SOS has other options that can either be added to the -c option or run in place of the -c option

-e XXXXXXX  (--exp)     rerun a single or range of exposures (eg 100-110) exposure (works with -c, -m, or -t)

-m XXXXX    (--mjd)     rerun for an MJD (these logs have _{mjd} added to the filename) (can be combined with -c to reprocess a night in the sos folder or -t in the sosredo folder)

-t  (--redoMode)    rerun and save in (these logs have _redo added to the filename) (in place of -c) -This mode is primarily for testing and not designed for use by the observers.

-j  (--joint)   will run both red and blue in parallel (overrides -r or -b, and can be run with any of the other options)

--no_reject     will override the calibration rejection codes (to force a faint calibration through) - This is a last resort if the lamp dies during the night and replacement is not an option, since it will likely cause difficulties in the main pipeline that require trial & error and hand holding to reduce. SOS will still provide the same warnings but will produce output for calibration (note the override in the night logs). This should NOT be forced if flat 80% < 200 (as reported by SOS) or 2 lamps are out. If the outage of a single lamp happens at the start or early in the night (or multiple lamps burn out), then immediate replacement (either then or at the end of that field) is highly preferred.


SOS Outputs
-----------
Tabulated Values
^^^^^^^^^^^^^^^^

Son-of-Spectro reduces four flavors of observations: bias/dark, flat, arc, & science. Select information is tabulated for each of these types of observations. These values are tabulated in yellow if they are going out of spec, and in red if they are very much out of spec. The values reported are:

`BIAS/DARK PERCENTILE98`:     The value in electrons of the 98-th percentile on the (overscan-corrected) image. For example, if this is 8, then 98% of the pixels are below 8 electrons. Bad regions on the CCD and saturated pixels are excluded from this evaluation. The CCDs can also accumulate charge during the day that may need several bias exposures to completely flush. If this does not work, then the CCD is probably warm.

`FLAT NGOODFIBER`: The number of illuminated fibers. This should be 320 if all the fibers are plugged and unbroken. Note that fibers that fall on particularly bad parts of the CCD can also be excluded from the good fiber list (the red CCDs have some of these cases). If this number is less than about 315, then you should be suspicious that some fibers have dropped out.

`FLAT XMIN (XMAX)`: The minimum (maximum) X position of the spectra on the CCD. If this is less than 0 (greater than 2047), then some of the spectra fall off the left (right) side of the CCD. This probably means that no one ran the Spectro Monthly Checkout. If it is off by more than 5 pix or so, you should probably run the Monthly Checkout (if it has been run contact |Contact| to reduce and update the specflat product)

`FLAT XSIGMA`: The profiles in the spatial (X) dimension have a gaussian fit with a width of this sigma. The median of this width is taken independently in each of 4 quadrants on a CCD, and the maximum of those 4 values reported. If the spectrographs are in focus, then this value should be about 1.0 pix. If it is larger, then the spectrographs may be out-of-focus, or the slit-heads may not be properly latched.

`ARC WAVEMIN (WAVEMAX)`: The minimum (maximum) wavelength (in Angstroms) of any spectra on this CCD. Because of the optical distortions, this is always for the central fiber. The edge fibers have less wavelength coverage. This probably means that no one ran the Spectro Monthly Checkout. If it is off by more than 5 pix or so, you should probably run the Monthly Checkout. (if it has been run contact |Contact| to reduce and update the specflat product)

`ARC BESTCORR`: The linear correlation coefficient between the arc spectra and a template arc spectrum. If they agree perfectly, then this is 1. The value of this correlation is typically 0.80. If it drops too much lower, then the arc spectra do not look as they should. This could happen if the arcs did not turn on, the flat field screens did not close, or some of the arc lines are missing -- if, for example, the Hg lamps all failed. If the correlation is less than 0.5 or so, then it's even possible that the code has found the incorrect wavelength solution.

`ARC NLAMPS`: The number of arc lines used to generate the wavelength solution. There are many more on the red CCD's because neon has many more lines there. If this drops to too few lines, then some of the lamps have not turned on properly, have not warmed up, or have warmed up too much.

`ARC WSIGMA`: The arc-line profiles in the wavelength (Y) dimension have a gaussian fit with a width of this sigma. The median of this width is taken independently in each of 4 quadrants on a CCD, and the maximum of those 4 values reported. If the spectrographs are in focus, then this value should be about 1.0 pix. If it is larger, then the spectrographs may be out-of-focus, or the slit-heads may not be properly latched. (One would expect that both XSIGMA and WSIGMA would go out-of-focus at the same time.)

`SCIENCE SKY/SEC`: The median sky counts in electrons per pixel. If this is too large, then there must be a light source near the telescope, or the night sky (the moon) is bright, or the CCD's are warming up and generating dark current.

`SCIENCE (S/N)2`: The signal-to-noise squared for objects at the fiducial magnitude limit. We choose (S/N)2 because it is an additive quantity with additional integration time. Exactly how this quantity is measured is described below.

`SCIENCE (S/N)2_v2`: The signal-to-noise squared for objects at a fiducial magnitude limit. We choose (S/N)2 because it is an additive quantity with additional integration time. Exactly how this quantity is measured is described below.

`SCIENCE Mag15 (S/N)2`: The signal-to-noise squared for objects at a fiducial magnitude of 15 limit. We choose (S/N)2 because it is an additive quantity with additional integration time. Exactly how this quantity is measured is described below.

`EXPTIME`: The exposure time (EXPTIME) header keyword from the first camera of this exposure to be reduced. The assumption is that all 2 cameras have the same EXPTIME.

`AIRTEMP`: The air temperature (AIRTEMP) header keyword from the first camera of this exposure to be reduced. The assumption is that all 2 cameras have the same AIRTEMP.

`UT`: The UT time computed from TAI in the header from the first camera of this exposure to be reduced.

`QUALITY`: This is the observer-input quality for this exposure. It can be set independently for each of the 2 cameras, but only one camera value (the first reduced frame) is reported in this table. The default value is "excellent" for everything except for dithered flats or spectro focus frames which are "test". The observers have the option of declaring exposures "excellent", "test", or "bad" using the sdr_hdrfix.py command.

The exact values of these yellow/red limits and further explanation can be found in the "idlspec2d" product in the file "examples/opLimits.par".

WARNING and ABORT messages
^^^^^^^^^^^^^^^^^^^^^^^^^^
There are a number of WARNING and ABORT messages that can appear if the pipeline runs into trouble when processing a frame. These messages appear at the bottom of each table. Each one-line message begins with the relevant file name, WARNING or ABORT, then a brief plain-text message.

Note that a single problem may cascade into a large number of warning messages. For example, an out-of-focus spectrograph will first produce the "Median spatial widths" message, probably followed by warnings about bad sky-residuals.



General Frame Messages (valid for flats/arcs/Science)
"""""""""""""""""""""""""""""""""""""""""""""""""""""
`bias exploding in b2 crazy quadrant`: This indicates a hardware electronics problem with the b2 bias levels. Keep observing, but note this prominently in the night log. The SOS b2 (S/N)2 values might be crazy -- if they look unusual, use the other channels to determine plate doneness, and note this in the log.

`Amp #... expected read noise = ..., measured = ... DN`: The number reported is the standard deviation in the bias region for either amplifier #2 (left side) or amplifier #3 (right side). This calculation is done clipping the half-percent of lowest and highest values. We trigger this warning if the value is ever 1.0 DN above the expected value. This can happen if there are a huge number of cosmic rays (if you've been integrating for hours), or if there is an electronics problem.

`Amp #... bias region difference at xxx-th-percentile =... DN:` This measures another statistic of the bias region for either amplifier #2 (left side) or amplifier #3 (right side). This measures the difference between the 16th-percentile and 84th-percentile, which should be equal to twice the read noise (e.g., 1-sigma). This test is done at the 68.2-percentile (1 sigma), 95.4-percentile (2 sigma), and 99.7-percentile (3 sigma). A warning is reported if this difference is either significantly too small or too large, as compared to what gaussian statistics dictate. This should catch the same sorts of electronics problems as listed for the above warning message.

`Amp #... way too many pixels (xxx%) below bias-5*sigma=... DN`: This test looks for anomalously low pixel values in the data region of the CCD. There should essentially never be any pixel values below 5-sigma less than the bias level, unless there is something wrong with the electronics.

`Fixing shifted rows (from electronics)`: This is very, very bad. The raw images have rows shifted to the right, sometimes by many pixels. Call José immediately.

`Fixing dropped-pixel rows (from electronics)`: This is very bad.  The raw images actually have some rows shifted by 1 or 2 pixels, usually more so near the bottom of the CCD (the first rows to be read). Call José immediately.

`Electronics shifted xxx rows by 2 pix`: This is indicative of an electronics problem with the BOSS electronics that shifts all or some of the rows in the raw images in all 4 amplifiers

`More than 10% of the image is rejected`: This is very bad. This can probably only happen if most of the CCD has saturated pixels, which probably means you're observing during twilight, the CCDs are warm, or the dome lights are on.

Science Frame Messages
""""""""""""""""""""""
* Unable to reduce science exposure:

  `Unable to reduce this science exposure (need flat)`: A valid reduced flat field frame is required for processing

  `Unable to reduce this science exposure (need arc)`: A valid reduced arc lamp frame is required for processing

  `Reject science ...`: A science exposure can be rejected if the header keywords indicate that the flat-field petals are closed, any flat-field or arc lamps are turned on, too many pixels are bad or saturated, or if the 25-th percentile of the image too large.

  `Reject science: Flat-field screens are closed!`

  `Reject science: Flat-field lamps turned on!`

  `Reject science as too bright: 25-th-percentile = ....`

  `Scattered light`: There was a high baseline count rate on the CCD, and appears even between fibers. This can most obviously occur if there are light sources in the CCD (such as the LED's we had for some time), or if the CCD's are warming up and generating dark current. There can also be a scattered light contribution just from a very bright sky, or if there are super-bright objects on some fibers that are scattering or bleeding light across the CCD. It's best to carefully inspect the raw images for problems.

  `Large flexure flat<->science`: There is a large shift (more than 1.00 pix) between the flat-field and the science exposure, presumably from flexure in the spectrographs. When this happens, another set of flat-fields (post-calibs) is recommended. (However, don't bother to take a set of calibrations on a different night. The Spectro-2D reductions never use calibrations from one night for science exposures on another.)

  `Whopping fiber ....`: The fibers listed have very bright objects that affect their neighbors on the CCD. If the objects are bright enough (12th mag?), then this can trigger other warnings such as "scattered light". I think it's safe to say that whopping fibers only show up when there has been a mistake made in the plate designs -- this is already recorded in PR 2471. You should check that this object does not saturate (> 30,000 ADU) the raw sdR image. Should it be saturating, reduce the exposure time to prevent saturation or move on to the next plate.


* The following warning messages are all based upon the quality of the sky-subtraction. Typically, we are able to model the sky spectrum (from the 16 sky fibers on each spectrograph) with a relative chi2 of around unity. At very bright sky lines, like 5577 Ang, the relative chi2 may be as large at 25. If the relative chi2 is large at other wavelengths, this means that there is excess light down the fibers that vary across the plate. This could be due to a light source near the telescope, or possibly a bright, non-uniform sky. Strong auroral activity is something that can do this, since the O I lines at 6300 and 6366 Ang are resolved on the sky.

  `Too few sky fibers to model sky-sub variance`: There are not enough good sky fibers on the CCD. This will only happen if for some reason there were far fewer than the mandated 16 sky fibers on a CCD, or most of those fibers just happened to be dead fibers. You should never see this message. If you do, the data is un-reducable. You could try one more exposure, but something is probably horribly wrong with this plate.

  `Median sky-residual chi2 = ... at ... Ang`: The median reduced chi2 for sky-subtraction is always around unity. If it is greater than 2, this message is triggered. If this message appears and you are not observing during twilight or with warm CCDs, then there must be a serious problem. Seriously out-of-focus spectrographs might trigger this, or lights on near the telescope.

  `Max sky-residual chi2 = ... at ... Ang (ignoring 5577)`: This is an informational message triggered if the reduced chi2 is greater than 100 anywhere other than at the 5577 Ang sky line. If there is auroral activity, then this could be triggered at a few specific lines like O I at 6300 and 6366 Ang. These O I lines has reduced chi2 values of about 100 on MJD 51999 during a solar storm (i.e. plate 336/51999). The Spectro-2D pipeline will automatically mask these wavelengths for any downstream analyses.

  `Bad sky residuals at ....`: This warning is triggered if there is a range of at least 25 Angstroms with a reduced chi in the sky-subtraction worse than 2. This could be due to a light source near the telescope, a bright, non-uniform sky, a warm CCD, scattered light, or out-of-focus spectrographs. For most of those cases, there should be other relevant warning messages preceding this one (like a warm CCD message). This message should be ignored at the edges of the CCDs wavelength coverage -- near 3800 or 6200 Ang for b1,b2 or near 5800 or 9200 Ang for r1,r2. If none of the above explanations apply, then there is some real problem.

  `Red Monster at ....`: This warning is triggered if we see bad sky residuals (above) that is specifically in the wavelength range of about 6400-6600 Ang. We think this happens when the handpaddle is still plugged in, in which case the observers should unplug it immediately. The threshold is set at reduced chi=2, and the worst that we have ever seen is at about the level chi=6. ii You should be able to see a bump in the extracted spectra, especially for the sky fibers.

Flat Frame Messages
"""""""""""""""""""
`Flat-field screens not closed`: The "FFS" keyword in the header indicates that at least one of the flat-field petals was not closed. When this happens, the flat or arc is not reduced.

`Flat-field lamps not turned on`: The "FF" keyword in the header indicates that at least one of the four flat-field lamps was not turned on. When this happens, the flat is not reduced.

`Reject flat (or arc) ... % bad pixels`: This condition is triggered when more than 2% of the non-masked pixels on the image are bad (saturated). When this happens, the flat (or arc) is not reduced. This probably happens if the CCDs are warm, the dome lights are on, or if for some reason the shutters were open too long.

`Reject flat (or arc) ... saturated rows`: This condition is triggered when there are more than 100 saturated rows on the image. When this happens, the flat (or arc) is not reduced. This probably happens if the CCDs are warm, the dome lights are on, or if for some reason the shutters were open too long.

`Reject flat as too faint`: This condition is triggered when the 80-th percentile of the image is less than 1000 electrons. When this happens, the flat is not reduced. Either the flat field screens were not closed, the lamps were not turned on, or the shutter didn't open.

`Possible Argon lines in superflat`: Emission lines are present in the quartz-halogen flat-field images, which are supposed to be featureless. When a number follows this message, that is a measure of the line strength -- the trigger is set to 0.01, but we usually see it as 0.1 to 0.5 when present. We have identified these rogue lines as argon. best guess is that these contaminating lines come from trace amounts of argon in the HgCd lamps, which must still have current running through them when they are supposed to be off.

`Median spatial widths = ...`: The spatial (X) widths of the fibers typically are gaussians with a sigma of 0.85 to 1.05 pixels. If the sigma is larger than 1.10 pixels in any of the 4 quadrants of a CCD (lower-left, lower-right, upper-left, upper-right), then this warning is triggered. It means that either the spectrographs are out of focus, or the slit-heads are not properly latched. Note that these widths are computed for both the flat and science exposures (but not for arcs -- we compute the widths in the dispersion dimension for arcs).

`Unable to reduce this flat exposure (need plug-map)`: There is something wrong with the plug-map  or confSummary file (You should not see this in the FPS era, but if you do contact José)

`All centroids rejecte`d`: Something went significantly wrong with the tracing. Try taking another flat and if this error appears again contact the |Contact| Immediately

Arc Frame Messages
""""""""""""""""""
`Flat-field screens not closed`: The "FFS" keyword in the header indicates that at least one of the flat-field petals was not closed. When this happens, the flat or arc is not reduced.

`Reject flat (or arc) ... % bad pixels`: This condition is triggered when more than 2% of the non-masked pixels on the image are bad (saturated). When this happens, the flat (or arc) is not reduced. This probably happens if the CCDs are warm, the dome lights are on, or if for some reason the shutters were open too long.

`Reject flat (or arc) ... saturated rows`: This condition is triggered when there are more than 100 saturated rows on the image. When this happens, the flat (or arc) is not reduced. This probably happens if the CCDs are warm, the dome lights are on, or if for some reason the shutters were open too long.

`X/4 .... Arc lamps are off`: The "NE", "HGCD", "HEAR" keywords in the header indicates that either the NE, HeAr, or HgCd lamps are not fully turned on. When this happens, the arc is not reduced.

`Arc lamps not turned on`: The "NE", "HGCD", "HEAR" keywords in the header indicates that either the NE, HeAr, or HgCd lamps are not fully turned on. When this happens, the arc is not reduced. (replaced by the message above)

`Neither Ne nor HeAR (HgCd) lamps turned on!`: The "NE", "HGCD", "HEAR" keywords in the header indicates that either the NE, HeAr, or HgCd lamps are not fully turned on. When this happens, the arc is not reduced.

`Best arc correlation = ...`: The cross-correlation of the arc spectrum with the template arc spectrum was bad. This problem also triggers BESTCORR as bad in the table (see above). This may be due to too little signal in the lamps, e.g. if the flat field petals did not close or the lamps did not turn on. If those are not the problems, then look at the raw image.

`Big wavelength gap`: Some arc lines were not found in the arc spectra, and there is a large gap in wavelength space without any lines. This will produce a poor wavelength solution. More arcs should be taken until this message does not appear.

`Median wavelength widths = ...`: The widths of the arc lines in the dispersion (Y) dimension typically are gaussians with a sigma of 0.90 to 1.10 pixels. If the sigma is larger than 1.10 pixels in any of the 4 quadrants of a CCD (lower-left, lower-right, upper-left, upper-right), then this warning is triggered. It means that either the spectrographs are out of focus, or the slit-heads are not properly latched.

`Arc exposure, waiting for flat before reducing`: Take a flat, and if that flat is valid then arc should be reduced, if not, then take another arc (or use SOS -j -c -e XXXX where XXXX is the exposure number to re-reduce the arc after getting the flat)

`Reject arc image too few lines`: The lamps were not properly on or warmed up. Take an additional Arc

`Cd I 3610 line missing (lamps not warm?)`: This particular arc line is the tie-down of the UV wavelength calibration and is either missing or has a poor fit. This is likely due to the HgCd lamps not being sufficiently warm before the exposure start. (depreciated)

`Wavelength mapping makes no sense!`: Something has going significantly wrong. Try taking another arc, if you get this again, visually check the arc to look for artifacts and the lamp status, and then check the focus. If all seems in order contact |Contact| for further investigation.

`Spline fit failed`: Something has gone wrong with the flat-arc pair. Try taking another pair of flat and arc.

Misc messages
"""""""""""""
`Airmass range = ...`: This warning is triggered if the airmass exceeds 2.5. At such large airmasses, the atmospheric refraction terms are getting large and the sky background is bright. The data is still perfectly usable, it's just not as good as taking data at lower airmass. Keep observing if you must, but the airmass will rapidly approach infinity!

`SOS disk is ...% full`: The specified disk is more than 95% full. Contact operations list to have Pipeline/Data teams clean out old data

`Sun above the horizon by ... deg for non-test exposure`: This warning message is to trap flat, arc, or science exposures taken during the day that have not been marked as either "test" or "bad" data. If these are test data, be sure to mark them as such. This is to prevent such data from being used later in the full reductions.

`Exposure number in header disagrees w/ filename`: This means that something is in a very confused state, and it is putting a different exposure number in the header (EXPOSURE keyword) from what is being used to generate file names. You should probably contact restart the BOSS(APO) or YAO(LCO) actors and contact José

`Error reading sdHdrFix file`: The sdHdrFix file (sdHdrFix-$MJD.par) in /home/sdss5/software/sdsscore/main/\*/sdHdrfix/ exists but is not a valid Yanny parameter file. This file is written by the procedure sdr_hdrfix.py if the observers have run that proc on any of the sdss5 machines, but can also be edited by hand. If the file is invalid, you should edit it to be valid, or delete it.

`Wrong number of elements for REDDEN_MED`: This is a warning that the reddening vector for this plate (reddeningMed in the plPlugMapM file) is not a 5-element vector, as expected). In this case, reddening values of zero are assumed. This means that the plate will not be observed as deep as it would have if non-zero reddening values were provided. (depreciated)


S/N Figures
^^^^^^^^^^^
A median signal-to-noise is computed for each object in the wavelength ranges [4000,5500] Angstroms (synthetic g-band) and [6910,8500] Angstroms (synthetic i-band). We plot these S/N values versus the PHOTO fiber magnitudes, which were measured in approximately a 3-arcsec diameter aperature (the same size as our fibers). If everything is working perfectly, then our S/N values should correlate very well with these PHOTO magnitudes.

We determine whether a plate is "done" based upon the signal-to-noise of the fainter objects on the plate. We do this by fitting a line to the (S/N)-vs.-magnitude plot in a specified wavelength range, then evaluating this fit at g=22.0 mag (blue CCDs), and i=21.0 mag (red CCDs). For the v2 S/N estimates the fits are evaluated at g=19 and r=19, and for the mag 15 S/N estimates the fits are evaluated at g=15 and r=15.

When the sky level is higher, we gain S/N more slowly at the fainter magnitudes where we are sky-limited rather than photon-limited.

Without moon, we have historically found ::
    
    log(S/N_g) = (zeropoint) - 0.31 * g
    log(S/N_i) = (zeropoint) - 0.31 * i
    
With partial moon, the slope steepens to -0.34 or worse.

However, in practice in SDSS-V we utilize ::

    log(S/N_g) = (zeropoint) - 0.32 * g
    log(S/N_i) = (zeropoint) - 0.36 * i
    
The fitting regions are denoted on the plot with vertical dotted lines. Arrows point to the evaluation of the fit on each of the 4 cameras, with the top panels corresponding to the blue CCD's (synthetic g-band) and the bottom panels corresponding to the red (synthetic r-band). In SDSS-V FPS, if less then 10 targets have magnitudes within the fit range, the fit range bright limit is extened to include all target.

The right-hand figures plot the residuals of each object from the fit. Symbol sizes on those right-hand plots indicate the magnitude of the deviation from the fit line. Symbol color is the same on the left as on the right, so green objects have more flux and red ones less. If the scale of the telescope is wrong, then you will see a radial drop-off in flux (red points on the edge of the plate). If you are observing too far over in air mass, then typically you correct to first order with a scale change, but a quadropole is left in these residuals. If one spectrograph has problems, then this will show up as red points in half of one of these figures.

Note that the (S/N)2 totals listed in the table and the figure might not exactly agree. This is because the fitting to (S/N)-vs.-magnitude is done on individual frames for the table, but on the summed S/N over all frames for the figure. The tabulated values are the ones we use to declare a plate done.


Throughput Summary Analysis Figures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the FPS era we added additional throughput analysis figures (linked via the SOS Summary Plots link at the top of the page). These figures show flux and S/N^2 as a functions of fiberid and magnitude. The Mag to Flux and Mag to S/N^2 figures include reference lines measured from a typically SDSS-IV/V plate to help understand throughput loss due to FPS fiber placements. This is also expressed in the x-focal vs y-focal vs log(ref_plate_flux/flux), which can be useful as a quick check of the overall fiber placement scale factors. The final panel is the extracted 1d spectra (though typically the scaling makes this figure hard to read).

Arc To Trace QA Figures
^^^^^^^^^^^^^^^^^^^^^^^
In the FPS era we have added the ability to use arc frames to tweak the traces as measured at the start of night operations. If this option is active, QA figures are produced of the trace shifts and linked via the Arc Shift Plots link at the top of the page. In the figures with 4 panels shows the original meaured shifts. The additional 3 panels shows the various transformations. Red vectors in the plots show discard points in the transformations. In the figure with 2 panels the top shows the net shift per fiber along the traces. And the bottom panel shows (if a trace flat is also taken for the field of the arc) the difference between the arc tweaked traces and the traces as measured by the flat for that field. In the lower panel, the dashed tails are portions of traces that are not used for the science frames.

Log Files
^^^^^^^^^
If the reduction of an exposure is catastrophically bad, it may not appear at all in the Son-of-Spectro table. However, there should still be a log file for this exposure on the data drive: /data/boss/sos/$MJD/splog-$CAMERA-$EXPOSURE.log

Reading this file should tell you what failed. The first and last lines of these files should contain "Started at" and "Finished at" followed by timestamps. If this does not provide you with any information you can check the latest process logs for the camera in /home/|SOS_user|/boss/sos/logs.


SOS Setup Requirement: Module
-----------------------------
At present, the SOS setup is managed via the idlspec2d module files at the observatories.
An example is shown below. Most of the requirements are the same as the main idlspec2d module,
however :code:`IDLSPEC2D_SOS`, :code:`BOSS_SPECTRO_DATA_N`, and :code:`BOSS_SPECTRO_DATA_S`
environmental variables should be set. :code:`IDLSPEC2D_SOS` tells the pipeline that it is running
in :code:`SOS` only mode. :code:`BOSS_SPECTRO_DATA_N` and :code:`BOSS_SPECTRO_DATA_S` are required always,
but at Utah they handled via other modules, so they should be set here manually.

```
#%Module5.0
# The first line of this file tells Modules that this is a module file.
# DO NOT ALTER IT!

proc ModulesHelp { } {
    global product version
    puts stderr "This module adds $product/$version to your environment."
}

set product idlspec2d
set version v6_2_0


module-whatis "Sets up $product/$version in your environment."


#
# DEPENDENCIES SECTION
#
# If your product requires other software to function, that should be declared
# here.  There are two types of dependencies: mandatory and optional.
# A mandatory dependency is a module load command followed by a prereq
# command.  An optional dependency is not followed by a prereq statement.
#
module unload specflat
module load specflat
module load sdsscore
module unload idlutils
module load idlutils/fps_boss
prereq idlutils/fps_boss
module load pyvista
#
# ENVIRONMENT SECTION
#
# The PRODUCT_ROOT and PRODUCT_DIR variables are used to set other
# environment variables, exported to the actual environment, by sdss4install
#
set PRODUCT_ROOT /home/sdss5/software
set PRODUCT_DIR $PRODUCT_ROOT/$product/$version
#
# This line creates an environment variable pointing to the install
# directory of your product.
#
setenv [string toupper $product]_DIR $PRODUCT_DIR
setenv [string toupper $product]_VER $version
setenv [string toupper $product]_SOS 1

# The lines below set various other environment variables.
setenv BOSS_SPECTRO_DATA_N /data/spectro/
setenv BOSS_SPECTRO_DATA_S /data/spectro/
setenv BOSS_DRP_DAILY_DIR /home/sdss5/boss/
setenv BOSS_QA_DIR /home/sdss5/boss/
setenv PYENV_VERSION idlspec2d-dev
setenv PYTHONUNBUFFERED 1

# Define SDHDRFIX_DIR
setenv SDHDRFIX_DIR /home/sdss5/software/sdsscore/main

append-path IDL_PATH  +$PRODUCT_DIR
prepend-path IDL_PATH +$PRODUCT_DIR/pro
prepend-path IDL_PATH +/usr/local/harris/idl88/lib
prepend-path IDL_PATH +/usr/local/harris/idl88/lib/obsolete
prepend-path IDL_PATH +/usr/local/harris/idl88/lib/graphics
prepend-path PATH $PRODUCT_DIR/bin
prepend-path PYTHONPATH $PRODUCT_DIR/python

```
