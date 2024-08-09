
Running a BOSS DRP Daily Run
============================


Daily Run
---------
In addition to running large batches of MJDs in one go, the pipeline is also designed to be run on a daily basis.

Crontab at Utah
^^^^^^^^^^^^^^^
The `crontab <https://man7.org/linux/man-pages/man5/crontab.5.html>`_ below has been designed so that two of these tasks can been written to allow for multiple version os the pipeline to be run in parallel (eg Master and |idlspec2d_version|). This is to ensure a stable reduction (only updated with new tags) and a development branch that will contain the leading edge developments, but will not contain a uniform reduction. In this mode, every time :ref:`uurundaily<uurundaily>` runs, it checks for a new MJD in the $BOSS_SPECTRO_DATA_N/S directories and if a new mjd exists, it updates the $daily_dir/etc/nextmjd.par file and starts run that mjd through the full pipeline. The numbers and asterisk, in the crontab below, before the command are when the script will run:

============  =====================
field         allowed values
============  =====================
minute        0-59
hour          0-23
day of month  1-31
month         1-12
day of week   0-7 (0 or 7 is Sunday
============  =====================

 ::

    15,45, * * * * cronrun.bash bhm/|idlspec2d_version| "uurundaily --module bhm/|idlspec2d_version| --lco --fast --merge3d --no_dither --monitor"
    20,20, * * * * cronrun.bash bhm/|idlspec2d_version| "uurundaily --module bhm/|idlspec2d_version| --apo --fast --merge3d --no_dither --monitor"
    0 2 * * * cronrun.bash bhm/|idlspec2d_version| "slurm_pySummary --module bhm/|idlspec2d_version| --full"
    30 22 * * * plot_QA.bash -m bhm/|idlspec2d_version| -n
    30 20 * * * plot_QA.bash -m bhm/|idlspec2d_version| -n -l

Manual Using uurundaily at Utah
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To manually use (:ref:`uurundaily<uurundaily>`) (which runs the full pipeline end-to-end), one of versions of the command below can be used

Run for any new MJDS (if a module/paths are preloaded)
""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

    uurundaily --lco --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor
    uurundaily --apo --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor

Run for a set MJD (if a module/paths are preloaded)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. note::
    This method does not update $daily_dir/etc/nextmjd.par
    
::

    uurundaily --lco --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor --mjd 60000
    uurundaily --apo --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor --mjd 60000

