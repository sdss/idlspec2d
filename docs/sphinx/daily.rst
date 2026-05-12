
Running a BOSS DRP Daily Run
============================


Daily Run
---------
In addition to running large batches of MJDs in one go, the pipeline is also designed to be run on a daily basis.

Crontab at Utah
^^^^^^^^^^^^^^^
The `crontab <https://man7.org/linux/man-pages/man5/crontab.5.html>`_ below has been designed so that two of these tasks can been written to allow for multiple version os the pipeline to be run in parallel (eg Master and |idlspec2d_version|). This is to ensure a stable reduction (only updated with new tags) and a development branch that will contain the leading edge developments, but will not contain a uniform reduction. In this mode, every time :ref:`boss_drp daily<boss_drp_daily_py>` (or :ref:`uurundaily<uurundaily>`) runs, it checks for a new MJD in the $BOSS_SPECTRO_DATA_N/S directories and if a new mjd exists, it updates the $daily_dir/etc/nextmjd.par file and starts run that mjd through the full pipeline. The numbers and asterisk, in the crontab below, before the command are when the script will run:

============  =====================
field         allowed values
============  =====================
minute        0-59
hour          0-23
day of month  1-31
month         1-12
day of week   0-7 (0 or 7 is Sunday
============  =====================

Updated Command Line Interface
""""""""""""""""""""""""""""""
Using the new unified command line interface in combination with the configuration files, the Crontab syntax and setup is simplified:

.. admonition:: Updated CLI
     
    ::
    
        SHELL=/usr/bin/bash
        DAILY_DIR=/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux/control/daily/
        BOSS_DPR_PIPE_CONFIG_PATH=$BOSS_SPECTRO_REDUX/config/boss_drp.yml
        BOSS_DPR_PIPE_QUEUE_PATH=$BOSS_SPECTRO_REDUX/config/queue.yml
        IDLSPEC2D_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/idlspec2d/daily"

        TMOD="bhm/|idlspec2d_version|"

        0,30 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash $TMOD "boss_drp daily --apo"
        5,35 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash $TMOD "boss_drp daily --lco"
        0    1  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash TMOD "boss_drp batch Summary --defaults"
        3    *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash $TMOD "boss_drp log"
        0    22 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash TMOD -n
        0    20 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash TMOD -l -n        
        0    2  * * 6 $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash $TMOD "boss_flatlib end2end -cnl --mjdstart -10"
        0    4  * * 6 $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash $TMOD "boss_flatlib end2end -cn --mjdstart -10"
        0    5  * * * $IDLSPEC2D_DIR/bin/cronrun.bash $TMOD "sas_mos_too boss -t all -d v6_2_0 --mjd -1"


Legacy Command Line Interface
"""""""""""""""""""""""""""""
Using the legacy command line interface, the different flags must be specified for each version:

.. admonition:: Master/main Daily Run
     
    ::
    
        SHELL=/usr/bin/bash
        DAILY_DIR=/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux/control/daily/
        IDLSPEC2D_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/idlspec2d/daily"
        TMOD="bhm/master"
        0,30 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "uurundaily --apo --daily"
        5,35 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "uurundaily --lco --daily"
        0    1  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "slurm_Summary --defaults"
        3    *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "daily_log"
        0    22 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash bhm/$TMOD -n
        0    20 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash bhm/$TMOD -l -n
        
.. admonition:: Tagged Daily Run
     
    ::
    
        SHELL=/usr/bin/bash
        DAILY_DIR=/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux/control/daily/
        IDLSPEC2D_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/idlspec2d/|idlspec2d_version|"
        TMOD="bhm/|idlspec2d_version|"
        0,30 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "uurundaily --apo --tagged"
        5,35 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "uurundaily --lco --tagged"
        0    1  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "slurm_Summary --defaults"
        3    *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "daily_log"
        0    22 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash bhm/$TMOD -n
        0    20 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash bhm/$TMOD -l -n
        0    2  * * 6 $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "flatlib end2end -cnl --mjdstart -10"
        0    4  * * 6 $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/$TMOD "flatlib end2end -cn --mjdstart -10"
        0    5  * * * $IDLSPEC2D_DIR/bin/cronrun.bash $TMOD "sas_mos_too boss -t all -d v6_2_0 --mjd -1"

.. admonition:: Development Daily Run
     
    ::
    
        SHELL=/usr/bin/bash
        DAILY_DIR=/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux/control/daily/
        IDLSPEC2D_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/idlspec2d/|idlspec2d_version|_alpha"
        TMOD="bhm/|idlspec2d_version|_alpha"
        0,30 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash work/$TMOD "uurundaily --apo --dev"
        5,35 *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash work/$TMOD "uurundaily --lco --dev"
        0    1  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash work/$TMOD "slurm_Summary --defaults"
        3    *  * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash work/$TMOD "daily_log"
        0    20 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash work/$TMOD -n
        0    20 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronplot_QA.bash work/$TMOD -l -n
.. admonition:: Run QA Plots for two run2ds at once

    ::
    
        55 23 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/|idlspec2d_version| "plot_qa --run2d master |idlspec2d_version| --mjds_low 60403 None --mjds_high None 60402 --lco"
        25 20 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/|idlspec2d_version| "plot_qa --run2d master |idlspec2d_version| --mjds_low 60403 None --mjds_high None 60402 --lco --html"
        35 22 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/|idlspec2d_version| "plot_qa --run2d master |idlspec2d_version| --mjds_low 60447 None --mjds_high None 60446"
        30 23 * * * $IDLSPEC2D_BASE/$TRUN2D/bin/cronrun.bash bhm/|idlspec2d_version| "plot_qa --run2d master |idlspec2d_version| --mjds_low 60447 None --mjds_high None 60446 --html"



Manual Using uurundaily at Utah
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To manually use (:ref:`boss_drp daily<boss_drp_daily_py>` or :ref:`uurundaily<uurundaily>`) (which runs the full pipeline end-to-end), one of versions of the command below can be used. The ``--daily``, ``--tagged``, and ``--dev`` tags can also be used, but the examples below give the explicite options.

Run for any new MJDS (if a module/paths are preloaded)
""""""""""""""""""""""""""""""""""""""""""""""""""""""
.. code-block:: shell

    boss_drp daily --lco
    boss_drp daily --apo

If you are using the legacy Command Line interface:

.. code-block:: shell

    uurundaily --lco --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor
    uurundaily --apo --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor



Run for a set MJD (if a module/paths are preloaded)
"""""""""""""""""""""""""""""""""""""""""""""""""""

.. note::
    This method does not update $daily_dir/etc/nextmjd.par

.. code-block:: shell

    boss_drp daily --lco --mjd 60000
    boss_drp daily --apo --mjd 60000

If you are using the legacy Command Line interface:

.. code-block:: shell

    uurundaily --lco --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor --mjd 60000
    uurundaily --apo --module bhm/|idlspec2d_version|  --fast --merge3d --no_dither --monitor --mjd 60000


