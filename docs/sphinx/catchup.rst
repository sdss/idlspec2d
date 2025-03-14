Running a BOSS DRP Catchup Run
==============================

Daily Coadds
^^^^^^^^^^^^
.. note::
    Several of the commands contain options to bundle the slurm jobs (nbundle). This is suggested if you have a large number of tasks, as each task gets added as a row to the slurm database. If bundeling is used, each bundle gets added instead. e.g. 10,000 tasks would take a while to load 10,000 rows appended to the task table, but if you set nbundle=10 then it will only create 1000 rows in the database, and each bundled set of task would then be treated as a single task in slurm. The main disadvantage of this option is that it coarse-grains the percent complete reported while monitoring the run, so it is best to keep nubundle small to prevent gross coarse graining of the percent complete.

build the spplan files
""""""""""""""""""""""
The BOSS pipeline operation centers on a set of plan files built with the command. ::

    spplan --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --sdssv --no_dither --quick --apo --log apo_plan.log
    spplan --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --sdssv --no_dither --quick --lco --log lco_plan.log

.. admonition:: Edit manual plans
        There are some situations where the automated proceedure to build the spplan2d files fails to build an optimal or function plan.
        In which case running the following command on within the ``$BOSS_SPECTRO_REDUX/$RUN2D`` of the previous RUN2D version will supply the list of manually edited files
        
    .. code-block:: shell
    
        grep -i "manual T   # Manually edited plan file (T: True, F: False)" */spPlan2d* fields/*/*/spPlan2d*
    

build fibermap files (*optional but recommended*)
"""""""""""""""""""""""""""""""""""""""""""""""""
During observations confSummary (FPS) or plPlugMapM (plates) files are created that map the targets to the BOSS spectrograph fibers. This step (:ref:`readfibermaps<readfibermaps>`; :ref:`slurm_readfibermap<slurm_readfibermap>`) reads the two types of files, converts them to a uniform format, and adds additional meta data (either from files or the SDSS5 internal database). If this step is skipped then they will be built with redux from uubatchpbs, but due to db connection limits it works better to prerun them.

.. note::
    if you are rerunning the pipeline of released data, you are encouraged to download the spfibermap files from the released data as this step uses the internal SDSSV Database.

.. code-block:: shell

    slurm_readfibermap --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --ppn 32 --apo
    slurm_readfibermap --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --ppn 32 --lco

build spTraceTab files
""""""""""""""""""""""
In the FPS operations era of SDSSV, a large emphasis was put on minimizing overheads. As part of this effort, the number of calibration frames has been reduced. In order to ensure proper tracing of the spectra, in light of observered flexure, the arc frames taken with each field are correlated with the arcs taken concurrently with trace flats at the start of evening observations. This step (:ref:`slurm_spTrace<slurm_spTrace>`) builds plan files of the calibration frames, traces the flat, and then builds the trace table (spTraceTab) files that are used by the pipelines inplace of the raw flat traces.

.. code-block:: shell
    
    spplan_trace --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --mjd_plans --mjdstart 59560 --apo --logfile apo_trace_plan.log
    spplan_trace --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --mjd_plans --mjdstart 60187 --lco --logfile lco_trace_plan.log

    slurm_spTrace --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --mjdstart 59560 --apo --skip_plan
    slurm_spTrace --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --mjdstart 60187 --lco --skip_plan

Run Daily Coadd
"""""""""""""""
This step (:ref:`uubatchpbs<uubatchpbs>`) takes the plan files built by :ref:`spplan<spplan>` and builds the redux-field-mjd script files. It then (if running at Utah) submits these redux-field-mjd scripts to the slurm queue. These scripts produce all of the field-mjd files.

.. code-block:: shell

    uubatchpbs --sdssv --obs apo --walltime "335:00:00" --nodes 7 --ppn 64 --merge3d
    uubatchpbs --sdssv --obs lco --walltime "335:00:00" --nodes 7 --ppn 64 --merge3d

Build Daily Summary Files
"""""""""""""""""""""""""
The final step of the pipeline is to take the individual field-mjd summary files and merge them in to final summary files.

.. code-block:: shell

    slurm_Summary --module bhm/|idlspec2d_version| --full --merge_only --walltime "335:00:00"

Field Epoch Coadds
^^^^^^^^^^^^^^^^^^
build the spplan files
""""""""""""""""""""""
Due to the nature of scheduling, weather, and engineering constraints, epochs are often split over several nights. In plate operations the epochs were typically defined by a single plugging of plate (with some additional constrains for specialed fields such as the RM fields). In the FPS operations, the epochs have a much more complicated defintion that rely on the cadences defined in the SDSSV internal database, where an epoch is defined by a set of designs with a defined maximum length between the first and last design. This step (:ref:`spplan_epoch<spplan_epoch>`) determines the exposures within an epoch and builds a plan file detailing the exposures to combine.

.. note::
    If you are rerunning the pipeline of released data, you are encouraged to download the spPlancombepoch files from the released data as this step uses the internal SDSSV Database.
    
.. code-block:: shell

    spplan_epoch --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --sdssv --apo --abandoned --logfile apo_epoch.log
    spplan_epoch --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --sdssv --lco --abandoned --logfile lco_epoch.log
    
.. note::
    If the run is being done for an IPL/DR Freeze include the "--started" flag to include epochs that have been started but not completed

.. code-block:: shell

    spplan_epoch --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --sdssv --apo --abandoned --logfile apo_epoch.log --started
    spplan_epoch --topdir $BOSS_SPECTRO_REDUX --run2d $RUN2D --sdssv --lco --abandoned --logfile lco_epoch.log --started


Run the epoch Coadd
"""""""""""""""""""
This step (:ref:`uubatchpbs<uubatchpbs>`) takes the plan files built by :ref:`spplan_epoch<spplan_epoch>` and builds the redux-field-mjd script files. It then (if running at Utah) submits these redux-field-mjd scripts to the slurm queue. These scripts produce all of the field-mjd files. The biggest difference between this and the daily version, is that the epoch redux scripts skip the initial extraction and calibration of the individual frames and uses those produced by the daily reduction.

.. code-block:: shell

    uubatchpbs --sdssv --walltime "335:00:00" --epoch --obs lco  --nodes 5 --ppn 64
    uubatchpbs --sdssv --walltime "335:00:00" --epoch --obs apo  --nodes 5 --ppn 64

Build Epoch Summary Files
"""""""""""""""""""""""""
The final step of the epoch pipeline is to take the individual field-mjd epoch summary files and merge them in to final summary files.
 
.. code-block:: shell

    slurm_Summary --module bhm/|idlspec2d_version| --full --epoch --merge_only --walltime "335:00:00"

Custom Coadds (eg. "allepoch")
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In SDSSV the variety of science programs (often sharing the same designs) required the development of an addition type of coadded spectra. In DR18 (v6_0_4), an early implementation of this was produced for the eFeds plates, where all of these plates, irregardless of the field centers and mjd were coadded together by running them through the full pipeline. In v6_1_1+ this strategy received a significant overhaul. Instead of reprocessing full plates from the start, we focus the effort on individual targets matching certain criteria and use the intermediate daily *spSpec* files (which includes the coadds of each each target on an exposure level) and produces an analogous to the spField files called spFullsky (because the original fields are not maintained). These files are then run through the 1d analysis and post production steps.

Managing the schema
"""""""""""""""""""
This step (:ref:`manage_coadd_Schema<manage_coadd_Schema>`) is to build the coadd schema model for the custom coadds.

.. code-block:: shell

    manage_coadd_Schema --name allepoch --DR  -r  -c  '*spiders*' '*bhm_gua*' '*bhm_csc*' '*mwm_erosita*' '*bhm_colr_galaxies*' -a

build the spplan files
""""""""""""""""""""""
Due to the nature of the FPS field designs, and the different requirements of different science programs, some coadding is needed on a target level.  This step (:ref:`spplan_target<spplan_target>`) uses the daily run summary file to determine the field and mjds of all observations of the selected targets, with the targets and cadences defined by the schema files (see :ref:`manage_coadd_Schema<manage_coadd_Schema>`). It then builds the a target level plan file. The coadded "MJD" is defined as the final observed MJD of each target and targets with the same "MJD" are grouped together for processing and analysis. If a "MJD" has less then 10 targets, they are grouped with the next largest MJD for operational efficiency.

.. code-block:: shell

    spplan_target --batch --DR --logfile lco_target_coadd_60280.log --lco
    spplan_target --batch --DR --logfile apo_target_coadd_60280.log --apo

Build the spFullSky files
"""""""""""""""""""""""""
This step (:ref:`uubatchpbs<uubatchpbs>`), similarly to the daily and epoch coadds, produces the redux script files and runs them. However, for the Custom Coadds, it initially only produces the spFullSky files, with the remaining steps run in the next step.

.. code-block:: shell

    uubatchpbs --sdssv --obs lco --nodes 1 --custom allepoch --allsky --coadd_only
    uubatchpbs --sdssv --obs apo --nodes 1 --custom allepoch --allsky --coadd_only

run 1d analysis and post processing steps
"""""""""""""""""""""""""""""""""""""""""
This step (:ref:`uubatchpbs<uubatchpbs>`), produces the redux script files and runs them for the 1D analysis and post processing steps.

.. code-block:: shell

    uubatchpbs --sdssv --obs lco --nodes 2 --custom allepoch --allsky --1dpost
    uubatchpbs --sdssv --obs apo --nodes 1 --custom allepoch --allsky --1dpost

Build Custom Coadd Summary Files
""""""""""""""""""""""""""""""""
The final step of the epoch pipeline is to take the individual Custom Coadded MJD summary files and merge them in to final summary files.

.. code-block:: shell

    slurm_Summary --module bhm/|idlspec2d_version| --full --custom allepoch --merge_only --walltime "335:00:00"
