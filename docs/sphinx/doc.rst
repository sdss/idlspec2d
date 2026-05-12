:tocdepth: 5

.. highlight:: none

Full Command Documention
========================

Documented below are the primary commands used to run the BOSS Data Reduction Pipeline. However, there are numerous other routines included in this package, which are called by these commands and have their own internal documentation.The legacy CLI interface is still included, documneted on :doc:`Legacy CLI<doc_legacy>`

Full Python Command Usage
-------------------------

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp:

boss_drp
^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_py:

.. code-block:: text

   Usage: boss_drp [OPTIONS] COMMAND [ARGS]...
   
   Options:
     --fullhelp  Show full help including all subcommands
     --help      Show this message and exit.
   
   Commands:
     version  Prints the IDLspec2D BOSS_DRP version
     config   BOSS DRP Config Commands
     plan     BOSS DRP Plan Commands
     daily    Run Pipeline Planning to Post
     batch    BOSS DRP Batch Cluster Tools
     run      Miscellaneous BOSS DRP Run Steps
     log      BOSS Pipeline Status Log
     clean    Tools to clean the BOSS DRP outputs
     tools    Miscellaneous BOSS DRP Tools

.. _boss_drp_version_py:

.. admonition:: boss_drp version
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp version [OPTIONS]
      
        Prints the IDLspec2D BOSS_DRP version
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.

.. _boss_drp_config_py:

.. admonition:: boss_drp config
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp config [OPTIONS]
      
        BOSS DRP Config Commands
      
      Options:
        --load_file TEXT      Config File to Load as the base
        --load_name TEXT      Config Name to Load as the base (eg. boss_drp,queue )
        --queue               Set this flag if you are loading a queue config
        --show                Print the Config
        --save TEXT           Location to Save the config to (required if you are
                              editing)
        -s, --save_to_envvar  Save to the ENVVAR (BOSS_DRP_PIPE_CONFIG_PATH or
                              BOSS_DRP_PIPE_CONFIG_PATH) if set
        --edit                Edit the Config
        --help                Show this message and exit.

boss_drp plan
"""""""""""""

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_plan_py:

.. admonition:: boss_drp plan
   :collapsible: open

   .. code-block:: text

      Usage: boss_drp plan [OPTIONS] COMMAND [ARGS]...
      
        BOSS DRP Plan Commands
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.
      
      Commands:
        daily        Produce the spPlan2d and spPlancomb files for the pipeline...
        trace        Produces spPlanTrace for the Use of Master Arc and Flat...
        epoch        Builds the spPlancombepoch files for the Epoch Coadd...
        CoaddSchema  Manage SDSSID/Catalogid Custom Coadds Schema
        target       Build SDSSID/CatalogID Custom Combine Plan

.. _boss_drp_plan_daily_py:

.. admonition:: boss_drp plan daily
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp plan daily [OPTIONS]
      
        Produce the spPlan2d and spPlancomb files for the pipeline run
      
      Options:
        --show_config                   Show the final config with out excecuting
                                        commands
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --skip1d / --no-skip1d          Skip spplan1d
        --skip2d / --no-skip2d          Skip spplan2d
        --override_manual / --no-override_manual
                                        Override/clobber manually edited plan
        -c, --clobber                   overwrites previous plan file
        --verbose TEXT                  Provide information about nonutlized frames
        --logfile TEXT                  Optional logfile (Including path)
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}
        --run2d TEXT                    Run2d to environmental variable
        --topdir TEXT                   Base run2d directory to BOSS_SPECTRO_REDUX
                                        environmental variable
        --remote / --no-remote          allow for remote access to data using sdss-
                                        access
        --release TEXT                  sdss_access data release
        --v_targ TEXT                   SDSS-V MOS Targeting Product Version (for no
                                        Database access use)
        --dither / --no-dither          Include Dither fields
        --commissioning / --no-commissioning
                                        Include SDSS-V FPS Commission Fields
        --sdssv / --no-sdssv            Include both SDSS-V Fields & Plates
        --fps / --no-fps                Include FPS Fields
        --plates / --no-plates          Include SDSS-V plates
        --legacy / --no-legacy          Include legacy (BOSS/eBOSS) plates
        --fieldend TEXT                 Ending Field
        --fieldstart TEXT               Starting Field
        --field TEXT                    Use data from these fields.
        --mjdend TEXT                   Ending MJD
        --mjdstart TEXT                 Starting MJD
        --mjd TEXT                      Use data from these MJDs.
        --manual_noarc / --no-manual_noarc
                                        if nomatched_arcs is False, builds spplan
                                        with unmatched arcs and mark as manual
        --multiple_arc / --no-multiple_arc
                                        Find all possible arc calibration frames
        --multiple_flat / --no-multiple_flat
                                        Find all possible flat calibration frames
        --minexp INTEGER                Min Science Exposures in Plan (default=1)
        --matched_arcs / --no-matched_arcs
                                        Allow Arc from another field/plate
        --matched_flats / --no-matched_flats
                                        Require Flat from a field/plate
        --quick / --no-quick            Use the list of new spPlan2d as a filter for
                                        fields
        --plate_epoch / --no-plate_epoch
                                        Use a variable max epoch length for plate
                                        coadd
        --help                          Show this message and exit.

.. _boss_drp_plan_trace_py:

.. admonition:: boss_drp plan trace
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp plan trace [OPTIONS]
      
        Produces spPlanTrace for the Use of Master Arc and Flat Frames to build
        Traces
      
      Options:
        --show_config                   Show the final config with out excecuting
                                        commands
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --override_manual / --no-override_manual
                                        Override/clobber manually edited plan
        -c, --clobber                   overwrites previous plan file
        --verbose TEXT                  Provide information about nonutlized frames
        --logfile TEXT                  Optional logfile (Including path)
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}
        --run2d TEXT                    Run2d to environmental variable
        --topdir TEXT                   Base run2d directory to BOSS_SPECTRO_REDUX
                                        environmental variable
        --remote / --no-remote          allow for remote access to data using sdss-
                                        access
        --release TEXT                  sdss_access data release
        --mjd_plans / --no-mjd_plans    Only build plans for MJDs with spPlan2d
        --mjdend TEXT                   Ending MJD
        --mjdstart TEXT                 Starting MJD
        --mjd TEXT                      Use data from these MJDs.
        --help                          Show this message and exit.

.. _boss_drp_plan_epoch_py:

.. admonition:: boss_drp plan epoch
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp plan epoch [OPTIONS]
      
        Builds the spPlancombepoch files for the Epoch Coadd Pipeline Runs
      
      Options:
        --show_config                   Show the final config with out excecuting
                                        commands
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --override_manual / --no-override_manual
                                        Override/clobber manually edited plan
        -c, --clobber                   overwrites previous plan file
        --logfile TEXT                  Optional logfile (Including path)
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}
        --run1d TEXT                    Run1d to environmental variable
        --run2d TEXT                    Run2d to environmental variable
        --topdir TEXT                   Base run2d directory to BOSS_SPECTRO_REDUX
                                        environmental variable
        --remote / --no-remote          allow for remote access to data using sdss-
                                        access
        --release TEXT                  sdss_access data release
        --v_targ TEXT                   SDSS-V MOS Targeting Product Version (for no
                                        Database access use)
        --sdssv / --no-sdssv            Include both SDSS-V Fields & Plates
        --fps / --no-fps                Include FPS Fields
        --plates / --no-plates          Include SDSS-V plates
        --legacy / --no-legacy          Include legacy (BOSS/eBOSS) plates
        --fieldend TEXT                 Ending Field
        --fieldstart TEXT               Starting Field
        --field TEXT                    Use data from these fields.
        --mjdend TEXT                   Ending MJD
        --mjdstart TEXT                 Starting MJD
        --mjd TEXT                      Use data from these MJDs.
        --min_epoch_len INTEGER         minimum length of epoch required to produce
                                        plan
        --started / --no-started        Create plans for started epochs (including
                                        unfinished)
        --abandoned / --no-abandoned    Create plans for abandoned epochs
        --minexp INTEGER                Min Science Exposures in Plan (default=1)
        --help                          Show this message and exit.

.. _boss_drp_plan_CoaddSchema_py:

.. admonition:: boss_drp plan CoaddSchema
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp plan CoaddSchema [OPTIONS]
      
        Manage SDSSID/Catalogid Custom Coadds Schema
      
      Options:
        -f, --coaddfile TEXT  File to store Coadding Schema (Default:
                              {topdir}/{run2d}/fields/SDSSV_BHM_COADDS.par)
        --topdir TEXT         Override value for the environment variable
                              $BOSS_SPECTRO_REDUX.
        --run2d TEXT          Override value for the environment variable $RUN2D
        --name TEXT           Name of Custom Coadd
        --DR                  DR/IPL Coadding
        -r, --rerun1d         Provides flag for coadd to be rerun though 1D analysis
        -a, --active          Activate (or deactivate) a Coadding Schema
        -c, --carton TEXT     list of cartons
        -i, --SDSSIDS TEXT    list of SDSS_IDS (or CatalogIDs if use_catid is set)
        -p, --program TEXT    list of programs
        -l, --legacy TEXT     list of Legacy Tags to include
        -u, --use_catid       Use CatalogIDs rather then SDSS_IDs
        --use_firstcarton     Use Firstcarton only for carton match (dont look at
                              db)
        -t, --cadence FLOAT   Number of days between coadd epochs
        -s, --show            Show Configurations
        --mjd TEXT            Use data from these MJDs.
        --help                Show this message and exit.

.. _boss_drp_plan_target_py:

.. admonition:: boss_drp plan target
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp plan target [OPTIONS]
      
        Build SDSSID/CatalogID Custom Combine Plan
      
      Options:
        --show_config                   Show the final config with out excecuting
                                        commands
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        -c, --clobber                   overwrites previous plan file
        --logfile TEXT                  Optional logfile (Including path)
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}
        --run1d TEXT                    Run1d to environmental variable
        --run2d TEXT                    Run2d to environmental variable
        --topdir TEXT                   Base run2d directory to BOSS_SPECTRO_REDUX
                                        environmental variable
        --mjdend TEXT                   Ending MJD
        --mjdstart TEXT                 Starting MJD
        --mjd TEXT                      Use data from these MJDs.
        --help                          Show this message and exit.

.. _boss_drp_daily_py:

.. admonition:: boss_drp daily
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp daily [OPTIONS]
      
        Plan, run Spectro-2D and Spectro-1D reductions, and run post pipeline steps
      
      Options:
        --module TEXT                   Module for daily run
        --apo                           Run apo
        --lco                           Run lco
        --mjd INTEGER                   Manually run for a single/list of mjd (does
                                        not update nextmjd.par)
        --range_mjd TEXT                Manually run for a range of mjds (does not
                                        update nextmjd.par)
        --dither / --no-dither          Skip Dither Engineering Fields
        --epoch / --no-epoch            Run Epoch Coadds
        --traceflat / --no-traceflat    Skip Building and using TraceFlats
        --force-arc2trace / --no-force-arc2trace
                                        Force Use of Arc2Trace
        --fibermap / --no-fibermap      Skip Pre-Run of readfibermap
        --no-prep / --prep              Skip building TraceFlats and spfibermaps
                                        before pipeline run
        --skip-plan [pipe|trace|all]    Skip the given plan
        --skip-plan-all                 Skip all plans
        --clobber [spplans|fibermap|trace|pipe|all]
                                        Clobber uubatchpbs + a combo of spPlan,
                                        fibermap, and TraceFlat run
        --clobber-all                   Clobber everything
        --healpix / --no-healpix        Turn off copy to healpix
        --summary / --no-summary        Build Summary Files
        --saveraw / --no-saveraw        save sdssproc outputs
        --debug / --no-debug            save extraction debug files
        --tagged                        sets --merge3d --sc tagged_daily --no-dither
                                        --monitor --allemail
        --daily                         sets --merge3d --sc fast_daily --monitor
                                        --allemail --no-healpix
        --dev                           sets --merge3d --sc tagged_daily --no-dither
                                        --monitor --no-healpix
        --topdir, --top-dir TEXT        Optional override value for the config
        --run1d TEXT                    Optional override value for config
        --run2d TEXT                    Optional override value for the config
        --nodist / --dist               unsets/sets --nodist and
                                        reactivates/deactivates the flux distortion
                                        corrections
        --bay15                         Set map3d to bayestar15 model
        --merge3d                       Set map3d to best 3d model
        --batch / --no-batch            run for multiple mjds in a single batch
        --nodb / --db                   skip Database operations
        --monitor / --no-monitor        Monitors pipeline status
        --pause INTEGER                 Pause time (s) in status updates
        --allemail / --no-allemail      Email intermediate log using all emails in
                                        $DAILY_DIR/etc/emails (defaults to first
                                        email only)
        --pipe_config, --pipe-config, --pc TEXT
                                        Queue Config name
        --pipe_config_file, --pipe-config-file, --pcf TEXT
                                        Queue Config File Path
        --queue_config, --queue-config, --qc TEXT
                                        Queue Config name
        --queue_config_file, --queue-config-file, --qcf TEXT
                                        Queue Config File Path
        --no-write / --write            skip writing and submitting job
        --nosubmit / --submit           Skip submitting uubatch job (ideal for
                                        allowing editting of plans)
        --walltime TEXT                 Wall time in hours
        --mem_per_cpu, --mem-per-cpu TEXT
                                        Memory allocated per CPU
        --nbundle INTEGER               Number of jobs to bundle
        --show_config                   Show the final config with out excecuting
                                        commands
        --help                          Show this message and exit.

boss_drp batch
""""""""""""""

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_batch_py:

.. admonition:: boss_drp batch
   :collapsible: open

   .. code-block:: text

      Usage: boss_drp batch [OPTIONS] COMMAND [ARGS]...
      
        BOSS DRP Batch Cluster Tools
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.
      
      Commands:
        Summary       Create daily field merge queue job
        pipe          Build idlspec2d redux and submit to the cluster queue.
        readfibermap  Create a batch readfibermap job.
        runfix        Check for failed runs and setup the runs to clean and...
        sos           Create SOS queue job.
        spTrace       Create spTrace Queue jobs

.. _boss_drp_batch_Summary_py:

.. admonition:: boss_drp batch Summary
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp batch Summary [OPTIONS]
      
        Create daily field merge queue job
      
      Options:
        -m, --module TEXT               module file to use (ex bhm/master[default]
                                        or bhm/v6_0_9)
        --topdir TEXT                   Optional override value for the config
        --run2d TEXT                    Optional override value for the config
        --run1d TEXT                    Optional override value for the config
        --epoch / --no-epoch            Run for epoch Coadds
        --custom TEXT                   Run for epoch Coadds
        --daily / --no-daily            only run if daily run has been run today
        --monitor / --no-monitor        Monitor job and send email at completion
                                        with the logs
        --fieldlist / --no-fieldlist    Running Fieldlist
        --backup INTEGER                Number of backups to keep, or None (or 0) to
                                        not create backup
        --n_iter INTEGER                number of iterations of field merge to run
        --ndays INTEGER                 Limit spAll update to last ndays
        --skip_specprimary              Skip calculation of Specprimary
        --update_specprimary            Only update new Specprimary
        --utah_daily / --no_utah        Load tagged daily run into
                                        Pipelines.boss_drp database table
        --verbose / --no-verbose        Run Fieldmerge with verbose
        --email_start / --no-email_start
                                        Send email at start of run
        --defaults                      Sets --merge_only  --backup 3  --monitor
                                        --update_specprimary --ndays 10 --qc
                                        summary_full
        --show_config                   Show the final config with out excecuting
                                        commands
        --queue_config_file, --qcf TEXT
                                        Queue Config File Path
        --queue_config, --qc TEXT       Queue Config name
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --walltime TEXT                 Wall time in hours
        --mem TEXT                      memory in bytes
        --ppn INTEGER                   Number of processors per node
        --nosubmit / --submit           Create queue job but do not submit it
        --help                          Show this message and exit.

.. _boss_drp_batch_pipe_py:

.. admonition:: boss_drp batch pipe
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp batch pipe [OPTIONS]
      
        Build idlspec2d redux and submit to the cluster queue.  Without access to
        the SDSS Slurm package, it prints the commands for manual execution
      
      Options:
        --allemail / --primary-email    Email intermediate log using all emails in
                                        $DAILY_DIR/etc/emails (defaults to first
                                        email only)
        --email / --no-email            Email log using $DAILY_DIR/etc/emails
        --1dpost / --all                Run 1d analysis and post processing only
        --coadd-only / --no-coadd-only  Run spspec_target_merge only
        --single-mjd / --no-single-mjd  Run Each Custom MJD coadd+1dpost as seperate
                                        job
        --allsky / --no-allsky          All Sky Coadds
        --custom TEXT                   Name of custom Coadd Schema
        --epoch / --no-epoch            Epoch Coadds
        --nbundle INTEGER               Number of jobs to bundle
        --nosubmit / --submit           Build, but not submit redux files
        --ppn INTEGER                   Number of processors per node
        --nodes INTEGER                 Number of Nodes
        --walltime TEXT                 Wall time in hours
        --mem-per-cpu TEXT              Memory allocated per CPU
        --no-write / --write            skip writing and submitting job
        --show_config                   Show the final config with out excecuting
                                        commands
        --queue_config_file, --qcf TEXT
                                        Queue Config File Path
        --queue_config, --qc TEXT       Queue Config name
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --mjdend INTEGER                Ending MJD
        --mjdstart INTEGER              Starting MJD
        --mjd INTEGER                   MJD dates to reduce; default="*"
        --fieldend TEXT                 End Field/Plate number
        --fieldstart TEXT               Starting Field/Plate number
        -f, --field TEXT                Plate/Field numbers to reduce default="*"
        --only1d / --no-only1d          Run spec1d step only (eg. spreduce1d_empca,
                                        XCSAO)
        --skip2d / --no-skip2d          Skip spreduce2d
        --no-healpix, --nohp / --healpix, --hp
                                        Turn off copy to healpix
        --clobber / --no-clobber        Clobber redux
        --a2t / --no-a2t                Force Use of Arc2Trace
        --v_targ TEXT                   SDSS-V MOS Targeting Product Version  (for
                                        no Database access use)
        --release TEXT                  sdss_access data release ...
        --fast_no_db TEXT               When using --no-db, streamlines process and
                                        only gets parallax from MOS target files
        --no-db / --db                  skip Database operations
        --debug / --no-debug            Save extraction debug files
        --saveraw / --no-saveraw        Save sdssproc outputs
        --fibermap_clobber / --no-fibermap_clobber
                                        Clobber spfibermap fits file
        --onestep_coadd / --no-onestep_coadd
                                        Use legacy one step version of coadd
        --update_specprimary            Only update new Specprimary
        --skip_specprimary              Skip Calculation of Specprimary
        --nodist / --dist               Unset --nodist and reactivate the flux
                                        distortion corrections
        --noxcsao / --xcsao             Skip pyXCSAO
        --map3d [bayestar15|bay15|merge3d]
                                        Name of 3d dustmap to use with MWM_fluxer
                                        (default=None)
        --MWM-fluxer, --mwm / --no-MWM-fluxer, --no-mwm
        --no-reject / --reject          Deactivate Rejection in Coadd
        --idlutils_1d TEXT              idlutils override version of spec1d
        --run2d TEXT                    Optional override value for the config
        --run1d TEXT                    Optional override value for the config
        --topdir TEXT                   Optional override value for the config
        --merge3d                       Set map3d to best 3d model
        --bay15                         Set map3d to bayestar15 model
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}. Can be repeated.
        --sdssv                         --mwm --no-reject --merge3d
        --help                          Show this message and exit.

.. _boss_drp_batch_readfibermap_py:

.. admonition:: boss_drp batch readfibermap
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp batch readfibermap [OPTIONS]
      
        Create a batch readfibermap job. Without access to the SDSS Slurm package,
        it prints the commands for manual execution
      
      Options:
        --topdir TEXT                   Boss Spectro Redux base directory
        --run2d TEXT                    Run2d
        --clobber / --no-clobber        Clobber spfibermaps
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}. Can be repeated.
        --v_targ TEXT                   SDSS-V MOS Targeting Product Version (for no
                                        Database access use)
        --mjdend INTEGER                Ending MJD
        --mjdstart INTEGER              Starting MJD
        --mjd INTEGER                   MJD dates to reduce; default="*"
        --show_config                   Show the final config with out excecuting
                                        commands
        --queue_config_file, --qcf TEXT
                                        Queue Config File Path
        --queue_config, --qc TEXT       Queue Config name
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --nbundle INTEGER               Number of jobs to bundle
        --ppn INTEGER                   Number of processors per node
        --nodes INTEGER                 Number of Nodes
        --walltime TEXT                 Wall time in hours
        --mem_per_cpu TEXT              Memory allocated per CPU
        --help                          Show this message and exit.

.. _boss_drp_batch_runfix_py:

.. admonition:: boss_drp batch runfix
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp batch runfix [OPTIONS]
      
        Check for failed runs and setup the runs to clean and rerun the crashed
        field-mjds
      
      Options:
        --full / --no-full              Rerun full pipeline regardless of crashed
                                        step
        --running / --no-running        Select Field-MJDs with running status
        --topdir TEXT                   Optional override value for the config
        --run2d TEXT                    Optional override value for the config
        --run1d TEXT                    Optional override value for the config
        --epoch / --daily               Run for epoch Coadds
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}. Can be repeated.
        --mjdend INTEGER                Ending MJD
        --mjdstart INTEGER              Starting MJD
        --mjd INTEGER                   MJD dates to reduce; default="*"
        --show_config                   Show the final config with out excecuting
                                        commands
        --queue_config_file, --qcf TEXT
                                        Queue Config File Path
        --queue_config, --qc TEXT       Queue Config name
        --pipe_config_file, --pcf TEXT  Queue Config File Path
        --pipe_config, --pc TEXT        Queue Config name
        --nbundle INTEGER               Number of jobs to bundle
        --ppn INTEGER                   Number of processors per node
        --nodes INTEGER                 Number of Nodes
        --walltime TEXT                 Wall time in hours
        --mem_per_cpu TEXT              Memory allocated per CPU
        --no-write / --write            skip writing and submitting job
        --nosubmit / --submit           Build, but not submit redux files
        --help                          Show this message and exit.

.. _boss_drp_batch_sos_py:

.. admonition:: boss_drp batch sos
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp batch sos [OPTIONS]
      
        Create SOS queue job. Without access to the SDSS Slurm package, it prints
        the commands for manual execution
      
      Options:
        --lco                           Run lco
        --apo                           Run apo
        --obs [apo|lco]                 Observatory {apo,lco}. Can be repeated.
        --mjdend INTEGER                Ending MJD
        --mjdstart INTEGER              Starting MJD
        --mjd INTEGER                   MJD dates to reduce; default="*"
        --no-reject                     Overrides the Calibration rejection (use
                                        with caution)
        -f, --clobber_fibermap          Clobbers the existing spfibermap files
        -n, --no-arc2trace              Skip Utilizing arc2trace refinements
        -o, --forcea2t                  Force arc2trace for all fields (even if flat
                                        exists for field)
        --sdssv-sn2 / --no-sdssv-sn2    Report a second set of SN2 values with
                                        updated fit parameters
        --sn2-15 / --no-sn2-15          Report a set of SN2 values with a fiducial
                                        mag of 15
        --bright                        Display BOSS_only Bright Time Operation
                                        SN2_15
        --show_config                   Show the final config with out excecuting
                                        commands
        --queue_config_file, --qcf TEXT
                                        Queue Config File Path
        --queue_config, --qc TEXT       Queue Config name
        --nbundle INTEGER               Number of jobs to bundle
        --ppn INTEGER                   Number of processors per node
        --nodes INTEGER                 Number of Nodes
        --walltime TEXT                 Wall time in hours
        --mem_per_cpu TEXT              Memory allocated per CPU
        --no-submit                     Build, but not submit redux files
        --help                          Show this message and exit.

.. _boss_drp_batch_spTrace_py:

.. admonition:: boss_drp batch spTrace
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp batch spTrace [OPTIONS]
      
      Options:
        --topdir TEXT                 Optional override value for the config
        --run2d TEXT                  Optional override value for the config
        --lco                         Run lco
        --apo                         Run apo
        --obs [apo|lco]               Observatory {apo,lco}. Can be repeated.
        --clobber / --no-clobber      Clobber the existing Plan files
        --debug / --no-debug          Save sdssproc outputs
        --saveraw / --no-saveraw      Clobber the existing Plan files
        --skip_plan / --no-skip_plan  Skip creating plans and use currently existing
                                      plans
        --mjdend INTEGER              Ending MJD
        --help                        Show this message and exit.

boss_drp run
""""""""""""

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_run_py:

.. admonition:: boss_drp run
   :collapsible: open

   .. code-block:: text

      Usage: boss_drp run [OPTIONS] COMMAND [ARGS]...
      
        Miscellaneous BOSS DRP Run Steps
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.
      
      Commands:
        readfibermap         Produces spfibermap file corresponding to a...
        boss_arcs_to_traces  Routine to transfer trace locations from an...
        run_PyXCSAO          Runs pyXCSAO to determine RVs
        fieldlist            Build/load BOSS Fieldlist
        fieldmerge           Build BOSS spAll Summary Files
        spSpec_reformat      Build Spec Files
        spcalib_qa           Compare photometric accuracy of standards
        update_flags         Update SDSSV Targeting flats in the summary files
        Plot_QA              Plot the SpectroPhotometry and SN2 QA plots

.. _boss_drp_run_readfibermap_py:

.. admonition:: boss_drp run readfibermap
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run readfibermap [OPTIONS]
      
        Produces spfibermap file corresponding to a spplan2d (or single confSummary
        file for SOS).
      
      Options:
        -p, --spplan2d TEXT  spplan2d file for idlspec2d run
        --topdir TEXT        Alternative output directory (defaults to location of
                             spplan2d file or /data/boss/sos/{mjd} for SOS)
        -c, --clobber        Overwrites previous spfibermap file
        --fast               When using --no-db, streamlines process and only gets
                             parallax from MOS target files
        --datamodel TEXT     Supply a datamodel file (defaults to
                             $IDLSPEC2D/datamodel/spfibermap_dm.par or
                             $IDLSPEC2D/datamodel/spfibermap_sos_dm.par for SOS)
        -s, --SOS            Produces spfibermap for SOS
        --release TEXT       sdss_access data release (defaults to sdsswork),
                             required if you do not have proprietary access
                             [default: sdsswork]
        --remote             Allow for remote access to data using sdss-access
        --v_targ TEXT        SDSS-V MOS Targeting Product Version (for no Database
                             access use)  [default: *]
        --confSummary TEXT   confSummary file for SOS (required with --SOS)
        --ccd [b2|r2|b1|r1]  CCD for SOS
        --mjd TEXT           MJD of observation
        --log                Creates log file in topdir
        --help               Show this message and exit.

.. _boss_drp_run_boss_arcs_to_traces_py:

.. admonition:: boss_drp run boss_arcs_to_traces
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run boss_arcs_to_traces [OPTIONS]
      
        Routine to transfer trace locations from an initial arc/flat pair given
        subsequent arc frames only
      
      Options:
        --mjd INTEGER      MJD to process  [required]
        --outdir TEXT      output directory
        --obs [lco|apo]    observatory  [default: lco]
        --vers TEXT        BOSS_SPECTRO_REDUX version  [default: master]
        --threads INTEGER  number of threads  [default: 8]
        --cams TEXT        Supply the camera for operation with SOS files
        --fitsname TEXT    Supply the FitsName for SOS error reporting
        --sosdir TEXT      Base SOS output directory
        --clobber          clobber?
        --no-hash          Skip updating the file hash
        --help             Show this message and exit.

.. _boss_drp_run_run_PyXCSAO_py:

.. admonition:: boss_drp run run_PyXCSAO
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run run_PyXCSAO [OPTIONS] FITSFILE
      
        Run PyXCSAO for a full spField FITS file using the phoenix_full1 template
        grid. The input file can be either a normal FITS or gzipped FITS file.
      
      Options:
        -r, --run1d TEXT  run1d name  [default: (env: RUN1D)]
        --epoch           Run for epoch coadds
        --custom TEXT     Name of custom coadd
        --help            Show this message and exit.

.. _boss_drp_run_fieldlist_py:

.. admonition:: boss_drp run fieldlist
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run fieldlist [OPTIONS]
      
        Build/load BOSS Fieldlist
      
      Options:
        -c, --create     Create Fieldlist
        --topdir TEXT    Optional override value for $BOSS_SPECTRO_REDUX  [default:
                         (env: BOSS_SPECTRO_REDUX)]
        --run1d TEXT     Optional override value for $RUN1D  [default: (env: RUN1D)]
        --run2d TEXT     Optional override value for $RUN2D  [default: (env: RUN2D)]
        --outdir TEXT    Optional output directory (defaults to topdir/$RUN2D)
        --skipcart TEXT  List of cartridges to skip
        --epoch          Produce FieldList for epoch coadds
        --abandoned      Produce FieldList for epoch coadds (including abondoned)
        --started        Produce FieldList for epoch coadds (including started)
        --basehtml TEXT  HTML path for figure (defaults relative to topdir)
        --logfile TEXT   Manually set logfile (including path)
        --debug          Print full python errors instead of simplified logger
                         messages
        --noplot         Skip updating the sky plots
        --help           Show this message and exit.

.. _boss_drp_run_fieldmerge_py:

.. admonition:: boss_drp run fieldmerge
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run fieldmerge [OPTIONS]
      
        Build BOSS spAll Summary Files
      
      Options:
        --run2d TEXT              Optional override value for the environment
                                  variable $RUN2D  [default: (env: RUN2D)]
        --indir TEXT              Optional override value for the environment
                                  variable $BOSS_SPECTRO_REDUX  [default: (env:
                                  BOSS_SPECTRO_REDUX)]
        --skip_line               Skip the generation of spAllLine.fits
        --include_bad             Include bad fields
        --legacy                  Include columns used by SDSS-IV and depreciated in
                                  SDSS-V
        --skip_specprimary        Skip creation of specprimary and associated
                                  columns
        --update_specprimary      Keep existing specprimary and associated columns
                                  and only update new row (and their secondaries)
        --lite                    Produce lite version of spAll file
        --include_XCSAO           Include XCSAO columns
        -f, --field TEXT          Run for a single Field
        -m, --mjd TEXT            Run for a single MJD
        --clobber                 Clobber all spAll-field-mjd files
        --bkup                    Backup existing spAll files
        --verbose                 Log columns not saved
        --logfile TEXT            Manually set logfile
        --epoch                   Produce spAll for epoch coadds
        --programs TEXT           List of programs to include. Repeat the option for
                                  multiple values.
        --datamodel TEXT          Supply a spAll datamodel file (defaults to
                                  $IDLSPEC2D/datamodel/spall_dm.par)
        --line_datamodel TEXT     Supply a spline datamodel file (defaults to
                                  $IDLSPEC2D/datamodel/spzline_dm.par)
        --outroot TEXT            Path and root of filename for output (defaults to
                                  spectra/full or summary)
        --allsky                  Build spAll for Allsky Custom Coadd
        --custom TEXT             Name of Custom Coadd
        --run1d TEXT              Optional override value for the environment
                                  variable $RUN1D (only for custom allsky coadds)
                                  [default: (env: RUN1D)]
        --ndays INTEGER           Limit update to last ndays
        --freeze_output           Freeze MJD limited parquet files
        --update_target_flags     Use the spTargeting file to update the summary
                                  file to the latest Targeting Flags
        --mjdstart INTEGER        Limit update to MJD on/after
        --mjdend INTEGER          Limit update to MJD on/before
        --MJD_dir TEXT            Location to save the MJD level temporary files
                                  (defaults to BOSS_SPECTRO_SCRATCH)
        --to_fits                 Dump Parquet to fits format
        --force, --force_rebuild  Rebuild Summary even if nothing changed
        --keep_active             Run "touch" on all intermediate files to keep them
                                  active
        --help                    Show this message and exit.

.. _boss_drp_run_spSpec_reformat_py:

.. admonition:: boss_drp run spSpec_reformat
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run spSpec_reformat [OPTIONS]
      
        Build Spec Files
      
      Options:
        -f, --field TEXT  Run for a single Field  [required]
        -m, --mjd TEXT    Run for a single MJD  [required]
        --topdir TEXT     Optional override value for the environment variable
                          $BOSS_SPECTRO_REDUX  [default: (env: BOSS_SPECTRO_REDUX)]
        --run2d TEXT      Optional override value for the environment variable
                          $RUN2D  [default: (env: RUN2D)]
        --run1d TEXT      Optional override value for the environment variable
                          $RUN1D  [default: (env: RUN1D)]
        --custom TEXT     Name of Custom Coadd schema
        -p, --plot        Create spec plots
        -e, --epoch       Run for epoch Coadds
        --lsdr10          Include Legacy Survey DR10 links on HTML
        --allsky          Reformat for Allsky Custom Coadd
        --help            Show this message and exit.

.. _boss_drp_run_spcalib_qa_py:

.. admonition:: boss_drp run spcalib_qa
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run spcalib_qa [OPTIONS]
      
        Compare photometric accuracy of standards
      
      Options:
        --run2d TEXT                    Optional override value for the enviro
                                        variable $RUN2D
        --bsr TEXT                      Optional override value for the enviro
                                        variable $BOSS_SPECTRO_REDUX
        -f, --field TEXT                Run for a single Field
        -m, --mjd TEXT                  Run for a single MJD
        -r, --rerun                     Rerun for all field-mjds in spAll
        -n, --nobkup                    Do not backup output and log file
        -e, --epoch                     run for epoch coadds
        -c, --catchup                   Run for missing field-mjds spAll
        --outdir TEXT                   Location to Save plots to (overrides the
                                        defaults)
        --run2d_alt TEXT                Alternative RUN2D for comparison
        -bsra, --boss_spectro_redux_alt TEXT
                                        Alternative BOSS_SPECTRO_REDUX for
                                        comparison
        --plot_only                     Create Plots but dont update summary file
        --help                          Show this message and exit.

.. _boss_drp_run_update_flags_py:

.. admonition:: boss_drp run update_flags
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run update_flags [OPTIONS]
      
        Update SDSSV Targeting flats in the summary files
      
      Options:
        --topdir TEXT  Optional override value for the environment variable
                       $BOSS_SPECTRO_REDUX  [default: (env: BOSS_SPECTRO_REDUX)]
        --run2d TEXT   Optional override value for the environment variable $RUN2D
                       [default: (env: RUN2D)]
        --custom TEXT  Name of Custom Coadd schema
        --clobber      Clobber spTargeting file
        --nobackup     Skip backup of existing summary files
        --help         Show this message and exit.

.. _boss_drp_run_Plot_QA_py:

.. admonition:: boss_drp run Plot_QA
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp run Plot_QA [OPTIONS]
      
        Plot the SpectroPhotometry and SN2 QA plots
      
      Options:
        -r, --run2d TEXT        List of run2ds  [default: (env: RUN2D)]
        --test TEXT             List of True/False test run2d (corresponding to
                                run2d)
        --test_path TEXT        test Run2d path modification  [default: /test/sean/]
        --mjds_low TEXT         List of mjd lower limits (use 'None' for no limit)
        --mjds_high TEXT        List of mjd upper limits (use 'None' for no limit)
        --clobber_lists         Clobber list of fieldIDs
        --lco / --apo           Flag for LCO vs APO
        --publish               Create publication version of plot
        --html                  Produces Plotly interactive HTML versions of the
                                plots
        --html_name TEXT        Name of HTML file (default = BOSS_QA-{obs}.html)
        -f, --fast_opsdb        Skips OpsDB queries for SOS SN2 (and uses cached if
                                available)
        -e, --epoch             Produce plots for epoch coadds
        -c, --cron              Produce cronlogs
        --fid, --fieldids TEXT  Limit to these FieldIDs
        --compare               Direct Comparison of the run2ds in list
        --start_mjd INTEGER     Limit to MJDs on or after
        --end_mjd INTEGER       Limit to MJDs on or before
        --output_dir TEXT       Overrides the output directory
        --help                  Show this message and exit.

.. _boss_drp_log_py:

.. admonition:: boss_drp log
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp log [OPTIONS]
      
        BOSS Pipeline Status Log
      
      Options:
        --obs [apo|lco]     Observatory for status update
        --apo               Run apo
        --lco               Run lco
        --mjd INTEGER       Update these MJDs
        --mjdstart INTEGER  Starting MJD
        --mjdend INTEGER    Ending MJD
        --epoch             Run for epoch Coadds
        --custom TEXT       Name of custom Coadd
        --topdir TEXT       Optional override value for the environment variable
                            $BOSS_SPECTRO_REDUX
        --run2d TEXT        Optional override value for the enviro variable $RUN2D
        --run1d TEXT        Optional override value for the enviro variable $RUN1D
        --email             Send each mjd status as email
        --fast              Skip updating index until end
        --refresh           Refresh all the existing Status logs for obs
        --refresh_error     Refresh existing Status logs for obs with errors
        --refresh_critical  Refresh existing Status logs for obs with critical
                            errors
        --force             Refresh Summaries pages
        --help              Show this message and exit.

boss_drp clean
""""""""""""""

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_clean_py:

.. admonition:: boss_drp clean
   :collapsible: open

   .. code-block:: text

      Usage: boss_drp clean [OPTIONS] COMMAND [ARGS]...
      
        Tools to clean the BOSS DRP outputs
      
      Options:
        --help  Show this message and exit.
      
      Commands:
        backups  Clean the Summary Table File Backups
        run      Clean pipeline prodcuts fro a given field, mjd, or field-mjd

.. _boss_drp_clean_backups_py:

.. admonition:: boss_drp clean backups
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp clean backups [OPTIONS]
      
        Clean the Summary Table File Backups
      
      Options:
        --topdir TEXT      Boss Spectro Redux base directory
        --run2d TEXT       Run2d
        --epoch            run for the epoch coadds
        --custom TEXT      Name of custom Coadd
        --backups INTEGER  Number of backups to keep
        --help             Show this message and exit.

.. _boss_drp_clean_run_py:

.. admonition:: boss_drp clean run
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp clean run [OPTIONS]
      
        Clean pipeline prodcuts fro a given field, mjd, or field-mjd
      
      Options:
        --clean_type, --clean [all|spec2d|comb|spec1d|post|merge|reformat|spcalib]
                                        Pipeline Step to start the cleaning
                                        [required]
        --topdir TEXT                   Optional override value for the environment
                                        variable $BOSS_SPECTRO_REDUX
        --run2d TEXT                    Optional override value for the environment
                                        variable $RUN2D
        --run1d TEXT                    Optional override value for the environment
                                        variable $RUN1D
        --epoch                         Clean up epoch run
        --reset                         if clean_type == all, then remove plans and
                                        redux
        --remove_redux                  if clean_type == all, then remove redux
        --dry                           Print Files to be removed rather then remove
        --verbose                       Print Files paths (with wildcards) to be
                                        removed
        -f, --field TEXT                Run for a single Field
        -m, --mjd TEXT                  Run for a single MJD
        --fmjd TEXT                     List of Field-MJDs to clean
        --help                          Show this message and exit.

boss_drp tools
""""""""""""""

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_tools_py:

.. admonition:: boss_drp tools
   :collapsible: open

   .. code-block:: text

      Usage: boss_drp tools [OPTIONS] COMMAND [ARGS]...
      
        Miscellaneous BOSS DRP Tools
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.
      
      Commands:
        sdR_hdrfix       Create the files used by the pipeline to fix the...
        flag_manual_cal  Build spManCal.par file to flag manual alternative...
        opFiber          Prints updated/refined values for opfibers using...

.. _boss_drp_tools_sdR_hdrfix_py:

.. admonition:: boss_drp tools sdR_hdrfix
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp tools sdR_hdrfix [OPTIONS] EXPID
      
        Create the files used by the pipeline to fix the header meta data of the
        BOSS exposures
      
      Options:
        -v, --value TEXT                updated header keyword value (required if
                                        key is set)
        -k, --key TEXT                  header keyword to update (required if value
                                        is set)
        --designid INTEGER              DesignID
        --confid INTEGER                ConfigureID
        --fieldid INTEGER               FieldID
        --cartid [FPS-S|FPS-N]          Cartridge Mounted
        --tai-beg FLOAT                 Starting time (tai) of exposure
        --exptime FLOAT                 Exposure length (s)
        --flavor [bias|dark|flat|arc|science|smear]
                                        Type/Flavor of exposure
        --quality [excellent|test|bad]  Set Quality flat of exposures
        --hartmann [out|right|left|closed]
                                        Hartmann Door Status
        --flat                          short cut to set FF = 1 1 1 1 & FFS =  1 1 1
                                        1 1 1 1 1
        --arc                           short cut to set all relevant arc lamps to 1
                                        1 1 1
        --HEAR [0|1]...                 HeAr arc Lamp
        --HGCD [0|1]...                 HeCd arc Lamp
        --NE [0|1]...                   Ne arc lamp
        --FFS [0|1]...                  Flat Field Screen
        --FF [0|1]...                   Flat Field Lamp
        -t, --test                      Flag as test
        -b, --bad                       Flag as bad
        --nogit                         Skip automatic git add
        -u, --no-update                 Skip updating SOS logs
        --cameras [b1|b2|r1|r2|??]      Cameras for hdr update
        --clobber                       Clobber sdHdrFix file
        --obs [apo|lco]                 Observatory  [required]
        -m, --mjd TEXT                  MJD of file (default: latest)
        --help                          Show this message and exit.

.. _boss_drp_tools_flag_manual_cal_py:

.. admonition:: boss_drp tools flag_manual_cal
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp tools flag_manual_cal [OPTIONS]
      
        Build spManCal.par file to flag manual alternative calibration frames for
        spPlan
      
      Options:
        -o, --observatory, --obs [apo|lco]
                                        Observatory  [required]
        -m, --mjd INTEGER               MJD  [required]
        -f, --field TEXT                FieldID  [required]
        -e, --expid INTEGER             Exposure ID to manually set the calibration
                                        frame exposure ID
        -t, --type [arc|flat]           Calibration Type  [required]
        --nogit                         Skip automatic git add
        --help                          Show this message and exit.

boss_drp tools opFiber
~~~~~~~~~~~~~~~~~~~~~~

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_drp_tools_opFiber_py:

.. admonition:: boss_drp tools opFiber
   :collapsible: open

   .. code-block:: text

      Usage: boss_drp tools opFiber [OPTIONS] COMMAND [ARGS]...
      
        Prints updated/refined values for opfibers using spFlat traces
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.
      
      Commands:
        Guess   Guess at the opFiberFPS parameters using a sdProc-XX-XXXXXXXX.fits
                file
        Refine  Refine the opFiberFPS parameters using a spFlat.

.. _boss_drp_tools_opFiber_Guess_py:

.. admonition:: boss_drp tools opFiber Guess
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp tools opFiber Guess [OPTIONS] SDPROCFILE
      
        Takes a processed image frame produced using the /sawraw flag in sdssproc
        (or indirecly via spreduce2d) as an input. It then uses either the
        bundlefiber list of number of fibers per bundle (supplied as input or via
        opFiberFPS) combined with the scipy peak finding algarithm to create a first
        guess of the peak fiber and bundle gaps. It uses the median flux of the 11
        central pixel (along the dispersion axis) to build the flux array
      
      Options:
        -b, --bundlefibers INTEGER  List of number of fibers per bundle. Use
                                    multiple times, for example: -b 2 -b 4 -b 4
        -m, --mjd INTEGER           MJD of new updated OpFiber Fiberparameter entry.
        -f, --plot                  Whether to plot the flux slice and detected
                                    peaks.
        --min_peak_sep FLOAT        Minimum separation between detected peaks.
                                    [default: 6]
        --min_peak_height FLOAT     Minimum flux level to be detected as a peak.
                                    [default: 5000]
        -p, --precision INTEGER     Precision of the reported fiberspacing and
                                    bundle gaps.  [default: 3]
        --help                      Show this message and exit.

.. _boss_drp_tools_opFiber_Refine_py:

.. admonition:: boss_drp tools opFiber Refine
   :collapsible: closed

   .. code-block:: text

      Usage: boss_drp tools opFiber Refine [OPTIONS] FITSFILE
      
        Refine the opFiberFPS parameters using a spFlat.
      
      Options:
        -p, --precision INTEGER  Precision of the reported fiberspacing and bundle
                                 gaps.  [default: 3]
        --help                   Show this message and exit.

.. _SOS:

SOS
^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _SOS_py:

.. code-block:: text

   Usage: SOS [OPTIONS] COMMAND [ARGS]...
   
     SOS process for reducing BOSS data on the Moutain
   
   Options:
     --fullhelp              Show full help including all subcommands
     --red                   Red Camera Process
     --blue                  Blue Camera Process
     --joint                 Both Camera Processes
     --catchup               Run Catchup on the night or (MJD)
     --redoMode              Save outputs of MJD or exposure to sosredo
     --test                  Save outputs and logs to sosredo/dev
     --unlock                Unlock Locked Files
     -e, --exp EXPID         exposure id (or range of exp id 500-510) (with or
                             without leading zeros)
     -m, --mjd MJD           MJD
     --nodb                  skip opsdb load
     --no-gz                 Overrides the requirement for '.gz' compressed files
                             (experimental)
     --no-reject             Overrides the Calibration rejection (use with
                             caution)
     -f, --clobber_fibermap  Clobbers the existing spfibermap files
     --no-sdssv-sn2          Report a second set of SN2 values with updated fit
                             parameters
     --no-sn2-15             Skip reporting a set of SN2 values with a fiducial
                             mag of 15 for engineering fields
     --bright                Display BOSS_only Bright Time Operation SN2_15 for
                             all fields
     -n, --no-arc2trace      Skip Utilizing arc2trace refinements
     -o, --forcea2t          Force arc2trace for all fields (even if flat exists
                             for field)
     --plot                  Produce Science Plots for each exposure
     -v, --verbose           prints the only (or red if joint) active SOS process
                             to terminal
     --help                  Show this message and exit.
   
   Commands:
     Tools  Tools to use with SOS

SOS Tools
"""""""""

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _SOS_Tools_py:

.. admonition:: SOS Tools
   :collapsible: open

   .. code-block:: text

      Usage: SOS Tools [OPTIONS] COMMAND [ARGS]...
      
        Tools to use with SOS
      
      Options:
        --fullhelp  Show full help including all subcommands
        --help      Show this message and exit.
      
      Commands:
        FiberQA              Create Fiber info Summary for SOS
        Log                  Build BOSS Exposure Log
        boss_arcs_to_traces  Routine to transfer trace locations from an...
        flag_manual_cal      Build spManCal.par file to flag manual alternative...
        hash                 Create or check the SOS file hash
        htmlIndex            Build sos/combined/index.html
        loadsn2              Load SOS SN2 values into OpsDB
        log2html             Create the HTML Logging Page for SOS
        parse_runtime        Process log file to calculate elapsed times for SOS
        plot                 Plot the Science frame for SOS
        robodamus            Plot the Robodamus predictions vs the SOS SN2...
        sdR_hdrfix           Create the files used by the pipeline to fix the...

.. _SOS_Tools_FiberQA_py:

.. admonition:: SOS Tools FiberQA
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools FiberQA [OPTIONS]
      
        Create Fiber info Summary for SOS
      
      Options:
        -s, --sosdir SOSDIR  Base SOS output directory  [required]
        -m, --mjd MJD        MJD  [required]
        -e, --exp EXPID      Exposure Name
        -n, --nocopy         Prevent copy to combined Directory
        --no-hash            Skip updating the file hash
        --help               Show this message and exit.

.. _SOS_Tools_Log_py:

.. admonition:: SOS Tools Log
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools Log [OPTIONS]
      
        Build BOSS Exposure Log
      
      Options:
        -s, --hide_summary              Hide data summary table
        -e, --hide_error                Hide SOS Error and Workings
        -r, --hart_raw                  Print raw form (instead of table form) of
                                        Hartmann Logs
        -c, --hide_hart, --hide_hartmann
                                        Hide cleaned version of Hartmann Logs as a
                                        table
        --new_ref                       Calculate new reference values in fratio and
                                        w_shift and show in place of fratio and
                                        w_shift
        -l, --long                      Long/detailed version of log
        -o, --observatory, --obs [apo|lco]
                                        Manually set observatory
        -y, --yesterday                 current mjd-1
        -m, --mjd TEXT                  MJD
        --help                          Show this message and exit.

.. _SOS_Tools_boss_arcs_to_traces_py:

.. admonition:: SOS Tools boss_arcs_to_traces
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools boss_arcs_to_traces [OPTIONS]
      
        Routine to transfer trace locations from an initial arc/flat pair given
        subsequent arc frames only
      
      Options:
        --mjd INTEGER      MJD to process  [required]
        --outdir TEXT      output directory
        --obs [lco|apo]    observatory  [default: lco]
        --vers TEXT        BOSS_SPECTRO_REDUX version  [default: master]
        --threads INTEGER  number of threads  [default: 8]
        --cams TEXT        Supply the camera for operation with SOS files
        --fitsname TEXT    Supply the FitsName for SOS error reporting
        --sosdir TEXT      Base SOS output directory
        --clobber          clobber?
        --no-hash          Skip updating the file hash
        --help             Show this message and exit.

.. _SOS_Tools_flag_manual_cal_py:

.. admonition:: SOS Tools flag_manual_cal
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools flag_manual_cal [OPTIONS]
      
        Build spManCal.par file to flag manual alternative calibration frames for
        spPlan
      
      Options:
        -o, --observatory, --obs [apo|lco]
                                        Observatory  [required]
        -m, --mjd INTEGER               MJD  [required]
        -f, --field TEXT                FieldID  [required]
        -e, --expid INTEGER             Exposure ID to manually set the calibration
                                        frame exposure ID
        -t, --type [arc|flat]           Calibration Type  [required]
        --nogit                         Skip automatic git add
        --help                          Show this message and exit.

.. _SOS_Tools_hash_py:

.. admonition:: SOS Tools hash
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools hash [OPTIONS]
      
        Create or check the SOS file hash
      
      Options:
        --mjd MJD   MJD to process  [required]
        -t, --redo  use sosredo directory (same as SOS -t or --redoMode option)
        -d, --test  use sosredo/dev directory (same as SOS -d or --test option)
        -u, --utah  use utah test directory (same as SOS --utah option)
        --create    Create the Hash file
        --check     Check the SOS Hash
        --transfer  Check the SOS Hash of a Utah transfer
        --lco       Build/check for lco at Utah
        --dummy     Create a dummy file to prevent an empty hash file
        --help      Show this message and exit.

.. _SOS_Tools_htmlIndex_py:

.. admonition:: SOS Tools htmlIndex
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools htmlIndex [OPTIONS]
      
        Build sos/combined/index.html
      
      Options:
        -s, --sosdir SOSDIR  Base SOS output directory  [required]
        -f, --force          Force Update of Index page
        --help               Show this message and exit.

.. _SOS_Tools_loadsn2_py:

.. admonition:: SOS Tools loadsn2
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools loadsn2 [OPTIONS]
      
        Load SOS SN2 values into OpsDB
      
      Options:
        --fits FITSFILE        The fits file is the science frame output from sos-
                               reduce  [required]
        --confSum CONFSUMMARY  confSummary-file  [required]
        -v, --verbose          verbose
        -u, --update           update (An error will occur if the exposure has
                               already been processed, unless set)
        --sdssv_sn2            Load sdssv_sn2
        --help                 Show this message and exit.

.. _SOS_Tools_log2html_py:

.. admonition:: SOS Tools log2html
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools log2html [OPTIONS]
      
        Create the HTML Logging Page for SOS
      
      Options:
        --mjd MJD                MJD of reduction  [required]
        --sosdir SOSDIR          Path to SOS save directory (ie. folder that
                                 contains logfile-?????.fits)  [required]
        -l, --logfile LOGFILE    Name of logfile (default: logfile-{mjd}.fits)
        -f, --htmlfile HTMLFILE  Name of output htmlfile (default:
                                 logfile-{mjd}.html)
        -o, --obs OBS            Observatory of observations (default: apo)
        -c, --copydir SAVEDIR    Where to save the htmls
        --fps                    build for FPS reductions
        --sdssv_sn2              Include SDSSV SN2 V2
        --sn2_15                 Include Mag 15 SN2
        --bright                 Include Mag 15 SN2 for all exposures
        --help                   Show this message and exit.

.. _SOS_Tools_parse_runtime_py:

.. admonition:: SOS Tools parse_runtime
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools parse_runtime [OPTIONS] LOGFILE
      
        Process log file to calculate elapsed times for SOS
      
      Options:
        -a, --all    Combine all daily logs of this format
        -s, --stamp  Add Date Stamp to output file
        --help       Show this message and exit.

.. _SOS_Tools_plot_py:

.. admonition:: SOS Tools plot
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools plot [OPTIONS]
      
        Plot the Science frame for SOS
      
      Options:
        --mjd MJD          MJD of reduction  [required]
        --expid EXPID      Exposure ID to plot  [required]
        --observatory OBS  Observatory (default: $OBSERVATORY)
        --ccd TEXT         CCDs to plot; defaults to both CCDs
        --redo             If set use sosredo rather then sos reductions
        --mask_end         Mask end of the spectra during plotting
        --ToOs             Plot only ToO fibers
        --assigned         Plot only fibers assigned to targets (includes ToOs)
        --science          Plot all fibers with science targets (includes ToOs and
                           assigned)
        --pdf              Plot all fibers into a single multi-panel PDF instead of
                           individual PNGs
        --help             Show this message and exit.

.. _SOS_Tools_robodamus_py:

.. admonition:: SOS Tools robodamus
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools robodamus [OPTIONS]
      
        Plot the Robodamus predictions vs the SOS SN2 Measurements
      
      Options:
        -m, --mjd MJD  SJD of reduction
        --help         Show this message and exit.

.. _SOS_Tools_sdR_hdrfix_py:

.. admonition:: SOS Tools sdR_hdrfix
   :collapsible: closed

   .. code-block:: text

      Usage: SOS Tools sdR_hdrfix [OPTIONS] EXPID
      
        Create the files used by the pipeline to fix the header meta data of the
        BOSS exposures
      
      Options:
        -v, --value TEXT                updated header keyword value (required if
                                        key is set)
        -k, --key TEXT                  header keyword to update (required if value
                                        is set)
        --designid INTEGER              DesignID
        --confid INTEGER                ConfigureID
        --fieldid INTEGER               FieldID
        --cartid [FPS-S|FPS-N]          Cartridge Mounted
        --tai-beg FLOAT                 Starting time (tai) of exposure
        --exptime FLOAT                 Exposure length (s)
        --flavor [bias|dark|flat|arc|science|smear]
                                        Type/Flavor of exposure
        --quality [excellent|test|bad]  Set Quality flat of exposures
        --hartmann [out|right|left|closed]
                                        Hartmann Door Status
        --flat                          short cut to set FF = 1 1 1 1 & FFS =  1 1 1
                                        1 1 1 1 1
        --arc                           short cut to set all relevant arc lamps to 1
                                        1 1 1
        --HEAR [0|1]...                 HeAr arc Lamp
        --HGCD [0|1]...                 HeCd arc Lamp
        --NE [0|1]...                   Ne arc lamp
        --FFS [0|1]...                  Flat Field Screen
        --FF [0|1]...                   Flat Field Lamp
        -t, --test                      Flag as test
        -b, --bad                       Flag as bad
        --nogit                         Skip automatic git add
        -u, --no-update                 Skip updating SOS logs
        --cameras [b1|b2|r1|r2|??]      Cameras for hdr update
        --clobber                       Clobber sdHdrFix file
        --obs [apo|lco]                 Observatory  [required]
        -m, --mjd TEXT                  MJD of file (default: latest)
        --help                          Show this message and exit.

.. _boss_flatlib:

boss_flatlib
^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _boss_flatlib_py:

.. code-block:: text

   Usage: boss_flatlib [OPTIONS] COMMAND [ARGS]...
   
     Build and analyze a library of flats to check for Fiber throughput Issues
   
   Options:
     --fullhelp  Show full help including all subcommands
     --help      Show this message and exit.
   
   Commands:
     reduce      Reduce/link the spFlats
     build       Build the flat library fits file
     plot        Plot Raw and Reduced Flat
     analyze     Run Full analysis on Flat library
     lowfiber    Check for Low fibers
     csv         Export CSV only
     timeSeries  Plot Throughout Time Series only
     end2end     Run full pipeline and plot time series (FPS only)

.. _boss_flatlib_reduce_py:

.. admonition:: boss_flatlib reduce
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib reduce [OPTIONS]
      
        Reduce/link the spFlats
      
      Options:
        --lco / --no-lco                Run for LCO data
        --run2d RUN2D                   Override $RUN2D
        -d, --dir FLATLIB_DIR           Flat Library Directory
        -c, --link_traceflat            Link all spTraceFlat files
        -a, --link_all                  Link all spFlat files regardless of
                                        spPlanTrace file
        --run                           Just link (if set), but do not run new
                                        spFlat files
        --submit / --no-submit          Submit the job to the queue
        --nodes NNODES                  Number of nodes to use
        --queue_config_file, --qcf QUEUE_CONFIG_PATH
                                        Queue Config File Path
        --queue_config, --qc QUEUE_CONFIG_NAME
                                        Queue Config name
        --deep / --no-deep              Check Pre-existing plans for completion
        --link / --no-link              Link Pre-existing spFlat Files
        -m, --mjd MJD                   MJDs to Run
        --fps                           Catch up FPS
        --plates                        Catch up Plates
        --legacy                        Catch up Legacy
        --help                          Show this message and exit.

.. _boss_flatlib_build_py:

.. admonition:: boss_flatlib build
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib build [OPTIONS]
      
        Build the flat library fits file
      
      Options:
        --lco / --no-lco       Run for LCO data
        --run2d RUN2D          Override $RUN2D
        -d, --dir FLATLIB_DIR  Flat Library Directory
        --help                 Show this message and exit.

.. _boss_flatlib_plot_py:

.. admonition:: boss_flatlib plot
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib plot [OPTIONS]
      
        Plot Raw and Reduced Flat
      
      Options:
        --lco / --no-lco       Run for LCO data
        --run2d RUN2D          Override $RUN2D
        -d, --dir FLATLIB_DIR  Flat Library Directory
        -s, --save SAVEDIR     Save Directory
        -m, --mjd MJD          List of mjds to plot
        -f, --flats FLATLIST   List of reduced flats to plot
        --help                 Show this message and exit.

.. _boss_flatlib_analyze_py:

.. admonition:: boss_flatlib analyze
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib analyze [OPTIONS]
      
        Run Full analysis on Flat library
      
      Options:
        --lco / --no-lco       Run for LCO data
        --run2d RUN2D          Override $RUN2D
        -d, --dir FLATLIB_DIR  Flat Library Directory
        -m, --mjd MJD          List of mjds to plot alone
        --plot                 Plot Flat
        --help                 Show this message and exit.

.. _boss_flatlib_lowfiber_py:

.. admonition:: boss_flatlib lowfiber
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib lowfiber [OPTIONS]
      
        Check for Low fibers
      
      Options:
        --lco / --no-lco           Run for LCO data
        --run2d RUN2D              Override $RUN2D
        -d, --dir FLATLIB_DIR      Flat Library Directory
        -m, --mjd MJD              List of mjds to plot alone
        -t, --threshold THRESHOLD  Threshold to flag lowfibers
        --help                     Show this message and exit.

.. _boss_flatlib_csv_py:

.. admonition:: boss_flatlib csv
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib csv [OPTIONS]
      
        Export CSV only
      
      Options:
        --lco / --no-lco       Run for LCO data
        --run2d RUN2D          Override $RUN2D
        -d, --dir FLATLIB_DIR  Flat Library Directory
        --help                 Show this message and exit.

.. _boss_flatlib_timeSeries_py:

.. admonition:: boss_flatlib timeSeries
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib timeSeries [OPTIONS]
      
        Plot Throughout Time Series only
      
      Options:
        --lco / --no-lco       Run for LCO data
        --run2d RUN2D          Override $RUN2D
        -d, --dir FLATLIB_DIR  Flat Library Directory
        -t, --TraceIDs         Label with Trace FiberIDs rather then slit FiberIDs
        --mjdstart MJDSTART    MJD to start reduction
        -m, --mjd MJD          List of mjds to plot alone
        --help                 Show this message and exit.

.. _boss_flatlib_end2end_py:

.. admonition:: boss_flatlib end2end
   :collapsible: closed

   .. code-block:: text

      Usage: boss_flatlib end2end [OPTIONS]
      
        Run full pipeline and plot time series (FPS only)
      
      Options:
        --lco / --no-lco                Run for LCO data
        --run2d RUN2D                   Override $RUN2D
        -d, --dir FLATLIB_DIR           Flat Library Directory
        -c, --link_traceflat            Link all spTraceFlat files
        -a, --link_all                  Link all spFlat files regardless of
                                        spPlanTrace file
        --run                           Just link (if set), but do not run new
                                        spFlat files
        --submit / --no-submit          Submit the job to the queue
        --nodes NNODES                  Number of nodes to use
        --queue_config_file, --qcf QUEUE_CONFIG_PATH
                                        Queue Config File Path
        --queue_config, --qc QUEUE_CONFIG_NAME
                                        Queue Config name
        --deep / --no-deep              Check Pre-existing plans for completion
        --link / --no-link              Link Pre-existing spFlat Files
        -t, --TraceIDs                  Label with Trace FiberIDs rather then slit
                                        FiberIDs
        --mjdstart MJDSTART             MJD to start reduction
        -m, --mjd MJD                   List of mjds to plot alone
        --help                          Show this message and exit.

.. _idlspec2d_version:

idlspec2d_version
^^^^^^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _idlspec2d_version_py:

.. code-block:: text

   Usage: idlspec2d_version [OPTIONS]
   
     Prints the IDLspec2D BOSS_DRP version
   
   Options:
     --fullhelp  Show full help including all subcommands
     --help      Show this message and exit.

.. _sdR_hdrfix:

sdR_hdrfix
^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _sdR_hdrfix_py:

.. code-block:: text

   Usage: sdR_hdrfix [OPTIONS] EXPID
   
     Create the files used by the pipeline to fix the header meta data of the
     BOSS exposures
   
   Options:
     -v, --value TEXT                updated header keyword value (required if
                                     key is set)
     -k, --key TEXT                  header keyword to update (required if value
                                     is set)
     --designid INTEGER              DesignID
     --confid INTEGER                ConfigureID
     --fieldid INTEGER               FieldID
     --cartid [FPS-S|FPS-N]          Cartridge Mounted
     --tai-beg FLOAT                 Starting time (tai) of exposure
     --exptime FLOAT                 Exposure length (s)
     --flavor [bias|dark|flat|arc|science|smear]
                                     Type/Flavor of exposure
     --quality [excellent|test|bad]  Set Quality flat of exposures
     --hartmann [out|right|left|closed]
                                     Hartmann Door Status
     --flat                          short cut to set FF = 1 1 1 1 & FFS =  1 1 1
                                     1 1 1 1 1
     --arc                           short cut to set all relevant arc lamps to 1
                                     1 1 1
     --HEAR [0|1]...                 HeAr arc Lamp
     --HGCD [0|1]...                 HeCd arc Lamp
     --NE [0|1]...                   Ne arc lamp
     --FFS [0|1]...                  Flat Field Screen
     --FF [0|1]...                   Flat Field Lamp
     -t, --test                      Flag as test
     -b, --bad                       Flag as bad
     --nogit                         Skip automatic git add
     -u, --no-update                 Skip updating SOS logs
     --cameras [b1|b2|r1|r2|??]      Cameras for hdr update
     --clobber                       Clobber sdHdrFix file
     --obs [apo|lco]                 Observatory  [required]
     -m, --mjd TEXT                  MJD of file (default: latest)
     --help                          Show this message and exit.

.. _BOSS_log:

BOSS_log
^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

.. _BOSS_log_py:

.. code-block:: text

   Usage: BOSS_log [OPTIONS]
   
     Build BOSS Exposure Log
   
   Options:
     -s, --hide_summary              Hide data summary table
     -e, --hide_error                Hide SOS Error and Workings
     -r, --hart_raw                  Print raw form (instead of table form) of
                                     Hartmann Logs
     -c, --hide_hart, --hide_hartmann
                                     Hide cleaned version of Hartmann Logs as a
                                     table
     --new_ref                       Calculate new reference values in fratio and
                                     w_shift and show in place of fratio and
                                     w_shift
     -l, --long                      Long/detailed version of log
     -o, --observatory, --obs [apo|lco]
                                     Manually set observatory
     -y, --yesterday                 current mjd-1
     -m, --mjd TEXT                  MJD
     --help                          Show this message and exit.

Full Bash Command Usage
-----------------------

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

cronplot_QA.bash
^^^^^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

::

   Usage: cronplot_QA.bash module [options]
   
   Description:
       Load the correct module and execute the QA plotting script.
   
   Options:
       -l          Use LCO observations (default is APO).
       -c          Include the --clobber_lists option.
       -n          Disable linking (default is linking enabled).
       -e          Include the --epoch option.
       -w          Generate HTML output.
       -u NAME     HTML output name (used with -w).
       -h          Display this help message and exit.
   
   Example:
       cronplot_QA.bash myModule -l -c -n -e -w -u test.html

cronrun.bash
^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

::

   usage: cronrun.bash module "command"

IDL Command Usage
-----------------

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

spreduce2d.pro
^^^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

::

   ;+
   ; NAME:
   ;   spreduce2d
   ;
   ; PURPOSE:
   ;   Calling script for SPREDUCE that reduces a night of data according
   ;   to a plan file.
   ;
   ; CALLING SEQUENCE:
   ;   spreduce2d, [ planfile, docams=, /do_telluric, /xdisplay, $
   ;    /writeflatmodel, /writearcmodel, /bbspec ]
   ;
   ; INPUTS:
   ;
   ; OPTIONAL INPUTS:
   ;   planfile   - Name(s) of output plan file; default to reducing all
   ;                plan files matching 'spPlan2d*.par'
   ;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2']
   ;   do_telluric- Passed to EXTRACT_OBJECT
   ;   xdisplay   - Send plots to X display rather than to plot file
   ;   writeflatmodel - passed to SPCALIB via SPREDUCE to trigger writing
   ;                    out of flat model info to file.
   ;   writearcmodel  - passed to SPCALIB via SPREDUCE to trigger writing
   ;                    out of arc model info to file.
   ;   bbspec         - use bbspec extraction code
   ;   noreject       - Override Bad calibration rejection (use with caution)
   ;
   ; Optional Keywords:
   ;   MWM_fluxer  - Utilize MWM optional settings (ie gaia reddening and different S/N cuts)
   ;
   ;
   ; OUTPUT:
   ;
   ; COMMENTS:
   ;   The following environment variables must be set:
   ;      BOSS_SPECTRO_DATA
   ;      SDSSCORE
   ;      SPECFLAT_DIR
   ;   Look for raw FITS data files in BOSS_SPECTRO_DATA/MJD.
   ;   Look for obsSummary files in SDSSCORE/MJD.
   ;   Look for spectroscopic flat files in SPECFLAT_DIR.
   ;
   ; EXAMPLES:
   ;
   ; BUGS:
   ;   This routine spawns the Unix command 'mkdir'.
   ;
   ; PROCEDURES CALLED:
   ;   cpbackup
   ;   idlspec2d_version()
   ;   idlutils_version()
   ;   splog
   ;   spreduce
   ;   yanny_free
   ;   yanny_par()
   ;   yanny_read
   ;
   ; INTERNAL SUPPORT ROUTINES:
   ;
   ; REVISION HISTORY:
   ;   02-Nov-1999  Written by David Schlegel, Princeton.
   ;      Apr-2010  Added "write[flat,arc]model" pass-through (A. Bolton, Utah)
   ;   15-Aug-2011  Added pass-through for spatial split of sky model (A. Bolton, Utah)
   ;   15-Nov-2018: Modified for use only one spectrograph for the BHM (HJIM)
   ;-

rm_combine_script.pro
^^^^^^^^^^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

::

   ;+
   ; NAME:
   ;   rm_combine_script
   ;
   ; PURPOSE:
   ;   Script to process epochs with the xyfit custom flux calibration
   ;
   ; CALLING SEQUENCE:
   ;
   ; INPUTS:
   ;   planfile   - Name(s) of output plan file
   ;
   ; OPTIONAL INPUTS:
   ;   run2d      - Name of the run2d
   ;   finaldir   - Additional subdirectory for output
   ;   xyfit      - Compute 2d flux corrections in the xy focal plane
   ;   bscore     - Fraction of best exposure score to use as a threshold for discarding exposures
   ;   minsn2     - Minimum S/N^2 to include science frame in coadd; default
   ;                to 0 to only include those with S/N > 0.
   ;                Note that all exposures with a score less than 0.2 times
   ;                the score of the best exposure are discarded; for those
   ;                purposes, the score used is the worst of all 4 cameras.
   ;
   ;
   ; Optional Keywords:
   ;   MWM_fluxer    - Utilize MWM optional settings (ie gaia reddening and different S/N cuts)
   ;   nofcorr       - Skip the step to generate and use the spFluxcorr* files
   ;   nodist        - Skip the step to generate and use the spFluxdistort* files
   ;   radec_coadd   - Coadd using ra-dec matching rather then catalogID matching
   ;   no_reject     - Turns off rejection in the coadding
   ;   onestep_coadd - Legacy algorithm for coadd. Coadding blue+red and all exposures
   ;                    at the the same time.
   ;   epoch         - Epoch Coadd flag for input and outputs
   ;   legacy        - Flag for Pre-SDSSV 2 Spectrograph data at APO
   ;   plates        - Flat for SDSSV 1 Spectrograph plate data at APO
   ;   loaddesi      - Load the DESI (JG) models for fluxing
   ;   skipfluxing   - Skip the step to generate spFluxcalib* files
   ;   skipfcorr     - Skip creation of flux-correction vectors and use prexisting spFluxcorr* files
   ;
   ; OUTPUT:
   ;
   ; COMMENTS:
   ; EXAMPLES:
   ;
   ; BUGS:
   ;   This routine spawns the Unix command 'mkdir'.
   ;
   ; PROCEDURES CALLED:
   ;   get_field_dir
   ;   djs_filepath
   ;   rm_spcombine_v5
   ;
   ;

spreduce1d_empca.pro
^^^^^^^^^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

::

   ;+
   ; NOTE: same as spreduce1d, but uses different QSO PCA templates
   ; NAME:
   ;   spreduce1d
   ;
   ; PURPOSE:
   ;   1-D reduction of spectra from 1 plate
   ;
   ; CALLING SEQUENCE:
   ;   spreduce1d, [ platefile, fiberid=, run1d=, /doplot, /debug, chop_data= ]
   ;
   ; INPUTS:
   ;
   ; OPTIONAL INPUTS:
   ;   platefile  - Plate file(s) from spectro-2D; default to all files
   ;                matching 'spPlate*.fits'
   ;   fiberid    - If specified, then only reduce these fiber numbers;
   ;                this must be a vector with unique values between 1 and
   ;                the number of fibers in the plate file
   ;   run1d      - Optional override value for the environment variable $RUN1D
   ;   doplot     - If set, then generate plots.  Send plots to a PostScript
   ;                file spDiagDebug1d-$PLATE-$MJD.ps unless /DEBUG is set.
   ;   debug      - If set, then send plots to the X display and wait for
   ;                a keystroke after each plot; setting /DEBUG forces /DOPLOT.
   ;   chop_data  - If set, then trim wavelength range to the specified range
   ;                in vacuum Ang (if a 2-element array), or to a default
   ;                trim range of [3600,10400] Ang.
   ;
   ; OUTPUTS:
   ;
   ; OPTIONAL OUTPUTS:
   ;
   ; COMMENTS:
   ;   Input files are read from the current directory.
   ;   Output files are written to the subdirectory $RUN1D.
   ;
   ;   Names of output files are derived from PLATEFILE.
   ;   For example, if PLATEFILE='spPlate-0306-51690.fits', then
   ;     ZALLFILE = 'spZall-0306-51690.fits'
   ;     ZBESTFILE = 'spZbest-0306-51690.fits'
   ;     ZLINEFILE = 'spZline-0306-51690.fits'
   ;
   ; EXAMPLES:
   ;
   ; BUGS:
   ;
   ; DATA FILES:
   ;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
   ;
   ; PROCEDURES CALLED:
   ;   cpbackup
   ;   dfpsclose
   ;   dfpsplot
   ;   djs_filepath()
   ;   elodie_best()
   ;   fileandpath()
   ;   filter_thru()
   ;   mrdfits()
   ;   mwrfits
   ;   qaplot_fcalibvec
   ;   splog
   ;   skymask()
   ;   speclinefit
   ;   star_dvelocity()
   ;   struct_addtags()
   ;   sxaddpar
   ;   sxdelpar
   ;   sxpar()
   ;   synthspec()
   ;   vdispfit
   ;   zfind()
   ;   zrefind()
   ;
   ; REVISION HISTORY:
   ;   28-Jun-2000  Written by D. Schlegel, Princeton
   ;   2010-2011: various template-related tweaks and Z_NOQSO, A. Bolton, Utah
   ;   01-Oct-2012: Adding ZNUM_NOQSO to the Z_NOQSO section, Joel Brownstein, Utah

spspec_target_merge.pro
^^^^^^^^^^^^^^^^^^^^^^^

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none

::

   ;+
   ; NAME:
   ;   spspec_target_merge
   ;
   ; PURPOSE:
   ;   To create spSpec and spFullsky target level coadds (independent of field-mjd) 
   ;
   ; CALLING SEQUENCE:
   ;
   ; INPUTS:
   ;   customplan - The spPlanCustom file for the coadd
   ;
   ; OPTIONAL KEYWORDS:
   ;   topdir - the daily coadd base directory
   ;
   ; OUTPUTS:
   ;
   ; OPTIONAL OUTPUTS:
   ;
   ; COMMENTS:
   ;
   ; EXAMPLES:
   ;
   ; BUGS:
   ;
   ; PROCEDURES CALLED:
   ;
   ; REVISION HISTORY:
   ;
   ;-

.. highlight:: defaults

.. End of document
