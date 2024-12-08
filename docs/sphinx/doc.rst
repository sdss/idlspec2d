:tocdepth: 2

.. highlight:: none

Full Command Documention
========================
Documented below are the primary commands used to run the BOSS Data Reduction Pipeline. However, there are numerous other routines included in this package, which are called by these commands and have their own internal documentation.

Full Bash and Python Command Usage
----------------------------------

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none


.. _BOSS_log:

BOSS_log
^^^^^^^^
::
 
    usage: BOSS_log [-h] [-m MJD] [-y] [-o {apo,lco}] [-l] [--new_ref] [-c] [-r]
                    [-e] [-s]
    
    Build BOSS Exposure Log
    
    optional arguments:
      -h, --help            show this help message and exit
      -m MJD, --mjd MJD     MJD
      -y, --yesterday       current mjd-1
      -o {apo,lco}, --observatory {apo,lco}, --obs {apo,lco}
                            Manually set observatory
      -l, --long            Long/detailed version of log
      --new_ref             Calculate new reference values in fratio and w_shift
                            and show in place of fratio and w_shift (edit to code
                            to save new value is required)
      -c, --hide_hart, --hide_hartmann
                            Hide cleaned version of Hartmann Logs as a table
      -r, --hart_raw        Print raw form (instead of table form) of Hartmann
                            Logs
      -e, --hide_error      Hide SOS Error and Workings
      -s, --hide_summary    Hide data summary table

.. _SOS:

SOS
^^^
::
 
    usage: SOS [-h] (-r | -b | -j) (-c | -t | -d) [-e EXP] [-m [MJD [MJD ...]]]
               [--nodb] [--no_gz] [--no_reject] [--clobber_fibermap]
               [--no_sdssv_sn2] [--no_sn2_15] [-n] [-o] [-v]
    
    SOS process for reducing BOSS data on the Moutain
    
    optional arguments:
      -h, --help            show this help message and exit
      -r, --red             Red Camera Process
      -b, --blue            Blue Camera Process
      -j, --joint           Both Camera Processes
      -c, --catchup         Run Catchup on the night or (MJD)
      -t, --redoMode        Save outputs of MJD or exposure to sosredo
      -d, --test            Save outputs and logs to sosredo/dev
      -e EXP, --exp EXP     exposure id (or range of exp id 500-510) (with or
                            without leading zeros)
      -m [MJD [MJD ...]], --mjd [MJD [MJD ...]]
                            MJD
      --nodb                skip opsdb load
      --no_gz               Overrides the requirement for '.gz' compressed files
                            (experimental)
      --no_reject           Overrides the Calibration rejection (use with caution)
      --clobber_fibermap, -f
                            Clobbers the existing spfibermap files
      --no_sdssv_sn2        Skip reporting a second set of SN2 values with updated
                            fit parameters
      --no_sn2_15           Skip reporting a set of SN2 values with a fiducial mag
                            of 15
      -n, --no_arc2trace    Skip Utilizing arc2trace refinements
      -o, --forcea2t        Force arc2trace for all fields (even if flat exists
                            for field)
      -v, --verbose         prints the only (or red if joint) active SOS process
                            to terminal

.. _boss_arcs_to_traces:

boss_arcs_to_traces
^^^^^^^^^^^^^^^^^^^
::
 
    Traceback (most recent call last):
      File "/Users/smorrison/Documents/Scripts/SDSSV_idlspec2D/TraceTweakClean/python/boss_drp/../../bin/boss_arcs_to_traces", line 7, in <module>
        from pyvista import boss
    ModuleNotFoundError: No module named 'pyvista'

.. _build_combined_html:

build_combined_html
^^^^^^^^^^^^^^^^^^^
::
 
    usage: build_combined_html [-h] [--force] sosdir
    
    build SOS combine index page
    
    positional arguments:
      sosdir      Base SOS output directory
    
    optional arguments:
      -h, --help  show this help message and exit
      --force     Force update

.. _cronrun.bash:

cronrun.bash
^^^^^^^^^^^^
::
 
    usage: cronrun.bash module 'script'
     

.. _daily_log:

daily_log
^^^^^^^^^
::
 
    usage: daily_log [-h] [--obs OBS [OBS ...]] [--mjd [MJD [MJD ...]]]
                     [--mjdstart MJDSTART] [--mjdend MJDEND] [--epoch]
                     [--custom CUSTOM] [--topdir TOPDIR] [--run1d RUN1D]
                     [--run2d RUN2D] [--email] [--fast]
    
    Build/load BOSS Pipeline Status Pages
    
    optional arguments:
      -h, --help            show this help message and exit
      --obs OBS [OBS ...]   Observatory for status update
      --mjd [MJD [MJD ...]]
                            Update these MJDs
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
      --epoch               Run for epoch Coadds
      --custom CUSTOM       Name of custom Coadd
      --topdir TOPDIR       Optional override value for the environment variable
                            $BOSS_SPECTRO_REDUX
      --run1d RUN1D         Optional override value for the enviro variable $RUN1D
      --run2d RUN2D         Optional override value for the enviro variable $RUN2D
      --email               Send each mjd status as email
      --fast                Skip updating index until end

.. _fieldlist:

fieldlist
^^^^^^^^^
::
 
    usage: fieldlist [-h] [--create] [--topdir TOPDIR]
                     [--run1d [RUN1D [RUN1D ...]]] [--run2d [RUN2D [RUN2D ...]]]
                     [--outdir OUTDIR] [--skipcart [SKIPCART [SKIPCART ...]]]
                     [--epoch] [--basehtml BASEHTML] [--logfile LOGFILE] [--debug]
                     [--noplot]
    
    Build/load BOSS Fieldlist
    
    optional arguments:
      -h, --help            show this help message and exit
      --create, -c          Create Fieldlist
      --topdir TOPDIR       Optional override value for the environment variable
                            $BOSS_SPECTRO_REDUX
      --run1d [RUN1D [RUN1D ...]]
                            Optional override value for the enviro variable $RUN1D
      --run2d [RUN2D [RUN2D ...]]
                            Optional override value for the enviro variable $RUN2D
      --outdir OUTDIR       Optional output directory (defaults to topdir/$RUN2D)
      --skipcart [SKIPCART [SKIPCART ...]]
                            Option list of cartridges to skip
      --epoch               Produce FieldList for epoch coadds
      --basehtml BASEHTML   html path for figure (defaults to relative from
                            topdir)
      --logfile LOGFILE     Manually Set logfile (including path)
      --debug               Overrides the logger of the simplified error messages
                            and prints standard python errors
      --noplot              Skips updating the sky plots

.. _fieldmerge:

fieldmerge
^^^^^^^^^^
::
 
    usage: fieldmerge [-h] [--run2d RUN2D] [--indir INDIR] [--skip_line]
                      [--include_bad] [--legacy] [--skip_specprimary] [--lite]
                      [--XCSAO] [--field FIELD] [--mjd MJD] [--clobber] [--bkup]
                      [--verbose] [--logfile LOGFILE] [--epoch]
                      [--programs [PROGRAMS [PROGRAMS ...]]]
                      [--datamodel DATAMODEL] [--line_datamodel LINE_DATAMODEL]
                      [--outroot OUTROOT] [--remerge_fmjd REMERGE_FMJD]
                      [--remerge_mjd REMERGE_MJD] [--merge_only] [--allsky]
                      [--custom CUSTOM] [--run1d RUN1D] [--limit LIMIT]
                      [--ndays MJDSTART]
    
    Build BOSS spAll Summary File
    
    optional arguments:
      -h, --help            show this help message and exit
      --run2d RUN2D         Optional override value for the enviro variable $RUN2D
      --indir INDIR         Optional override value for the environment variable
                            $BOSS_SPECTRO_REDUX
      --skip_line           skip the generation of spAllLine.fits
      --include_bad         include bad fields
      --legacy              Include columns used by SDSS-IV and depreciated in
                            SDSS-V
      --skip_specprimary    Skip creation of specprimary and associated columns
      --lite                Produce lite version of spAll file
      --XCSAO               Include XCSAO columns
      --field FIELD, -f FIELD
                            Run for a single Field
      --mjd MJD, -m MJD     Run for a single MJD
      --clobber             Clobber all spAll-field-mjd files
      --bkup                Backup existing spAll files
      --verbose             Log columns not saved
      --logfile LOGFILE     Manually set logfile
      --epoch               Produce spAll for epoch coadds
      --programs [PROGRAMS [PROGRAMS ...]]
                            List of programs to include
      --datamodel DATAMODEL
                            Supply a spAll datamodel file (defaults to
                            $IDLSPEC2D/datamodel/spall_dm.par
      --line_datamodel LINE_DATAMODEL
                            Supply a spline datamodel file (defaults to
                            $IDLSPEC2D/datamodel/spzline_dm.par
      --outroot OUTROOT     Path and root of filename for output (defaults to
                            spectra/full or summary)
      --remerge_fmjd REMERGE_FMJD, -r REMERGE_FMJD
                            Field-MJD to replace in spAll
      --remerge_mjd REMERGE_MJD
                            MJD to replace in spAll
      --merge_only, -o      Skip Building new spAll-Field-MJD files and just merge
                            existing
      --allsky              Build spAll for Allsky Custom Coadd
      --custom CUSTOM       Name of Custom Coadd
      --run1d RUN1D         Optional override value for the enviro variable $RUN1D
                            (only for custom allsky coadds)
      --limit LIMIT         Limit number of Field-MJD spAll files to read before
                            save
      --ndays MJDSTART      Limit update to last ndays

.. _filecheck:

filecheck
^^^^^^^^^
::
 
    usage: filecheck [-h] cmd file
    
            Check File (uncompressed or gz) favor/instrument/quality
            
            science:
              return "true" if the fits file is a science frame.  This
              is determined by flavor=science in the header.  If flavor
              is not in the header, "false" is returned.
        
            test:
              return "true" if the fits file is a test frame.  This is
              determined by quality=test in the header.  If quality
              is not in the header, "false" is returned
        
            excellent:
              return "true" if the fits file is a excellent frame.  This is
              determined by quality=excellent in the header.  If quality
              is not in the header, "true" is returned
        
            boss:
              return "true" if the plPlugMapM file is a boss frame.
              this is determined by instrument=boss in the header.
              If instrument is not in the header, "false" is returned.
            
    
    positional arguments:
      cmd         file check command
      file        fits file
    
    optional arguments:
      -h, --help  show this help message and exit

.. _fluxcorr_prior:

fluxcorr_prior
^^^^^^^^^^^^^^
::
 
    usage: fluxcorr_prior [-h] [--xythrucorr] planfile
    
    Try solving with a prior that fluxcorr = 1
    
    positional arguments:
      planfile      name of the plan file
    
    optional arguments:
      -h, --help    show this help message and exit
      --xythrucorr  Apply XY throughput Correction

.. _idlspec2d_version:

idlspec2d_version
^^^^^^^^^^^^^^^^^
::
 
    usage: idlspec2d_version [-h]
    
    Prints the IDLspec2D BOSS_DRP version
    
    optional arguments:
      -h, --help  show this help message and exit

.. _loadSN2Value:

loadSN2Value
^^^^^^^^^^^^
::
 
    usage: loadSN2Value [-h] [-v] [-u] [--sdssv_sn2] fits confSum
    
    Load SOS SN2 values into OpsDB
    
    positional arguments:
      fits           The fits file is the science frame output from sos-reduce
      confSum        confSummary-file
    
    optional arguments:
      -h, --help     show this help message and exit
      -v, --verbose  verbose
      -u, --update   update (An error will occur if the exposure has already been
                     processed, unless set)
      --sdssv_sn2    Load sdssv_sn2

.. _manage_coadd_Schema:

manage_coadd_Schema
^^^^^^^^^^^^^^^^^^^
::
 
    usage: manage_coadd_Schema [-h] [--coaddfile COADDFILE] [--topdir TOPDIR]
                               [--run2d RUN2D] [--name NAME] [--DR] [--rerun1d]
                               [--active] [--carton [CARTON [CARTON ...]]]
                               [--SDSSIDS [SDSSIDS [SDSSIDS ...]]]
                               [--program [PROGRAM [PROGRAM ...]]]
                               [--legacy [LEGACY [LEGACY ...]]] [--use_catid]
                               [--use_firstcarton] [--cadence CADENCE] [--show]
                               [--mjd [MJD [MJD ...]]]
    
    Manage Custom Coadds
    
    optional arguments:
      -h, --help            show this help message and exit
      --coaddfile COADDFILE, -f COADDFILE
                            File to store Coadding Schema (Default:
                            {topdir}/{run2d}/fields/SDSSV_BHM_COADDS.par)
      --topdir TOPDIR       Override value for the environment variable
                            $BOSS_SPECTRO_REDUX.
      --run2d RUN2D         Override value for the environment variable $RUN2D
      --name NAME           Name of Custom Coadd
      --DR                  DR/IPL Coadding
      --rerun1d, -r         Provides flag for coadd to be rerun though 1D analysis
      --active, -a          Activate (or deactivate) a Coadding Schema
      --carton [CARTON [CARTON ...]], -c [CARTON [CARTON ...]]
                            list of cartons
      --SDSSIDS [SDSSIDS [SDSSIDS ...]], -i [SDSSIDS [SDSSIDS ...]]
                            list of SDSS_IDS (or CatalogIDs if use_catid is set)
      --program [PROGRAM [PROGRAM ...]], -p [PROGRAM [PROGRAM ...]]
                            list of programs
      --legacy [LEGACY [LEGACY ...]], -l [LEGACY [LEGACY ...]]
                            list of Legacy Tags to include
      --use_catid, -u       Use CatalogIDs rather then SDSS_IDs
      --use_firstcarton     Use Firstcarton only for carton match (dont look at
                            db)
      --cadence CADENCE, -t CADENCE
                            Number of days between coadd epochs
      --show, -s            Show Configurations
      --mjd [MJD [MJD ...]]
                            Use data from these MJDs.

.. _parse_runtime:

parse_runtime
^^^^^^^^^^^^^
::
 
    usage: parse_runtime [-h] [-a] [-s] file_path [file_path ...]
    
    Process log file to calculate elapsed times for SOS.
    
    positional arguments:
      file_path    Path to the log file
    
    optional arguments:
      -h, --help   show this help message and exit
      -a, --all    Combine all daily logs of this format
      -s, --stamp  Add Date stamp to output file

.. _plot_qa:

plot_qa
^^^^^^^
::
 
    usage: plot_qa [-h] [-r [RUN2D [RUN2D ...]]] [-t [TEST [TEST ...]]]
                   [--test_path TEST_PATH] [--mjds_low [MJDS_LOW [MJDS_LOW ...]]]
                   [--mjds_high [MJDS_HIGH [MJDS_HIGH ...]]] [--clobber_lists]
                   [--lco] [--publish] [-e]
    
    Plot QA
    
    optional arguments:
      -h, --help            show this help message and exit
      -r [RUN2D [RUN2D ...]], --run2d [RUN2D [RUN2D ...]]
                            List of run2ds
      -t [TEST [TEST ...]], --test [TEST [TEST ...]]
                            List of True/False test run2d (corresponding to run2d)
      --test_path TEST_PATH
                            test Run2d path modification
      --mjds_low [MJDS_LOW [MJDS_LOW ...]]
                            List of mjd lower limits - use None for no limit
                            (corresponding to run2d)
      --mjds_high [MJDS_HIGH [MJDS_HIGH ...]]
                            List of mjd upper limits - use None for no limit
                            (corresponding to run2d)
      --clobber_lists       Clobber list of fieldIDs
      --lco                 Flag for LCO vs APO
      --publish             create publication version of plot
      -e, --epoch           produce plots for epoch coadds

.. _read_sos:

read_sos
^^^^^^^^
::
 
    usage: read_sos [-h] [--exp EXP] [--nocopy] directory mjd
    
    Create Fiber info Summary for SOS
    
    positional arguments:
      directory          SOS Directory
      mjd                mjd
    
    optional arguments:
      -h, --help         show this help message and exit
      --exp EXP, -e EXP  Exposure Name
      --nocopy, -n       Prevent copy to combined Directory

.. _readfibermaps:

readfibermaps
^^^^^^^^^^^^^
::
 
    usage: readfibermaps [-h] [-p SPPLAN2D] [--topdir TOPDIR] [-c] [--fast]
                         [--datamodel DATAMODEL] [-s] [--release RELEASE]
                         [--remote] [--dr19] [--confSummary CONFSUMMARY]
                         [--ccd {b2,r2,b1,r1}] [--mjd MJD] [--log]
    
    Produces spfibermap file corresponding to a spplan2d (or single confSummary
    file for SOS)
    
    optional arguments:
      -h, --help            show this help message and exit
      -p SPPLAN2D, --spplan2d SPPLAN2D
                            spplan2d file for idlspec2d run
      --topdir TOPDIR       Alternative output directory (defaults to location of
                            spplan2d file or /data/boss/sos/{mjd} for SOS)
      -c, --clobber         overwrites previous spfibermap file
      --fast                When using --no_db, streamlines process and only gets
                            parallax from MOS target files
      --datamodel DATAMODEL
                            Supply a datamodel file (defaults to
                            $IDLSPEC2D/datamodel/spfibermap_dm.par or
                            $IDLSPEC2D/datamodel/spfibermap_sos_dm.par for SOS)
      -s, --SOS             produces spfibermap for SOS
      --release RELEASE     sdss_access data release (defaults to sdsswork),
                            required if you do not have proprietary access,
                            otherwise see https://sdss-
                            access.readthedocs.io/en/latest/auth.html#auth
      --remote              allow for remote access to data using sdss-access
      --dr19                Limit targeting flags to DR19 cartons
    
    SOS:
      Options of use with SOS only
    
      --confSummary CONFSUMMARY
                            confSummary file for SOS (required for with --SOS)
      --ccd {b2,r2,b1,r1}   CCD for SOS
      --mjd MJD             MJD of observation
      --log                 creates log file in topdir

.. _run_PyXCSAO:

run_PyXCSAO
^^^^^^^^^^^
::
 
    WARNING: pyxcsao is not installed
    usage: run_PyXCSAO [-h] [--run1d RUN1D] [--epoch] fitsfile
    
    Runs pyXCSAO to determine RVs
    
    positional arguments:
      fitsfile              fits file
    
    optional arguments:
      -h, --help            show this help message and exit
      --run1d RUN1D, -r RUN1D
                            run1d name
      --epoch               run for epoch Coadds

.. _sdR_hdrfix:

sdR_hdrfix
^^^^^^^^^^
::
 
    usage: sdR_hdrfix [-h] [--mjd MJD] --obs {APO,LCO} [--clobber]
                      [--cameras {b1,b2,r1,r2,??}] [--bad] [--test]
                      [--FF {0,1} {0,1} {0,1} {0,1}]
                      [--FFS {0,1} {0,1} {0,1} {0,1} {0,1} {0,1} {0,1} {0,1}]
                      [--NE {0,1} {0,1} {0,1} {0,1}]
                      [--HGCD {0,1} {0,1} {0,1} {0,1}]
                      [--HEAR {0,1} {0,1} {0,1} {0,1}] [--arc] [--flat]
                      [--hartmann {Out,Right,Left,Closed}]
                      [--quality {excellent,test,bad}]
                      [--flavor {bias,dark,flat,arc,science,smear}]
                      [--exptime EXPTIME] [--tai-beg TAI_BEG]
                      [--cartid {FPS-S,FPS-N}] [--fieldid FIELDID]
                      [--confid CONFIGID] [--designid DESIGNID] [--key KEY]
                      [--value VALUE]
                      expid
    
    Create the files used by the pipeline to fix the header meta data of the BOSS
    exposures
    
    positional arguments:
      expid                                                  Exposure ID
    
    optional arguments:
      -h, --help                                             show this help
                                                             message and exit
      --mjd MJD, -m MJD                                      mjd of file (default:
                                                             latest MJD)
      --clobber                                              clobber sdHdrFix file
      --cameras {b1,b2,r1,r2,??}                             Cameras for hdr
                                                             update (?? for all
                                                             cameras) [default:??]
    
    Required arguments:
      --obs {APO,LCO}                                        Set Observatory
    
    Optional Quality Update (exclusive)
        At current only use if still exposing or don't run SOS after for Science Frames
         (skip and note in Night Log (and/or email) if uncertain):
      --bad, -b                                              flag as quality=bad
      --test, -t                                             flag as quality=test
    
    Optional lamp/screen keys to Update (1:on, 0:off):
      --FF {0,1} {0,1} {0,1} {0,1}                           Flat Field Lamp
      --FFS {0,1} {0,1} {0,1} {0,1} {0,1} {0,1} {0,1} {0,1}  Flat Field Screen
      --NE {0,1} {0,1} {0,1} {0,1}                           Ne arc lamp
      --HGCD {0,1} {0,1} {0,1} {0,1}                         HeCd arc Lamp
      --HEAR {0,1} {0,1} {0,1} {0,1}                         HeAr arc Lamp
      --arc                                                  short cut to set all
                                                             relevant arc lamps to
                                                             1 1 1 1
      --flat                                                 short cut to set FF =
                                                             1 1 1 1 & FFS = 1 1 1
                                                             1 1 1 1 1
      --hartmann {Out,Right,Left,Closed}                     Hartmann Door Status
    
    Optional Common keys to Update
        At current only use if still exposing or don't run SOS after for Science Frames
         (skip and note in Night Log (and/or email) if uncertain):
      --quality {excellent,test,bad}                         Set Quality flat of
                                                             exposures
    
    Optional Specialized Keys to Update 
        At current only use if still exposing or don't run SOS after
         (skip and note in Night Log (and/or email) if uncertain):
      --flavor {bias,dark,flat,arc,science,smear}            Type/Flavor of
                                                             exposure
      --exptime EXPTIME                                      Exposure length (s)
      --tai-beg TAI_BEG                                      Starting time (tai)
                                                             of exposure
      --cartid {FPS-S,FPS-N}                                 Cartridge Mounted
      --fieldid FIELDID                                      FieldID
      --confid CONFIGID                                      ConfigureID
      --designid DESIGNID                                    DesignID
    
    Manually update a key 
        At current only use if still exposing or don't run SOS after
         (skip and note in Night Log (and/or email) if uncertain):
      --key KEY, -k KEY                                      header keyword to
                                                             update (required if
                                                             value is set)
      --value VALUE, -v VALUE                                updated header
                                                             keyword value
                                                             (required if key is
                                                             set)
    
    one or more update options are required

.. _slurm_Summary:

slurm_Summary
^^^^^^^^^^^^^
::
 
    usage: slurm_Summary [-h] [--module MODULE] [--topdir TOPDIR] [--run2d RUN2D]
                         [--run1d RUN1D] [--walltime WALLTIME] [--fast]
                         [--mem MEM] [--daily] [--epoch] [--custom CUSTOM]
                         [--full] [--monitor] [--no_submit] [--merge_only]
                         [--no_fieldlist] [--backup BACKUP] [--limit LIMIT]
                         [--n_iter N_ITER] [--log2daily] [--email_start]
                         [--skip_specprimary] [--verbose]
    
    Create daily field merge slurm job
    
    optional arguments:
      -h, --help            show this help message and exit
      --module MODULE, -m MODULE
                            module file to use (ex bhm/master[default] or
                            bhm/v6_0_9)
      --topdir TOPDIR       Boss Spectro Redux base directory
      --run2d RUN2D         Run2d
      --run1d RUN1D         Run1d
      --walltime WALLTIME, -w WALLTIME
                            Job wall time (format hh:mm:ss) default = "40:00:00"
      --fast                use fast allocation
      --mem MEM             memory in bytes
      --daily               only run if daily run has been run today
      --epoch               run for the epoch coadds
      --custom CUSTOM       Name of custom Coadd
      --full                Use a full cluster node
      --monitor             Monitor job and send email at completion with the logs
      --no_submit           Create slurm job but do not submit it
      --merge_only          Run fieldmerge in merge_only mode
      --no_fieldlist        Skip Running Fieldlist
      --backup BACKUP       Number of backups to keep, or None to not create
                            backup
      --limit LIMIT         Limit number of new field-mjds to update
      --n_iter N_ITER       number of iterations of field merge to run
      --log2daily           Save Logs to $DAILY_DIR/logs/Summary/
      --email_start         Send email at start of run
      --skip_specprimary    Skip building specprimary in fieldmerge
      --verbose             Run Fieldmerge with verbose

.. _slurm_readfibermap:

slurm_readfibermap
^^^^^^^^^^^^^^^^^^
::
 
    usage: slurm_readfibermap [-h] [--module MODULE] [--topdir TOPDIR]
                              [--run2d RUN2D] [--clobber] [--apo] [--lco] [--dr19]
                              [--mjd [MJD [MJD ...]]] [--mjdstart MJDSTART]
                              [--mjdend MJDEND] [--mem_per_cpu MEM_PER_CPU]
                              [--walltime WALLTIME] [--ppn PPN]
    
    Create daily field merge slurm job. Without access to the SDSS Slurm package,
    it prints the commands for manual execution
    
    optional arguments:
      -h, --help            show this help message and exit
      --module MODULE, -m MODULE
                            module file to use (ex bhm/master or bhm/v6_0_9)
      --topdir TOPDIR       Boss Spectro Redux base directory
      --run2d RUN2D         Run2d
      --clobber             Clobber spfibermaps
      --apo                 run apo
      --lco                 run lco
      --dr19                Limit targeting flags to DR19 cartons
    
    Select MJDs:
      --mjd [MJD [MJD ...]]
                            MJD dates to reduce; default="*"
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
    
    Slurm Options:
      --mem_per_cpu MEM_PER_CPU
                            Memory allocated per CPU
      --walltime WALLTIME   Wall time in hours
      --ppn PPN             Number of processors per node

.. _slurm_sos:

slurm_sos
^^^^^^^^^
::
 
    usage: slurm_sos [-h] [--apo] [--lco] [--mjd [MJD [MJD ...]]]
                     [--mjdstart MJDSTART] [--mjdend MJDEND] [--no_reject]
                     [--clobber_fibermap] [--no_sdssv_sn2] [-n] [-o]
                     [--mem_per_cpu MEM_PER_CPU] [--walltime WALLTIME]
                     [--nodes NODES] [--ppn PPN] [--no_submit]
    
    Create SOS slurm job. Without access to the SDSS Slurm package, it prints the
    commands for manual execution
    
    optional arguments:
      -h, --help            show this help message and exit
      --apo                 run apo
      --lco                 run lco
    
    Select MJDs:
      --mjd [MJD [MJD ...]]
                            MJD dates to reduce; default=Today
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
    
    SOS Options:
      --no_reject           Overrides the Calibration rejection (use with caution)
      --clobber_fibermap, -f
                            Clobbers the existing spfibermap files
      --no_sdssv_sn2        Skip reporting a second set of SN2 values with updated
                            fit parameters
      -n, --no_arc2trace    Skip Utilizing arc2trace refinements
      -o, --forcea2t        Force arc2trace for all fields (even if flat exists
                            for field)
    
    Slurm Options:
      --mem_per_cpu MEM_PER_CPU
                            Memory allocated per CPU
      --walltime WALLTIME   Wall time in hours
      --nodes NODES         Number of nodes to use; default=1
      --ppn PPN             Number of processors per node
      --no_submit           Skip submitting process to queue

.. _slurm_spTrace:

slurm_spTrace
^^^^^^^^^^^^^
::
 
    usage: slurm_spTrace [-h] [--module MODULE] [--topdir TOPDIR] [--run2d RUN2D]
                         [--mjd [MJD [MJD ...]]] [--mjdstart MJDSTART]
                         [--mjdend MJDEND] [--lco] [--clobber] [--debug]
                         [--skip_plan] [--nodes NODES]
    
    Create spTrace slurm jobs. Without access to the SDSS Slurm package, it prints
    the commands for manual execution.
    
    optional arguments:
      -h, --help            show this help message and exit
      --module MODULE, -m MODULE
                            module file to use (ex bhm/master or bhm/v6_0_9)
      --topdir TOPDIR       Boss Spectro Redux base directory
      --run2d RUN2D         Run2d
      --mjd [MJD [MJD ...]]
                            Use data from these MJDs.
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
      --lco                 Build Run files for LCO
      --clobber             Clobber the existing Plan files
      --debug               Run in debug mode
      --skip_plan           Skip creating plans and use currently existing plans
      --nodes NODES         Number of nodes to use to run arc2trace

.. _sos_command:

sos_command
^^^^^^^^^^^
::
 
    usage: sos_command -f name -i path -p name -l path -s path -m 00000 [-d -e]
     
       -f    Fits file name
       -i    Fits file directory path
       -p    Plugmap file name
       -l    Plugmap file directory path
       -s    SOS Directory
       -m    MJD
     
       -d    Dry run.
       -e    FPS mode
       -a    no cal mode
       -n    no OpsDB upload
       -v    calculate SN2_v2 (SDSS-V)
     
    All parameters except -d, -a, and -e are required, FPS mode is set by default.
    Normally sos_command will be called by sos_runnerd.

.. _spSpec_reformat:

spSpec_reformat
^^^^^^^^^^^^^^^
::
 
    usage: spSpec_reformat [-h] --field FIELD --mjd MJD [--topdir TOPDIR]
                           [--run2d RUN2D] [--run1d RUN1D] [--custom CUSTOM]
                           [--plot] [--epoch] [--lsdr10] [--allsky]
    
    Build Spec Files
    
    optional arguments:
      -h, --help            show this help message and exit
      --field FIELD, -f FIELD
                            Run for a single Field
      --mjd MJD, -m MJD     Run for a single MJD
      --topdir TOPDIR       Optional override value for the environment variable
                            $BOSS_SPECTRO_REDUX
      --run2d RUN2D         Optional override value for the enviro variable $RUN2D
      --run1d RUN1D         Optional override value for the enviro variable $RUN2D
      --custom CUSTOM       Name of Custom Coadd schema
      --plot, -p            Create spec plots
      --epoch, -e           Run for epoch Coadds
      --lsdr10              Include Legacy Survey DR10 links on HTML
      --allsky              Reformat for Allsky Custom Coadd

.. _spplan:

spplan
^^^^^^
::
 
    usage: spplan [-h] [--skip2d] [--skip1d] [--module MODULE] [--topdir TOPDIR]
                  [--run2d RUN2D] [--lco] [--logfile LOGFILE] [--verbose VERBOSE]
                  [-c] [--release RELEASE] [--remote] [--override_manual]
                  [--mjd [MJD [MJD ...]]] [--mjdstart MJDSTART] [--mjdend MJDEND]
                  [--field [FIELD [FIELD ...]]] [--fieldstart FIELDSTART]
                  [--fieldend FIELDEND] [--legacy] [--plates] [--fps] [--sdssv]
                  [--no_commissioning] [--no_dither] [--matched_flats]
                  [--nomatched_arcs] [--minexp MINEXP] [--single_flat]
                  [--multiple_arc] [--manual_noarc] [--plate_epoch] [--quick]
    
    Produce the spPlan2d and spPlancomb files for the pipeline run
    
    optional arguments:
      -h, --help            show this help message and exit
    
    General:
      General Setup Options
    
      --skip2d              Skip spplan2d
      --skip1d              Skip spplan1d
      --module MODULE       Module file to load for run
      --topdir TOPDIR       Base run2d directory to override module or
                            environmental variable
      --run2d RUN2D         Run2d to override module or environmental variable
      --lco                 Build Run files for LCO
      --logfile LOGFILE     Optional logfile (Including path)
      --verbose VERBOSE     Provide information about nonutlized frames
      -c, --clobber         overwrites previous plan file
      --release RELEASE     sdss_access data release (defaults to sdsswork),
                            required if you do not have proprietary access,
                            otherwise see https://sdss-
                            access.readthedocs.io/en/latest/auth.html#auth
      --remote              allow for remote access to data using sdss-access
      --override_manual     Override/clobber manually edited plan
    
    MJD/Field Filtering:
      MJD/Field Filtering Options
    
      --mjd [MJD [MJD ...]]
                            Use data from these MJDs.
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
      --field [FIELD [FIELD ...]]
                            Use data from these fields.
      --fieldstart FIELDSTART
                            Starting Field
      --fieldend FIELDEND   Ending Field
      --legacy              Include legacy (BOSS/eBOSS) plates
      --plates              Include SDSS-V plates
      --fps                 Include FPS Fields
      --sdssv               Include both SDSS-V Fields & Plates
      --no_commissioning    Exclude SDSS-V FPS Commission Fields
      --no_dither           Exclude Dither fields
    
    RUN2D:
      spPlan2d Setup Options
    
      --matched_flats       Require Flat from a field/plate
      --nomatched_arcs      Allow Arc from another field/plate
      --minexp MINEXP       Min Science Exposures in Plan (default=1)
      --single_flat         Only find the closest flat calibration frame
      --multiple_arc        Find all possible arc calibration frames
      --manual_noarc        if nomatched_arcs is False, builds spplan with
                            unmatched arcs and mark as manual
    
    RUN1D:
      spPlancomb Setup Options
    
      --plate_epoch         Use a variable max epoch length for plate coadd
      --quick               Use the list of new spPlan2d as a filter for fields

.. _spplan_epoch:

spplan_epoch
^^^^^^^^^^^^
::
 
    usage: spplan_epoch [-h] [--module MODULE] [--topdir TOPDIR] [--run2d RUN2D]
                        [--run1d RUN1D] [--mjd MJD] [--mjdstart MJDSTART]
                        [--mjdend MJDEND] [--field FIELD] [--fieldst FIELDSTART]
                        [--fieldend FIELDEND] [--fps] [--sdssv] [--clobber]
                        [--minexp MINEXP] [--lco] [--logfile LOGFILE]
                        [--abandoned] [--started] [--min_epoch_len MIN_EPOCH_LEN]
                        [--release RELEASE] [--remote]
    
    Builds the spPlancombepoch files
    
    optional arguments:
      -h, --help            show this help message and exit
      --module MODULE       Module file to load for run
      --topdir TOPDIR       Override value for the environment variable
                            $BOSS_SPECTRO_REDUX.
      --run2d RUN2D         Override value for the environment variable $RUN2D
      --run1d RUN1D         Override value for the environment variable $RUN1D
      --mjd MJD             Use data from these MJDs.
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
      --field FIELD         Look for the input data files in topdir/fieldid;
                            default to search all subdirectories. Note that this
                            need not be integer-valued, but could be for example
                            '0306_test'.
      --fieldst FIELDSTART  Starting fieldid
      --fieldend FIELDEND   Ending fieldid
      --fps                 Only produce epoch coadds for FPS Fields
                            (Fields>16000)
      --sdssv               Only produce epoch coadds for SDSS-V Fields
                            (Fields>15000)
      --clobber             If set, then over-write conflicting plan files
      --minexp MINEXP       Set minimum number of Science Frames for plan creation
      --lco                 Create Plans for LCO
      --logfile LOGFILE, -l LOGFILE
                            File for logging
      --abandoned           Create plans for abandoned epochs
      --started             Create plans for started epochs (including unfinished)
      --min_epoch_len MIN_EPOCH_LEN
                            minimum length of epoch required to produce plan
      --release RELEASE     sdss_access data release (defaults to sdsswork),
                            required if you do not have proprietary access,
                            otherwise see https://sdss-
                            access.readthedocs.io/en/latest/auth.html#auth
      --remote              allow for remote access to data using sdss-access

.. _spplan_target:

spplan_target
^^^^^^^^^^^^^
::
 
    usage: spplan_target [-h] (--manual | --batch) [--module MODULE] [--name NAME]
                         [--coaddfile COADDFILE] [--topdir TOPDIR] [--run2d RUN2D]
                         [--run1d RUN1D] [--clobber] [--logfile LOGFILE] [--DR]
                         [--cartons [CARTONS [CARTONS ...]]]
                         [--catalogids [CATALOGIDS [CATALOGIDS ...]]]
                         [--program [PROGRAM [PROGRAM ...]]]
                         [--mjd [MJD [MJD ...]]] [--mjdstart MJDSTART]
                         [--mjdend MJDEND] [--coadd_mjdstart COADD_MJDSTART]
                         [--rerun1d] [--use_catid] [--use_firstcarton] [--useDB]
                         [--lco] [--apo]
    
    Build CatalogID Combine Plan
    
    optional arguments:
      -h, --help            show this help message and exit
      --manual              Manaully run a Coadd Schema (from coaddfile if only
                            name is set)
      --batch               Batch run all active Coadd Schema in batch file
                            located {topdir}/{run2d}/{name}
      --module MODULE       Module file to load for run
      --name NAME           Name of Custom Coadd
      --coaddfile COADDFILE
                            File of store Coadding Schema
      --topdir TOPDIR       Override value for the environment variable
                            $BOSS_SPECTRO_REDUX.
      --run2d RUN2D         Override value for the environment variable $RUN2D
      --run1d RUN1D         Override value for the environment variable $RUN1D
      --clobber             If set, then over-write conflicting plan files
      --logfile LOGFILE     File for logging
      --DR                  DR/IPL Batch Coadding
      --cartons [CARTONS [CARTONS ...]]
                            list of cartons
      --catalogids [CATALOGIDS [CATALOGIDS ...]]
                            list of sdss_ids (or catalogids)
      --program [PROGRAM [PROGRAM ...]]
                            list of programs
      --mjd [MJD [MJD ...]]
                            Use data from these MJDs.
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
      --coadd_mjdstart COADD_MJDSTART
                            First Coadd MJD to include
      --rerun1d             Provides flag for coadd to be rerun though 1D analysis
      --use_catid, -u       Uses CatalogID rather then sdss_id
      --use_firstcarton     Use Firstcarton only for carton match (dont look at
                            db)
      --useDB               Use sdss targetdb instead of the Semaphore targeting
                            flag (if not use_firstcarton)
      --lco                 Create Plans for LCO
      --apo                 Create Plans for APO

.. _spplan_trace:

spplan_trace
^^^^^^^^^^^^
::
 
    usage: spplan_trace [-h] [--module MODULE] [--topdir TOPDIR] [--run2d RUN2D]
                        [--lco] [--logfile LOGFILE] [--verbose] [-c]
                        [--release RELEASE] [--remote] [--override_manual]
                        [--mjd [MJD [MJD ...]]] [--mjdstart MJDSTART]
                        [--mjdend MJDEND]
    
    Produces spPlanTrace
    
    optional arguments:
      -h, --help            show this help message and exit
    
    General:
      General Setup Options
    
      --module MODULE       Module file to load for run
      --topdir TOPDIR
      --run2d RUN2D         Run2d to override module or environmental variable
      --lco                 Build Run files for LCO
      --logfile LOGFILE     Optional logfile (Including path)
      --verbose             Provide information about nonutlized frames
      -c, --clobber         overwrites previous plan file
      --release RELEASE     sdss_access data release (defaults to sdsswork),
                            required if you do not have proprietary access,
                            otherwise see https://sdss-
                            access.readthedocs.io/en/latest/auth.html#auth
      --remote              allow for remote access to data using sdss-access
      --override_manual     Override/clobber manually edited plan
    
    MJD/Field Filtering:
      MJD/Field Filtering Options
    
      --mjd [MJD [MJD ...]]
                            Use data from these MJDs.
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD

.. _sxpar.py:

sxpar.py
^^^^^^^^
::
 
    usage: sxpar.py [-h] [-v] fitsfile keyword
    
    Simply parse a fits header
    
    positional arguments:
      fitsfile       The fits file to read
      keyword        Header keyword to parse
    
    optional arguments:
      -h, --help     show this help message and exit
      -v, --verbose  verbose

.. _sxpar_retry.py:

sxpar_retry.py
^^^^^^^^^^^^^^
::
 
    usage: sxpar_retry.py [-h] [-v] fitsfile keyword
    
    Simply parse a fits header, retrying if failed
    
    positional arguments:
      fitsfile       The fits file to read
      keyword        Header keyword to parse
    
    optional arguments:
      -h, --help     show this help message and exit
      -v, --verbose  verbose

.. _update_flags:

update_flags
^^^^^^^^^^^^
::
 
    usage: update_flags [-h] [--run2d RUN2D] [--topdir TOPDIR] [--clobber]
                        [--custom [CUSTOM [CUSTOM ...]]] [--nobackup]
    
    Update SDSSV Targeting flats inn the summary files
    
    optional arguments:
      -h, --help            show this help message and exit
      --run2d RUN2D         idlspec2d Run2d version
      --topdir TOPDIR       idlspec2d Run2d topdir
      --clobber             Clobber spTargeting file
      --custom [CUSTOM [CUSTOM ...]]
                            List of name of custom coadd schema
      --nobackup            Skip backup of existing

.. _uubatchpbs:

uubatchpbs
^^^^^^^^^^
::
 
    usage: uubatchpbs [-h] [--sdssv] [--sdssv_fast] [--sdssv_noshare] [--apo]
                      [--lco] [--bay15] [--merge3d] [--obs [OBS [OBS ...]]]
                      [--topdir TOPDIR] [--run1d RUN1D] [--run2d RUN2D]
                      [--idlutils_1d IDLUTILS_1D] [--no_reject] [--MWM_fluxer]
                      [--map3d {bayestar15,bay15,merge3d}] [--no_healpix]
                      [--noxcsao] [--skip_specprimary] [--no_merge_spall]
                      [--skip2d] [--only1d] [--onestep_coadd] [--fibermap_clobber]
                      [--saveraw] [--debug] [--no_db] [--fast_no_db FAST_NO_DB]
                      [--release RELEASE] [--dr19] [--a2t]
                      [--field [FIELD [FIELD ...]]] [--fieldstart FIELDSTART]
                      [--fieldend FIELDEND] [--mjd [MJD [MJD ...]]]
                      [--mjdstart MJDSTART] [--mjdend MJDEND] [--no_write]
                      [--shared] [--fast] [--mem_per_cpu MEM_PER_CPU]
                      [--walltime WALLTIME] [--nodes NODES] [--ppn PPN]
                      [--nosubmit] [--clobber] [--epoch] [--custom CUSTOM]
                      [--allsky] [--coadd_only] [--1dpost] [--email]
    
    Build idlspec2d redux and submit to slurm. Without access to the SDSS Slurm
    package, it prints the commands for manual execution
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Short cuts:
      --sdssv               --mwm --no_merge_spall --no_reject --shared
      --sdssv_fast          --sdssv --fast --shared
      --sdssv_noshare       --sdssv (without --shared)
      --apo                 Run apo only
      --lco                 Run lco only
      --bay15               Set map3d to bayestar15 model
      --merge3d             Set map3d to best 3d model
    
    idlspec2d Run options:
      --obs [OBS [OBS ...]]
                            Observatory {apo,lco}
      --topdir TOPDIR       Optional override value for the environment variable
                            $BOSS_SPECTRO_REDUX
      --run1d RUN1D         Optional override value for the enviro variable $RUN1D
      --run2d RUN2D         Optional override value for the enviro variable $RUN2D
      --idlutils_1d IDLUTILS_1D
                            idlutils override version of spec1d
      --no_reject           Deactivate Rejection in Coadd
      --MWM_fluxer, --mwm
      --map3d {bayestar15,bay15,merge3d}
                            Name of 3d dustmap to use with MWM_fluxer
                            (default=None)
      --no_healpix, --nohp  Turn off copy to healpix
      --noxcsao             Skip pyXCSAO
      --skip_specprimary    Skip Calculation of Specprimary
      --no_merge_spall      Skip building full SpAll File
      --skip2d              Skip spreduce2d
      --only1d              run spec1d step only (eg. spreduce1d_empca, XCSAO)
      --onestep_coadd       Use legacy one step version of coadd
      --fibermap_clobber    Clobber spfibermap fits file
      --saveraw             Save sdssproc outputs
      --debug               Save extraction debug files
      --no_db               skip Database operations
      --fast_no_db FAST_NO_DB
                            When using --no_db, streamlines process and only gets
                            parallax from MOS target files
      --release RELEASE     sdss_access data release (defaults to sdsswork),
                            required if you do not have proprietary access,
                            otherwise see https://sdss-
                            access.readthedocs.io/en/latest/auth.html#auth
      --dr19                Limit targeting flags to DR19 cartons
      --a2t                 Force Use of Arc2Trace
    
    Select Fields:
      --field [FIELD [FIELD ...]], -f [FIELD [FIELD ...]]
                            Plate/Field numbers to reduce default="*"
      --fieldstart FIELDSTART
                            Starting Field/Plate number
      --fieldend FIELDEND   End Field/Plate number
    
    Select MJDs:
      --mjd [MJD [MJD ...]], -m [MJD [MJD ...]]
                            MJD dates to reduce; default="*"
      --mjdstart MJDSTART   Starting MJD
      --mjdend MJDEND       Ending MJD
    
    Slurm Options:
      --no_write            skip writing and submitting job
      --shared              Node sharing
      --fast                Use SDSS fast queue
      --mem_per_cpu MEM_PER_CPU
                            Memory allocated per CPU
      --walltime WALLTIME   Wall time in hours
      --nodes NODES         Number of Nodes
      --ppn PPN             Number of processors per node
      --nosubmit            Build, but not submit redux files
      --clobber             Clobber redux
    
    Custom Coadd Options:
      --epoch               Epoch Coadds
      --custom CUSTOM       Name of custom Coadd Schema
      --allsky              All Sky Coadds
      --coadd_only          Run spspec_target_merge only
      --1dpost              Run 1d analysis and post processing only
    
    Email outputs:
      --email               Email log using $DAILY_DIR/etc/emails

.. _uurundaily:

uurundaily
^^^^^^^^^^
::
 
    usage: uurundaily [-h] [--module MODULE] [--apo] [--lco]
                      [--mjd [MJD [MJD ...]]] [--range_mjd RANGE_MJD]
                      [--no_dither] [--epoch] [--no_merge3d] [--summary]
                      [--no_traceflat] [--no_prep] [--no_fibermap]
                      [--skip_plan [{pipe,trace,True,all}]]
                      [--clobber [{spPlans,fibermap,trace,True,all} [{spPlans,fibermap,trace,True,all} ...]]]
                      [--saveraw] [--debug] [--fast] [--nosubmit] [--noslurm]
                      [--batch] [--nodb] [--monitor] [--allemail] [--pause PAUSE]
                      [--walltime WALLTIME] [--mem_per_cpu MEM_PER_CPU]
    
    Process the BOSS data for a single MJD end-to-end (including plan files)
    
    optional arguments:
      -h, --help            show this help message and exit
      --module MODULE       Module for daily run
      --no_merge3d          Skip using prototype 3D Dustmap (in merge mode)
    
    Field-MJD Selection:
      Arguments to control the Field-MJD Selection to run
    
      --apo                 Run for APO Only
      --lco                 Run for LCO Only
      --mjd [MJD [MJD ...]]
                            Manually run for a single/list of mjd (does not update
                            nextmjd.par)
      --range_mjd RANGE_MJD
                            Manually run for a range of mjds (does not update
                            nextmjd.par)
      --no_dither           Skip Dither Engineering Fields
      --epoch               Run Epoch Coadds
    
    Pipeline Steps:
      Arguments to control which steps of the full pipeline are run
    
      --summary             Build Summary Files
      --no_traceflat        Skip Building and using TraceFlats
      --no_prep             Skip building TraceFlats and spfibermaps before
                            pipeline run
      --no_fibermap         Skip Pre-Run of readfibermap
      --skip_plan [{pipe,trace,True,all}]
                            Skip the given plan {pipe,trace,all (flagging
                            --skip_plan with name will default to all)}
      --clobber [{spPlans,fibermap,trace,True,all} [{spPlans,fibermap,trace,True,all} ...]]
                            Clobber uubatchpbs + a combo of spPlan, fibermap, and
                            TraceFlat run {fibermap,trace, all (flagging --clobber
                            with name will default to all)}
    
    Debug:
      Arguments to saving of optional debugging files
    
      --saveraw             save sdssproc outputs
      --debug               save extraction debug files
    
    Pipeline Options:
      Arguments to set the misc pipeline options
    
      --fast                turn on fast user for slurm
      --nosubmit            Skip submitting uubatch job (ideal for allowing
                            editting of plans)
      --noslurm             Skip creating uubatch job
      --batch               run for multiple mjds in a single batch
      --nodb                skip Database operations
      --monitor             Monitors pipeline status
      --allemail            Email intermediate log using all emails in
                            $DAILY_DIR/etc/emails (defaults to first email only)
      --pause PAUSE         Pause time (s) in status updates
      --walltime WALLTIME   Wall time in hours
      --mem_per_cpu MEM_PER_CPU
                            Memory allocated per CPU

IDL Command Usage
-----------------

.. contents::
    :depth: 3
    :local:
    :class: this-will-duplicate-information-and-it-is-still-useful-here
    :backlinks: none


.. _spreduce2d.pro:

spreduce2d.pro
^^^^^^^^^^^^^^
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

.. _rm_combine_script.pro:

rm_combine_script.pro
^^^^^^^^^^^^^^^^^^^^^
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

.. _spreduce1d_empca.pro:

spreduce1d_empca.pro
^^^^^^^^^^^^^^^^^^^^
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

.. _spcalib_qa.pro:

spcalib_qa.pro
^^^^^^^^^^^^^^
::
 
    ; NAME:
    ;   spcalib_qa
    ;
    ; PURPOSE:
    ;   Compare photometric accuracy of standards
    ;
    ; CALLING SEQUENCE:
    ;   SpCalib_QA, [run2d=, fieldid=, mjd=]
    ;
    ; INPUTS:
    ;
    ; OPTIONAL INPUTS:
    ;   field       - field to include
    ;   mjd         - MJD to include
    ;   run2d       - RUN2D version of reduction
    ;
    ; OUTPUTS:
    ;
    ; OPTIONAL OUTPUTS:
    ;
    ; COMMENTS:
    ;   Depends on the spAll files (either full run2d version or field-mjd version)
    ;
    ; EXAMPLES:
    ;
    ; BUGS:
    ;
    ; DATA FILES:
    ;
    ; Function Called:
    ;   mpfitfun
    ;   field_to_string
    ;   djs_filepath
    ;   mrdfits
    ;   sdss_flagval
    ;
    ; External PROCEDURES CALLED:
    ;   plot
    ;   XYOUTS
    ;   cpbackup
    ;
    ; Internal PROCEDURES CALLED:
    ;   std_hist
    ;
    ; REVISION HISTORY:
    ;   21-June-2022  Written by S. Morrison (UIUC)

.. _spspec_target_merge.pro:

spspec_target_merge.pro
^^^^^^^^^^^^^^^^^^^^^^^
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
